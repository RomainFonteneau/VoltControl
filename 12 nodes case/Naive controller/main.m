clear all

%% Parameters definition
seed = 1234;
nt = 50;  % Number of simulation steps
N  = 3;   % MPC prediction horizon (steps)

% MPC cost weights: min sum_i [ ||V_load - V_ref||^2 + alpha*||Q_gen||^2 ]
% Voltage tracking term has coefficient 1; alpha penalises reactive power deviation.
alpha = 0.001;
u_max = 0.2;    % Maximum generator voltage increment per step (p.u.)

% Gaussian noise standard deviations
sigma_rate = 0;  % Rate used to compute std deviations from nominal values; 0 = no noise

% Jacobian recalculation period: 0 = every step, n = every n steps
t_jac = 0;

%----------Capacitor (lowest priority actuator)------------
add_cap_bool   = true;
Q_lim          = 0.2;   % Soft reactive power limit (p.u.)
t_cap          = 5;     % Minimum steps between two switches on the same bus (cooldown)
var_v_cap_last = 0.01;  % Max allowed voltage variation (p.u.) over the past t_cap_last steps
t_cap_last     = 3;
var_v_cap_next = 0.01;  % Max allowed reference variation (p.u.) over the next t_cap_next steps
t_cap_next     = 3;
cap            = 5;     % Capacitor/coil size (Mvar)

%---------OLTC (naive: direct reference voltage)------------
% The OLTC is commanded by passing the voltage reference of each secondary bus
% directly to run_pf_oltc_step. The tap moves one step per call whenever
% |V_actual - V_ref| > tol, without any enumeration, cost function or cooldown.
% run_pf_oltc_step parameters — must match the function internals
step_size = 0.00625;  % Tap step size (p.u.)
tol       = 0.005;    % Voltage deadband

%% Preparation
rng(seed);

mpc = case9xx_Bsh();
mpopt = mpoption('verbose',0,'out.all',0,'pf.enforce_q_lims',1);
define_constants;
baseMVA = mpc.baseMVA;
fullJac = true;

% Branch indices of the two 400/63 kV OLTC transformers (bus 10->5 and bus 11->6)
transformer_indices = [16; 17];
n_oltc = length(transformer_indices);

% Bus index sets
pv_idx    = find(mpc.bus(:,BUS_TYPE)==PV);
slack_idx = find(mpc.bus(:,BUS_TYPE)==REF);
load_idx  = find(mpc.bus(:,BUS_TYPE)==PQ);
target_idx= find(mpc.bus(:,PD)~=0);

n_load   = size(load_idx,  1);
n_gen    = size(pv_idx,    1);
n_slack  = size(slack_idx, 1);
n_bus    = n_load + n_gen + n_slack;
n_target = size(target_idx,1);
n_state  = n_load + n_gen + n_gen;

% gen_number(i) = row index in mpc.gen of the generator at bus i
gen_number = zeros(n_bus, 1);
for i = 1:n_gen+n_slack
    gen_number(mpc.gen(i,1)) = i;
end

% Hard limits
Q_gen_max  = mpc.gen(pv_idx, QMAX) / baseMVA;
Q_gen_min  = mpc.gen(pv_idx, QMIN) / baseMVA;
V_gen_max  = mpc.bus(pv_idx,   VMAX);
V_gen_min  = mpc.bus(pv_idx,   VMIN);
V_load_max = mpc.bus(load_idx, VMAX);
V_load_min = mpc.bus(load_idx, VMIN);

% Nominal values for noise injection
Pd_nom = mpc.bus(target_idx, PD);
Qd_nom = mpc.bus(target_idx, QD);
Pg_nom = mpc.gen(gen_number(pv_idx), PG);

sigma_load_P = Pd_nom * sigma_rate;
sigma_load_Q = Qd_nom * sigma_rate;
sigma_gen_P  = Pg_nom * sigma_rate;

cap_availability = zeros(n_bus, 1);

% Secondary (LV) buses and their position within target_idx.
% oltc_in_target(j): row index in Yc_V of transformer j's secondary bus voltage reference.
reg_buses = mpc.branch(transformer_indices, T_BUS);
oltc_in_target = zeros(n_oltc, 1);
for j = 1:n_oltc
    idx = find(target_idx == reg_buses(j));
    if ~isempty(idx)
        oltc_in_target(j) = idx;
    end
end

%% Objective
V_target=[1.03;1.02;0.98];% Fixed reference for buses 5, 6 and 8
Yc_V = repmat(V_target, 1, nt + max(N, t_cap_next) + 1);
Yc_Q = zeros(n_gen, nt + max(N, t_cap_next) + 1);
Yc   = repmat([V_target; zeros(n_gen,1)], nt + max(N, t_cap_next) + 1, 1);

%% Constant MPC matrices (computed once)
[GAMMA_cst,f1_cst,f2_cst,b1,b2,Hu,Hx,PSI,C_tilde,A,C,J1] = cal_mat_cst(mpc,N,u_max,alpha);

%% Initial power flow
results = runpf(mpc, mpopt);

x = [results.bus(load_idx, VM); results.gen(gen_number(pv_idx), QG)/baseMVA; results.bus(pv_idx, VM)];
JAC = makeJac(results, fullJac);

% Logging
X_log     = zeros(n_state, nt+1);   X_log(:, 1) = x;
u_log     = zeros(n_gen,   nt);
J_log     = zeros(4,       nt+1);
tap_log   = zeros(n_oltc,  nt+2);   tap_log(:, 1:2) = repmat(mpc.branch(transformer_indices, TAP), 1, 2);
cap_log   = zeros(n_bus, 1);
Pd_log    = zeros(n_target, nt);
Qd_log    = zeros(n_target, nt);
Pg_log    = zeros(n_gen,    nt);
V_sec_log = zeros(n_oltc,   nt);  % target voltages given to run_pf_oltc_step

[J_Q, J_V, J_VQ] = Cal_criterion(results, C*x, Yc_V, Yc_Q, 1, alpha);
J_log(:, 1) = [0; J_Q; J_V; J_VQ];

%% Control loop
for k = 0:nt-1

    % State-dependent MPC matrices from the current Jacobian
    [Aineq, H, f1, f2, B] = Cal_mat_var(mpc, N, GAMMA_cst, Hu, Hx, PSI, C_tilde, f1_cst, f2_cst, JAC);

    % No tap state correction: the naive controller does not anticipate the in-flight tap.
    % The MPC is solved directly from the measured state x.

    % --- Priority 1: MPC (generator voltage) ---
    Yc_k = Yc(k*(n_target+n_gen)+1 : (k+N)*(n_target+n_gen), 1);
    [U, J_qp, exitflag] = quadprog(H, f1*x - f2*Yc_k, Aineq, b1 - b2*x, ...
        [], [], [], [], [], optimset('Display', 'off'));
    if exitflag <= 0
        warning('QP failed at step %d — applying zero input', k);
        U = zeros(n_gen*N, 1);
    end
    u_k  = U(1:n_gen);
    u_k1 = U(n_gen+1 : 2*n_gen);
    u_log(:, k+1) = u_k;

    e_free = J1*x - Yc_k;
    JMPC   = J_qp + e_free' * PSI * e_free;

    % State predictions for capacitor logic
    x_estimated_k1 = A * x + B * u_k;
    x_estimated_k2 = A * x_estimated_k1 + B * u_k1;

    % --- Naive OLTC: pass the voltage reference of each secondary bus directly ---
    % run_pf_oltc_step will move the tap one step whenever |V_actual - target_v| > tol.
    target_v = Yc_V(oltc_in_target, k+1);
    V_sec_log(:, k+1) = target_v;

    % --- Priority 2: Capacitor/coil decision ---
    use_add_cap = ~(all(x_estimated_k2(n_load+1:n_load+n_gen) < Q_lim) && ...
                    all(x_estimated_k2(n_load+1:n_load+n_gen) > -Q_lim));
    if k < t_cap_last || ~use_add_cap || ~add_cap_bool
        cap_log          = [cap_log, cap_log(:, end)];
        cap_availability = max(cap_availability - ones(n_bus,1), zeros(n_bus,1));
    else
        [mpc, cap_availability, cap_log] = add_cap(mpc, JAC, x_estimated_k1, x_estimated_k2, ...
            cap_availability, cap_log, t_cap, var_v_cap_last, t_cap_last, var_v_cap_next, ...
            t_cap_next, cap, k, X_log(:, 1:k+1), Yc_V, Q_lim);
    end

    % Apply MPC command: increment generator voltage setpoints, clamped to physical limits
    mpc.gen(gen_number(pv_idx), VG) = min(max(mpc.gen(gen_number(pv_idx), VG) + u_k, V_gen_min), V_gen_max);

    % Inject Gaussian noise on loads and non-slack generation
    mpc.bus(target_idx, PD) = Pd_nom + sigma_load_P .* randn(length(target_idx), 1);
    mpc.bus(target_idx, QD) = Qd_nom + sigma_load_Q .* randn(length(target_idx), 1);
    mpc.gen(gen_number(pv_idx), PG) = Pg_nom + sigma_gen_P .* randn(n_gen, 1);

    % Run power flow and advance each OLTC tap one step toward target_v
    [mpc, results, final_taps] = run_pf_oltc_step(mpc, target_v, transformer_indices);

    x = [results.bus(load_idx, VM); results.gen(gen_number(pv_idx), QG)/baseMVA; results.bus(pv_idx, VM)];

    X_log(:, k+2)   = x;
    tap_log(:, k+3) = final_taps;
    Pd_log(:, k+1)  = results.bus(target_idx, PD);
    Qd_log(:, k+1)  = results.bus(target_idx, QD);
    Pg_log(:, k+1)  = results.gen(gen_number(pv_idx), PG);

    [J_Q, J_V, J_VQ] = Cal_criterion(results, C*x, Yc_V, Yc_Q, k+1, alpha);
    J_log(:, k+2) = [JMPC; J_Q; J_V; J_VQ];

    % Update Jacobian (capacitor included; tap not yet applied)
    if t_jac == 0 || mod(k+1, t_jac) == 0
        JAC = makeJac(results, fullJac);
    end
end

%% DISPLAY

% %% Active and reactive powers
% figure;
% 
% % 1. Active load power
% subplot(3,1,1);
% hold on; grid on;
% hX = stairs(1:nt, Pd_log');
% set(hX, 'LineWidth', 1.5);
% hN = plot([1 nt], repmat(Pd_nom', 2, 1), '--', 'LineWidth', 1);
% for i = 1:n_target; hN(i).Color = hX(i).Color; end
% legend(hX(:)', compose("Bus %d", target_idx'), 'Location', 'northeast');
% title('Active load power');
% xlabel('k'); ylabel('P_{load} (MW)');
% 
% % 2. Reactive load power
% subplot(3,1,2);
% hold on; grid on;
% hX = stairs(1:nt, Qd_log');
% set(hX, 'LineWidth', 1.5);
% hN = plot([1 nt], repmat(Qd_nom', 2, 1), '--', 'LineWidth', 1);
% for i = 1:n_target; hN(i).Color = hX(i).Color; end
% legend(hX(:)', compose("Bus %d", target_idx'), 'Location', 'northeast');
% title('Reactive load power');
% xlabel('k'); ylabel('Q_{load} (Mvar)');
% 
% % 3. Active generation power (non-slack)
% subplot(3,1,3);
% hold on; grid on;
% hX = stairs(1:nt, Pg_log');
% set(hX, 'LineWidth', 1.5);
% hN = plot([1 nt], repmat(Pg_nom', 2, 1), '--', 'LineWidth', 1);
% for i = 1:n_gen; hN(i).Color = hX(i).Color; end
% legend(hX(:)', compose("Gen %d", pv_idx'), 'Location', 'northeast');
% title('Active generation power (non-slack)');
% xlabel('k'); ylabel('P_{gen} (MW)');


%% Summary
figure;

% 1. Voltage in targeted loads
subplot(3,3,1);
zoom on; hold on; grid on;
idx = find(ismember(load_idx, target_idx));
hX = stairs(0:nt, X_log(idx,:)');
set(hX, 'LineWidth', 1.5)
hY = stairs(1:nt, Yc_V(1:n_target,1:nt)', 'x');
for i = 1:n_target
    hY(i).Color = hX(i).Color;
    hY(i).LineStyle = 'none';
    hY(i).Marker = 'x';
    hY(i).LineWidth = 0.8;
end
hZ = stairs(1:nt,[repmat(V_load_max(idx),1,nt);repmat(V_load_min(idx),1,nt)]');
for i = 1:n_target
    hZ(i).Color = hX(i).Color;
    hZ(i+n_target).Color = hX(i).Color;
    hZ(i).LineStyle = 'none';
    hZ(i+n_target).LineStyle = 'none';
    hZ(i).Marker = 'diamond';
    hZ(i+n_target).Marker = 'diamond';
    hZ(i).LineWidth = 0.8;
    hZ(i+n_target).LineWidth = 0.8;
end
legendHandles = [hX(:)',hY(1),hZ(1)];
legendStrings = [compose("Load %d", target_idx'),"Target","Constraints"];
legend(legendHandles, legendStrings, 'Location', 'northeast');
title('Voltage in targeted loads');
xlabel('k'); ylabel('Voltage (p.u)')

% 2. Reactive power in generators
subplot(3,3,2);
zoom on; hold on; grid on;
hX = stairs(0:nt,X_log(n_load+1:n_load+n_gen,:)');
set(hX, 'LineWidth', 1.5)
hZ = stairs(0:nt,[repmat(Q_gen_max,1,nt+1);repmat(Q_gen_min,1,nt+1)]');
for i = 1:n_gen
    hZ(i).Color = hX(i).Color;
    hZ(i+n_gen).Color = hX(i).Color;
    hZ(i).LineStyle = 'none';
    hZ(i+n_gen).LineStyle = 'none';
    hZ(i).Marker = 'diamond';
    hZ(i+n_gen).Marker = 'diamond';
    hZ(i).LineWidth = 0.8;
    hZ(i+n_gen).LineWidth = 0.8;
end
if add_cap_bool
    hY = stairs(0:nt,[-ones(n_gen,nt+1)*Q_lim;ones(n_gen,nt+1)*Q_lim]');
    for i = 1:n_gen
        hY(i).Color = hX(i).Color;
        hY(i+n_gen).Color = hX(i).Color;
        hY(i).LineStyle = 'none';
        hY(i+n_gen).LineStyle = 'none';
        hY(i).Marker = 'x';
        hY(i+n_gen).Marker = 'x';
        hY(i).LineWidth = 2;
        hY(i+n_gen).LineWidth = 2;
    end
    legendHandles = [hX(:)',hY(1),hZ(1)];
    legendStrings = [compose("Gen %d", pv_idx'),"Targeted limits","Constraints"];
else
    legendHandles = [hX(:)',hZ(1)];
    legendStrings = [compose("Gen %d", pv_idx'),"Constraints"];
end
legend(legendHandles, legendStrings, 'Location', 'northeast');
title('Reactive power in generators');
xlabel('k'); ylabel('Reactive power (p.u)')

% 3. Voltage in generators
subplot(3,3,3);
zoom on; hold on; grid on;
hX = stairs(0:nt,X_log(n_load+n_gen+1:end,:)');
set(hX, 'LineWidth', 1.5)
hZ = stairs(0:nt,[repmat(V_gen_max,1,nt+1);repmat(V_gen_min,1,nt+1)]');
for i = 1:n_gen
    hZ(i).Color = hX(i).Color;
    hZ(i+n_gen).Color = hX(i).Color;
    hZ(i).LineStyle = 'none';
    hZ(i+n_gen).LineStyle = 'none';
    hZ(i).Marker = 'diamond';
    hZ(i+n_gen).Marker = 'diamond';
    hZ(i).LineWidth = 0.8;
    hZ(i+n_gen).LineWidth = 0.8;
end
legendHandles = [hX(:)',hZ(1)];
legendStrings = [compose("Gen %d", pv_idx'),"Constraints"];
legend(legendHandles, legendStrings, 'Location', 'northeast');
title('Voltage in generators');
xlabel('k'); ylabel('Voltage (p.u)')

% 4. Voltage in loads or capacitors
subplot(3,3,4)
zoom on; hold on; grid on;
if ~add_cap_bool
    hX = stairs(0:nt, X_log(1:n_load,:)');
    set(hX, 'LineWidth', 1.5)
    hZ = stairs(1:nt,[repmat(V_load_max,1,nt);repmat(V_load_min,1,nt)]');
    for i = 1:n_load
        hZ(i).Color = hX(i).Color;
        hZ(i+n_load).Color = hX(i).Color;
        hZ(i).LineStyle = 'none';
        hZ(i+n_load).LineStyle = 'none';
        hZ(i).Marker = 'diamond';
        hZ(i+n_load).Marker = 'diamond';
        hZ(i).LineWidth = 0.8;
        hZ(i+n_load).LineWidth = 0.8;
    end
    legendHandles = [hX(:)',hZ(1)];
    legendStrings = [compose("Load %d", load_idx'),"Constraints"];
    legend(legendHandles, legendStrings, 'Location', 'northeast');
    title('Voltage in loads');
    xlabel('k'); ylabel('Voltage (p.u)')
else
    hX = stairs(0:nt,cap_log(target_idx,:)');
    set(hX, 'LineWidth', 1.5)
    legendStrings = [compose("Load %d", target_idx')];
    legend(legendStrings{:});
    title('Capacitors');
    xlabel('k'); ylabel('Capacitors')
end

% 5. Case
subplot(3,3,5)
img = imread('case.png');
imshow(img);
title('Case considered');

% 6. Voltage input (generators)
subplot(3,3,6)
zoom on; hold on; grid on;
hX = stairs(0:(nt-1),u_log');
set(hX, 'LineWidth', 1.5)
hZ = stairs(0:(nt-1),[-u_max*ones(N,nt);u_max*ones(N,nt)]');
for i = 1:n_gen
    hZ(i).Color = hX(i).Color;
    hZ(i+n_gen).Color = hX(i).Color;
    hZ(i).LineStyle = 'none';
    hZ(i+n_gen).LineStyle = 'none';
    hZ(i).Marker = 'diamond';
    hZ(i+n_gen).Marker = 'diamond';
    hZ(i).LineWidth = 0.8;
    hZ(i+n_gen).LineWidth = 0.8;
end
legendHandles = [hX(:)',hZ(1)];
legendStrings = [compose("Gen %d", pv_idx'),"Constraints"];
legend(legendHandles, legendStrings, 'Location', 'northeast');
title('Voltage input');
xlabel('k'); ylabel('Voltage input (p.u)')

% 7. Tap
subplot(3,3,7);
zoom on; hold on; grid on;
hX = stairs(0:nt+1, tap_log');
set(hX, 'LineWidth', 1.5)
legendHandles = hX(:)';
legendStrings = [compose("transformer %d-%d (branch %d)", [mpc.branch(transformer_indices,1:2),transformer_indices])];
legend(legendHandles, legendStrings, 'Location', 'northeast');
title('Tap')
xlabel('k'); ylabel('Tap')

% 8. Criterion on V
subplot(3,3,8)
zoom on; hold on; grid on;
hX = stairs(0:nt, J_log(3,:)');
set(hX, 'LineWidth', 1.5);
x_last = nt;
y_last = J_log(3,nt);
plot(x_last, y_last, 'ro','MarkerFaceColor','r');
text(x_last, y_last, sprintf(' %.3e', y_last), ...
    'VerticalAlignment','bottom','HorizontalAlignment','left');
title('Criterion on V');
xlabel('k'); ylabel('sum (V-V_{target})^2')

% 9. Capacitors
if add_cap_bool
    subplot(3,3,9);
    zoom on; hold on; grid on;
    hX = stairs(0:nt,cap_log(target_idx,:)');
    set(hX, 'LineWidth', 1.5)
    legendStrings = [compose("Load %d", target_idx')];
    legend(legendStrings{:});
    title('Capacitors');
    xlabel('k'); ylabel('Capacitors')
end

% %% DISPLAY REPORT — Axis 1: controller comparison
% 
% fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures', 'ComparaisonVersion');
% if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
% prefix = 'naive';
% 
% idx_tgt = find(ismember(load_idx, target_idx));
% reg_buses_rpt = mpc.branch(transformer_indices, T_BUS);
% [~, reg_pos_rpt] = ismember(reg_buses_rpt, load_idx);
% 
% % ---- Figure 1: Target voltages / Reactive power ----
% fig1 = figure('Units','centimeters','Position',[2 2 16 12]);
% 
% subplot(2,1,1); hold on; grid on;
% hV = stairs(0:nt, X_log(idx_tgt,:)');
% set(hV, 'LineWidth', 1.5);
% for i = 1:n_target
%     plot([0 nt], Yc_V(i,1)*[1 1], '--', 'Color', hV(i).Color, 'LineWidth', 1);
%     plot([0 nt], V_load_max(idx_tgt(i))*[1 1], ':', 'Color', [.6 .6 .6], 'LineWidth', .8, 'HandleVisibility','off');
%     plot([0 nt], V_load_min(idx_tgt(i))*[1 1], ':', 'Color', [.6 .6 .6], 'LineWidth', .8, 'HandleVisibility','off');
% end
% hRef = plot(nan, nan, 'k--', 'LineWidth', 1);
% hCst = plot(nan, nan, ':', 'Color', [.6 .6 .6], 'LineWidth', .8);
% legend([hV(:)', hRef, hCst], [compose("Bus %d", target_idx'), "Reference", "Limits"], 'Location','northeast');
% title('Target bus voltages'); xlabel('k'); ylabel('Voltage (p.u.)');
% 
% subplot(2,1,2); hold on; grid on;
% hQ = stairs(0:nt, X_log(n_load+1:n_load+n_gen,:)');
% set(hQ, 'LineWidth', 1.5);
% plot([0 nt], [ Q_lim  Q_lim], 'k--', 'LineWidth', 1, 'HandleVisibility','off');
% plot([0 nt], [-Q_lim -Q_lim], 'k--', 'LineWidth', 1, 'HandleVisibility','off');
% for i = 1:n_gen
%     plot([0 nt], Q_gen_max(i)*[1 1], ':', 'Color',[.6 .6 .6], 'LineWidth',.8, 'HandleVisibility','off');
%     plot([0 nt], Q_gen_min(i)*[1 1], ':', 'Color',[.6 .6 .6], 'LineWidth',.8, 'HandleVisibility','off');
% end
% legend(hQ(:)', compose("Gen %d", pv_idx'), 'Location','northeast');
% title('Generator reactive power'); xlabel('k'); ylabel('Q_{gen} (p.u.)');
% 
% exportgraphics(fig1, fullfile(fig_dir, [prefix '_volt_reactive.pdf']), 'ContentType','vector');
% 
% % ---- Figure 2: Reactive power / Capacitors ----
% fig2 = figure('Units','centimeters','Position',[2 2 16 12]);
% 
% subplot(2,1,1); hold on; grid on;
% hQ = stairs(0:nt, X_log(n_load+1:n_load+n_gen,:)');
% set(hQ, 'LineWidth', 1.5);
% plot([0 nt], [ Q_lim  Q_lim], 'k--', 'LineWidth', 1, 'HandleVisibility','off');
% plot([0 nt], [-Q_lim -Q_lim], 'k--', 'LineWidth', 1, 'HandleVisibility','off');
% for i = 1:n_gen
%     plot([0 nt], Q_gen_max(i)*[1 1], ':', 'Color',[.6 .6 .6], 'LineWidth',.8, 'HandleVisibility','off');
%     plot([0 nt], Q_gen_min(i)*[1 1], ':', 'Color',[.6 .6 .6], 'LineWidth',.8, 'HandleVisibility','off');
% end
% legend(hQ(:)', compose("Gen %d", pv_idx'), 'Location','northeast');
% title('Generator reactive power'); xlabel('k'); ylabel('Q_{gen} (p.u.)');
% 
% subplot(2,1,2); hold on; grid on;
% hC = stairs(0:nt, cap_log(target_idx,:)');
% set(hC, 'LineWidth', 1.5);
% legend(hC(:)', compose("Bus %d", target_idx'), 'Location','northeast');
% title('Capacitors / reactors'); xlabel('k'); ylabel('Number of units');
% 
% exportgraphics(fig2, fullfile(fig_dir, [prefix '_reactive_cap.pdf']), 'ContentType','vector');
% 
% % ---- Figure 3: Tap / Target voltage given to OLTC ----
% fig3 = figure('Units','centimeters','Position',[2 2 16 12]);
% 
% subplot(2,1,1); hold on; grid on;
% hT = stairs(0:nt+1, tap_log');
% set(hT, 'LineWidth', 1.5);
% legend(hT(:)', compose("Transfo %d to %d", mpc.branch(transformer_indices,1:2)), 'Location','northeast');
% title('OLTC tap position'); xlabel('k'); ylabel('Tap ratio (p.u.)');
% 
% subplot(2,1,2); hold on; grid on;
% hAct = stairs(0:nt, X_log(reg_pos_rpt,:)');
% set(hAct, 'LineWidth', 1.5);
% hTgt = stairs(1:nt, V_sec_log', '--');
% for i = 1:n_oltc
%     hTgt(i).Color = hAct(i).Color;
%     hTgt(i).LineWidth = 1;
% end
% hDum1 = plot(nan, nan, 'k-',  'LineWidth', 1.5);
% hDum2 = plot(nan, nan, 'k--', 'LineWidth', 1);
% legend([hDum1, hDum2], {'Actual secondary voltage', 'Target given to OLTC'}, 'Location','northeast');
% title('Secondary bus voltage and OLTC target'); xlabel('k'); ylabel('Voltage (p.u.)');
% 
% exportgraphics(fig3, fullfile(fig_dir, [prefix '_tap_voltcall.pdf']), 'ContentType','vector');
% 
% % ---- Figure 4: Criterion J_V ----
% fig4 = figure('Units','centimeters','Position',[2 2 16 6]);
% hold on; grid on;
% stairs(0:nt, J_log(3,:)', 'LineWidth', 1.5);
% plot(nt, J_log(3,end), 'ro', 'MarkerFaceColor','r');
% text(nt, J_log(3,end), sprintf('  %.3e', J_log(3,end)), ...
%     'VerticalAlignment','bottom', 'HorizontalAlignment','left');
% title('Voltage criterion J_V'); xlabel('k'); ylabel('\Sigma(V - V_{ref})^2');
% 
% exportgraphics(fig4, fullfile(fig_dir, [prefix '_criterion.pdf']), 'ContentType','vector');
