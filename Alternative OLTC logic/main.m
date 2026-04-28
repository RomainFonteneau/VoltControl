clear all

%% Parameters definition
seed = 1234;
nt = 50;  % Number of simulation steps
N  = 3;   % MPC prediction horizon (steps)

% MPC cost weights: min sum_i [ ||V_load - V_ref||^2 + alpha*||Q_gen||^2 ]
% Voltage tracking term has coefficient 1; alpha penalises reactive power deviation.
% Hessian is scaled by 1/alpha inside cal_mat_cst to avoid numerical issues.
alpha = 0.001;  % Weight on Q_gen
u_max = 0.2;    % Maximum generator voltage increment per step: |u(k)| <= u_max (p.u.)

% Gaussian noise standard deviations
sigma_rate = 0;  % Rate used to compute std deviations from nominal values; 0 = no noise

% Jacobian recalculation period: 0 = every step, n = every n steps
% The Jacobian is shared by all sensitivity functions (Cal_CvCq, Cal_Cv_cap, Cal_Cv_tap)
% It is recomputed after capacitor switching but before the power flow,
% so it reflects the current network topology (capacitors) but not yet the delayed tap.
t_jac = 0;

%----------Capacitor (lowest priority actuator)------------
add_cap_bool   = true;  % Set to false to disable capacitor/coil switching entirely
Q_lim          = 0.2;   % Soft reactive power limit (p.u.): capacitors/coils act to keep Q_gen in [-Q_lim, Q_lim]
                        % Stricter than the hard limits QMAX/QMIN enforced by the QP
t_cap          = 5;     % Minimum number of steps between two switches on the same bus (cooldown)
var_v_cap_last = 0.01;  % Max allowed voltage variation (p.u.) at the target bus over the past t_cap_last steps
t_cap_last     = 3;     % Number of past steps used for the voltage stability check
var_v_cap_next = 0.01;  % Max allowed reference variation (p.u.) over the next t_cap_next steps
t_cap_next     = 3;     % Number of future steps used for the reference stability check
cap            = 5;     % Capacitor/coil size (Mvar); adding a coil is equivalent to cap = -cap

%---------OLTC (intermediate priority actuator)------------
% OLTC enumeration criterion: J = ||V_load - V_ref||^2 + alpha_tap*||Q_gen_out||^2
% Independent from the MPC weights; evaluated at k+2 (two-step tap delay).
alpha_tap = 0.001;  % Weight on reactive power in the OLTC criterion

step_size = 0.00625;  % Tap step size per action (p.u.); must match run_pf_oltc_step
tol       = 0.005;    % Voltage deadband used to encode tap direction; must match run_pf_oltc_step
tap_min   = 0.9;      % Physical tap lower bound; must match run_pf_oltc_step
tap_max   = 1.1;      % Physical tap upper bound; must match run_pf_oltc_step

% Same transient-avoidance logic as for capacitors, but applied only to the
% secondary bus of each transformer (not all load buses).
t_oltc          = 1;     % Minimum number of steps between two tap changes on the same transformer (cooldown)
var_v_oltc_last = 0.01;  % Max allowed voltage variation (p.u.) at the secondary bus over the past t_oltc_last steps
t_oltc_last     = 3;     % Number of past steps used for the voltage stability check
var_v_oltc_next = 0.01;  % Max allowed reference variation (p.u.) over the next t_oltc_next steps
t_oltc_next     = 3;     % Number of future steps used for the reference stability check

%% Preparation
rng(seed);

mpc = case9xx_Bsh(); % Load network data
mpopt = mpoption('verbose',0,'out.all',0,'pf.enforce_q_lims',1);
define_constants;
baseMVA = mpc.baseMVA;
fullJac = true; % Full Jacobian (angles + magnitudes) required by makeJac for sensitivity extraction

% Branch indices of the two 400/63 kV OLTC transformers (bus 10->5 and bus 11->6)
transformer_indices = [16; 17];
n_oltc = length(transformer_indices);


% Bus index sets
pv_idx    = find(mpc.bus(:,BUS_TYPE)==PV);  % Controlled generators (PV buses)
slack_idx = find(mpc.bus(:,BUS_TYPE)==REF); % Slack generator (reference bus, uncontrolled)
load_idx  = find(mpc.bus(:,BUS_TYPE)==PQ);  % Load buses (all PQ buses, including transit buses)
target_idx= find(mpc.bus(:,PD)~=0);         % Subset of load buses with active consumption (voltage tracking targets)

n_load   = size(load_idx,  1);
n_gen    = size(pv_idx,    1);
n_slack  = size(slack_idx, 1);
n_bus    = n_load + n_gen + n_slack;
n_target = size(target_idx,1);
n_state  = n_load + n_gen + n_gen;   % size of state vector x = [V_load; Q_gen; V_gen]

% gen_number(i) = row index in mpc.gen of the generator at bus i (0 if no generator)
gen_number = zeros(n_bus, 1);
for i = 1:n_gen+n_slack
    gen_number(mpc.gen(i,1)) = i;
end

% Hard reactive power limits (from network data, enforced as QP constraints)
Q_gen_max = mpc.gen(pv_idx, QMAX) / baseMVA;
Q_gen_min = mpc.gen(pv_idx, QMIN) / baseMVA;

% Voltage limits for generators and loads (from network data, enforced as QP constraints)
V_gen_max  = mpc.bus(pv_idx,   VMAX);
V_gen_min  = mpc.bus(pv_idx,   VMIN);
V_load_max = mpc.bus(load_idx, VMAX);
V_load_min = mpc.bus(load_idx, VMIN);

% Nominal load and generation values (noise is re-drawn from these each step)
Pd_nom = mpc.bus(target_idx, PD);
Qd_nom = mpc.bus(target_idx, QD);
Pg_nom = mpc.gen(gen_number(pv_idx), PG);

sigma_load_P = Pd_nom * sigma_rate;  % Std of active load noise (MW)
sigma_load_Q = Qd_nom * sigma_rate;  % Std of reactive load noise (Mvar)
sigma_gen_P  = Pg_nom * sigma_rate;  % Std of active generation noise (MW, non-slack only)

% Cooldown counters
cap_availability  = zeros(n_bus,  1);
oltc_availability = zeros(length(transformer_indices), 1);

%% Objective
V_load_min_target=mpc.bus(target_idx,VMIN);
V_load_max_target=mpc.bus(target_idx,VMAX);

% Voltage reference for the target buses (buses with active load)
% Option 1: fixed reference (uncomment below)
V_target=[1.03;1.02;0.98];% Fixed reference for buses 5, 6 and 8
% Option 2: random reference drawn uniformly within each bus's [Vmin, Vmax]
%V_target=V_load_min_target+diag(rand(size(target_idx)))*(V_load_max_target-V_load_min_target);

% Reference trajectories over time (constant here; structure supports time-varying references)
% Yc_V(i,k): voltage reference for target bus i at step k (1-indexed)
% Yc_Q(i,k): reactive power reference for generator i at step k (zero = minimise Q_gen)
% Extra columns beyond nt allow action_oltc and add_cap to look ahead up to t_cap_next/t_oltc_next steps
Yc_V=repmat(V_target,1,nt+max(N,t_cap_next)+1);
Yc_Q=zeros(n_gen,nt+max(N,t_cap_next)+1);

% Yc: reference vector stacked for the QP over the horizon
% At each step k, the relevant block is Yc(k*(n_target+n_gen)+1 : (k+N)*(n_target+n_gen))
% Format per step: [V_target (n_target x 1); zeros(n_gen, 1)]
Yc=repmat([V_target;zeros(n_gen,1)],nt+max(N,t_cap_next)+1,1);

%% Constant MPC matrices (independent of operating point, computed once)
[GAMMA_cst,f1_cst,f2_cst,b1,b2,Hu,Hx,PSI,C_tilde,A,C,J1] = cal_mat_cst(mpc,N,u_max,alpha);

%% Initial power flow
results = runpf(mpc, mpopt);

% State vector: x = [V_load (n_load x 1); Q_gen (n_gen x 1); V_gen (n_gen x 1)]
x = [results.bus(load_idx, VM); results.gen(gen_number(pv_idx), QG)/baseMVA; results.bus(pv_idx, VM)];

% Initial Jacobian
JAC = makeJac(results, fullJac);

% Logging (pre-allocated; column k+1 = state after step k, column 1 = initial state)
X_log     = zeros(n_state, nt+1);   X_log(:, 1) = x;
u_log     = zeros(n_gen,   nt);
J_log     = zeros(4,       nt+1);
tap_log   = zeros(length(transformer_indices), nt+2);
tap_log(:, 1:2) = repmat(mpc.branch(transformer_indices, TAP), 1, 2);
cap_log   = zeros(n_bus, 1);        % dynamic: add_cap appends one column per step
Pd_log    = zeros(n_target, nt);
Qd_log    = zeros(n_target, nt);
Pg_log    = zeros(n_gen,    nt);
V_sec_log = zeros(length(transformer_indices), nt);  % target voltages given to run_pf_oltc_step

[J_Q, J_V, J_VQ] = Cal_criterion(results, C*x, Yc_V, Yc_Q, 1, alpha);
J_log(:, 1) = [0; J_Q; J_V; J_VQ];

%% Control loop
% oltc_idx and direction carry the tap decision from the previous step.
% Initialised to 0 so that no tap correction is applied at k=0 (no prior tap action).
oltc_idx  = 0;
direction = 0;

for k = 0:nt-1

    % State-dependent MPC matrices from the current Jacobian
    [Aineq, H, f1, f2, B] = Cal_mat_var(mpc, N, GAMMA_cst, Hu, Hx, PSI, C_tilde, f1_cst, f2_cst, JAC);

    % Tap correction: the tap decided at k-1 was advanced by run_pf_oltc_step AFTER the power
    % flow, so x(k) reflects the network without that tap increment. x_corrected anticipates
    % the state the MPC will actually face once the in-flight tap takes effect.
    [Cv_tap, Cq_tap] = Cal_Cv_tapCq_tap(mpc, x, transformer_indices, JAC);
    if oltc_idx ~= 0
        x_corrected = x + [Cv_tap(load_idx, transformer_indices(oltc_idx)); ...
                            Cq_tap(pv_idx,   transformer_indices(oltc_idx)); ...
                            zeros(n_gen, 1)] * step_size * direction;
    else
        x_corrected = x;
    end

    % --- Priority 1: MPC (generator voltage) ---
    % QP solved from x_corrected so the MPC accounts for the pending tap effect.
    Yc_k = Yc(k*(n_target+n_gen)+1 : (k+N)*(n_target+n_gen), 1);
    [U, J_qp, exitflag] = quadprog(H, f1*x_corrected - f2*Yc_k, Aineq, b1 - b2*x_corrected, ...
        [], [], [], [], [], optimset('Display', 'off'));
    if exitflag <= 0
        warning('QP failed at step %d — applying zero input', k);
        U = zeros(n_gen*N, 1);
    end
    u_k  = U(1:n_gen);
    u_k1 = U(n_gen+1 : 2*n_gen);
    u_log(:, k+1) = u_k;

    % Reconstruct full MPC cost for logging
    % J_full = e_free'*PSI*e_free + J_qp_unscaled, where quadprog uses 0.5*U'HU convention.
    % Here we approximate: J_log stores the quadprog value plus the free-response term.
    e_free = J1*x - Yc_k;
    JMPC   = J_qp + e_free' * PSI * e_free;

    % State predictions (used by OLTC and capacitor logic)
    x_estimated_k1 = A * x_corrected + B * u_k;
    x_estimated_k2 = A * x_estimated_k1 + B * u_k1;

    % --- Priority 2: OLTC tap decision ---
    % Evaluates all candidate tap actions at k+2 and selects the best one.
    % If a tap moves, x_estimated_k2 is updated so the capacitor logic sees a consistent state.
    [oltc_idx, direction, target_v, x_estimated_k2, oltc_availability] = action_oltc( ...
        mpc, JAC, x_estimated_k2, transformer_indices, step_size, tol, alpha_tap, Q_lim, ...
        Yc_V, k, tap_min, tap_max, oltc_availability, t_oltc, var_v_oltc_last, t_oltc_last, ...
        var_v_oltc_next, t_oltc_next, X_log(:, 1:k+1));

    % --- Priority 3: Capacitor/coil decision ---
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

    % Read actual state
    x = [results.bus(load_idx, VM); results.gen(gen_number(pv_idx), QG)/baseMVA; results.bus(pv_idx, VM)];

    % Log
    X_log(:, k+2)      = x;
    tap_log(:, k+3)    = final_taps;
    V_sec_log(:, k+1)  = target_v;
    Pd_log(:, k+1)     = results.bus(target_idx, PD);
    Qd_log(:, k+1)     = results.bus(target_idx, QD);
    Pg_log(:, k+1)     = results.gen(gen_number(pv_idx), PG);

    [J_Q, J_V, J_VQ] = Cal_criterion(results, C*x, Yc_V, Yc_Q, k+1, alpha);
    J_log(:, k+2) = [JMPC; J_Q; J_V; J_VQ];


    % Update Jacobian (capacitor included; tap not yet applied)
    if t_jac == 0 || mod(k+1, t_jac) == 0
        JAC = makeJac(results, fullJac);
    end
end

%% DISPLAY

% %% Case 
% figure;
% img = imread('case.png');
% imshow(img);
% title('Case considered');
% 
% %% Voltage in loads
% figure;
% zoom on; hold on; grid on;
% 
% hX = stairs(0:nt, X_log(1:n_load,:)');%Measured
% set(hX, 'LineWidth', 1.5)
% 
% hY = stairs(1:nt, Yc_V(1:n_target,1:nt)', 'x');%Target
% for i =1:n_target
%     j=find(load_idx==target_idx(i));
%     hY(i).Color = hX(j).Color;
%     hY(i).LineStyle = 'none';   
%     hY(i).Marker = 'x';        
%     hY(i).LineWidth = 0.8;
% end
% 
% hZ = stairs(1:nt,[repmat(V_load_max,1,nt);repmat(V_load_min,1,nt)]');%Constraints
% for i = 1:n_load
%     hZ(i).Color = hX(i).Color;
%     hZ(i+n_load).Color = hX(i).Color;
%     hZ(i).LineStyle = 'none';  
%     hZ(i+n_load).LineStyle = 'none';   
%     hZ(i).Marker = 'diamond'; 
%     hZ(i+n_load).Marker = 'diamond';        
%     hZ(i).LineWidth = 0.8;
%     hZ(i+n_load).LineWidth = 0.8;
% end
% 
% legendHandles = [hX(:)',hY(1),hZ(1)];       
% legendStrings = [compose("Load %d", load_idx'),"Target","Constraints" ];
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% 
% title('Voltage in loads');
% xlabel('k');
% ylabel('Voltage (p.u)')
% 
% %% Voltage in targeted loads
% figure;
% zoom on; hold on; grid on;
% 
% idx = find(ismember(load_idx, target_idx));
% hX = stairs(0:nt, X_log(idx,:)');%Measured
% set(hX, 'LineWidth', 1.5)
% 
% hY = stairs(1:nt, Yc_V(1:n_target,1:nt)', 'x');%Target
% for i =1:n_target
%     hY(i).Color = hX(i).Color;
%     hY(i).LineStyle = 'none';   
%     hY(i).Marker = 'x';        
%     hY(i).LineWidth = 0.8;
% end
% 
% hZ = stairs(1:nt,[repmat(V_load_max(idx),1,nt);repmat(V_load_min(idx),1,nt)]');%Constraints
% for i = 1:n_target
%     hZ(i).Color = hX(i).Color;
%     hZ(i+n_target).Color = hX(i).Color;
%     hZ(i).LineStyle = 'none';  
%     hZ(i+n_target).LineStyle = 'none';   
%     hZ(i).Marker = 'diamond'; 
%     hZ(i+n_target).Marker = 'diamond';        
%     hZ(i).LineWidth = 0.8;
%     hZ(i+n_target).LineWidth = 0.8;
% end
% 
% legendHandles = [hX(:)',hY(1),hZ(1)];       
% legendStrings = [compose("Load %d", target_idx'),"Target","Constraints" ];
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% 
% title('Voltage in targeted loads');
% xlabel('k');
% ylabel('Voltage (p.u)')
% 
% %% Reactive power in generators
% figure;
% zoom on; hold on; grid on;
% 
% hX = stairs(0:nt,X_log(n_load+1:n_load+n_gen,:)');%Measured
% set(hX, 'LineWidth', 1.5)
% 
% hZ = stairs(0:nt,[repmat(Q_gen_max,1,nt+1);repmat(Q_gen_min,1,nt+1)]');%Constraints
% for i = 1:n_gen
%     hZ(i).Color = hX(i).Color;
%     hZ(i+n_gen).Color = hX(i).Color;
%     hZ(i).LineStyle = 'none';  
%     hZ(i+n_gen).LineStyle = 'none';   
%     hZ(i).Marker = 'diamond'; 
%     hZ(i+n_gen).Marker = 'diamond';        
%     hZ(i).LineWidth = 0.8;
%     hZ(i+n_gen).LineWidth = 0.8;
% end
% 
% if add_cap_bool
%     hY=stairs(0:nt,[-ones(n_gen,nt+1)*Q_lim;ones(n_gen,nt+1)*Q_lim]'); %Targeted limits
%     for i = 1:n_gen
%         hY(i).Color = hX(i).Color;
%         hY(i+n_gen).Color = hX(i).Color;
%         hY(i).LineStyle = 'none';  
%         hY(i+n_gen).LineStyle = 'none';   
%         hY(i).Marker = 'x'; 
%         hY(i+n_gen).Marker = 'x';        
%         hY(i).LineWidth = 2;
%         hY(i+n_gen).LineWidth = 2;
%     end 
% 
%     legendHandles = [hX(:)',hY(1),hZ(1)];       
%     legendStrings = [compose("Gen %d", pv_idx'),"Targeted limits","Constraints"];
% else
%     legendHandles = [hX(:)',hZ(1)];       
%     legendStrings = [compose("Gen %d", pv_idx'),"Constraints"];
% 
% end
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% 
% title('Reactive power in generators');
% xlabel('k');
% ylabel('Reactive power (p.u)')
% 
% %% Voltage in generators
% figure;
% zoom on; hold on ;grid on;
% 
% hX = stairs(0:nt,X_log(n_load+n_gen+1:end,:)'); %Measured
% set(hX, 'LineWidth', 1.5)
% 
% hZ = stairs(0:nt,[repmat(V_gen_max,1,nt+1);repmat(V_gen_min,1,nt+1)]');%Constraints
% for i = 1:n_gen
%     hZ(i).Color = hX(i).Color;
%     hZ(i+n_gen).Color = hX(i).Color;
%     hZ(i).LineStyle = 'none';  
%     hZ(i+n_gen).LineStyle = 'none';   
%     hZ(i).Marker = 'diamond'; 
%     hZ(i+n_gen).Marker = 'diamond';        
%     hZ(i).LineWidth = 0.8;
%     hZ(i+n_gen).LineWidth = 0.8;
% end
% 
% legendHandles = [hX(:)',hZ(1)];       
% legendStrings = [compose("Gen %d", pv_idx'),"Constraints" ];
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% 
% title('Voltage in generators');
% xlabel('k');
% ylabel('Voltage (p.u)')
% 
% %% Voltage input
% figure;
% zoom on; hold on ;grid on;
% 
% hX=stairs(0:(nt-1),u_log');%Measured
% set(hX, 'LineWidth', 1.5)
% 
% hZ=stairs(0:(nt-1),[-u_max*ones(N,nt);u_max*ones(N,nt)]'); %Constraints
% for i = 1:n_gen
%     hZ(i).Color = hX(i).Color;
%     hZ(i+n_gen).Color = hX(i).Color;
%     hZ(i).LineStyle = 'none';  
%     hZ(i+n_gen).LineStyle = 'none';   
%     hZ(i).Marker = 'diamond'; 
%     hZ(i+n_gen).Marker = 'diamond';        
%     hZ(i).LineWidth = 0.8;
%     hZ(i+n_gen).LineWidth = 0.8;
% end
% 
% legendHandles = [hX(:)',hZ(1)];       
% legendStrings = [compose("Gen %d", pv_idx'),"Constraints" ];
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% 
% title('Voltage input');
% xlabel('k');
% ylabel('Voltage input (p.u)')
% 
% 
% %% Capacitor
% if add_cap_bool
%     figure;
%     zoom on; hold on ;grid on;
% 
%     hX=stairs(0:nt,cap_log(target_idx,:)');
%     set(hX, 'LineWidth', 1.5)
% 
%     legendStrings = [compose("Load %d", target_idx')];
%     legend(legendStrings{:});
% 
%     title('Capacitors');
%     xlabel('k');
%     ylabel('Capacitors')
% end
% 
% %% Criterion
% figure;
% 
% % 1. Criterion MPC
% subplot(2,2,1);
% zoom on; hold on; grid on;
% 
% hX = stairs(1:nt, J_log(1,2:end)');
% set(hX, 'LineWidth', 1.5);
% 
% title('MPC Criterion');
% xlabel('k'); ylabel('Criterion MPC');
% 
% % 2. Criterion Q^2
% subplot(2,2,2);
% zoom on; hold on; grid on;
% 
% hX = stairs(0:nt, J_log(2,:)');
% set(hX, 'LineWidth', 1.5);
% 
% title('Criterion on Q :');
% xlabel('k'); 
% ylabel('sum Q^2')
% 
% % 3. Criterion v^2
% subplot(2,2,3);
% zoom on; hold on; grid on;
% 
% hX = stairs(0:nt, J_log(3,:)');
% set(hX, 'LineWidth', 1.5);
% 
% title('Criterion on V :');
% xlabel('k'); 
% ylabel('sum (V-V_{target})^2')
% 
% % 4. Criterion on V and Q
% subplot(2,2,4);
% zoom on; hold on; grid on;
% 
% hX = stairs(0:nt, J_log(4,:)');
% set(hX, 'LineWidth', 1.5);
% 
% title('Criterion on V and Q :');
% xlabel('k'); 
% ylabel('$\sum (V - V_{\mathrm{target}})^2 + \alpha \sum Q^2$', ...
%        'Interpreter','latex')
% 
% %% Tap
% figure
% zoom on; hold on; grid on;
% 
% hX = stairs(0:nt+1, tap_log');%Measured
% set(hX, 'LineWidth', 1.5)
% 
% legendHandles = hX(:)';       
% legendStrings = [compose("transformer %d-%d (branch %d)", [mpc.branch(transformer_indices,1:2),transformer_indices])];
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% 
% title('Tap')
% xlabel('k')
% ylabel('Tap')

% %% Active and reactive powers (noise visualisation)
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
hX = stairs(0:nt, X_log(idx,:)');%Measured
set(hX, 'LineWidth', 1.5)

hY = stairs(1:nt, Yc_V(1:n_target,1:nt)', 'x');%Target
for i =1:n_target
    hY(i).Color = hX(i).Color;
    hY(i).LineStyle = 'none';   
    hY(i).Marker = 'x';        
    hY(i).LineWidth = 0.8;
end

hZ = stairs(1:nt,[repmat(V_load_max(idx),1,nt);repmat(V_load_min(idx),1,nt)]');%Constraints
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
legendStrings = [compose("Load %d", target_idx'),"Target","Constraints" ];
legend(legendHandles, legendStrings, 'Location', 'northeast');

title('Voltage in targeted loads');
xlabel('k');
ylabel('Voltage (p.u)')

% 2. Reactive in generators
subplot(3,3,2);
zoom on; hold on; grid on;

hX = stairs(0:nt,X_log(n_load+1:n_load+n_gen,:)');%Measured
set(hX, 'LineWidth', 1.5)

hZ = stairs(0:nt,[repmat(Q_gen_max,1,nt+1);repmat(Q_gen_min,1,nt+1)]');%Constraints
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
        hY=stairs(0:nt,[-ones(n_gen,nt+1)*Q_lim;ones(n_gen,nt+1)*Q_lim]'); %Targeted limits
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
xlabel('k');
ylabel('Reactive power (p.u)')

% 3. Voltage in generators
subplot(3,3,3);
zoom on; hold on; grid on;

hX = stairs(0:nt,X_log(n_load+n_gen+1:end,:)'); %Measured
set(hX, 'LineWidth', 1.5)

hZ = stairs(0:nt,[repmat(V_gen_max,1,nt+1);repmat(V_gen_min,1,nt+1)]');%Constraints
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
legendStrings = [compose("Gen %d", pv_idx'),"Constraints" ];
legend(legendHandles, legendStrings, 'Location', 'northeast');

title('Voltage in generators');
xlabel('k');
ylabel('Voltage (p.u)')

%. 4. Voltage in loads
subplot(3,3,4)
zoom on; hold on; grid on;

hX = stairs(0:nt, X_log(1:n_load,:)');%Measured
set(hX, 'LineWidth', 1.5)

  hZ = stairs(1:nt,[repmat(V_load_max,1,nt);repmat(V_load_min,1,nt)]');%Constraints
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
legendStrings = [compose("Load %d", load_idx'),"Constraints" ];
legend(legendHandles, legendStrings, 'Location', 'northeast');

title('Voltage in loads');
xlabel('k');
ylabel('Voltage (p.u)')

% 5. Case
subplot(3,3,5)

img = imread('case.png');
imshow(img);
title('Case considered');

%6. Voltage input
subplot(3,3,6)

zoom on; hold on ;grid on;

hX=stairs(0:(nt-1),u_log');%Measured
set(hX, 'LineWidth', 1.5)

hZ=stairs(0:(nt-1),[-u_max*ones(N,nt);u_max*ones(N,nt)]'); %Constraints
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
legendStrings = [compose("Gen %d", pv_idx'),"Constraints" ];
legend(legendHandles, legendStrings, 'Location', 'northeast');

title('Voltage input');
xlabel('k');
ylabel('Voltage input (p.u)')

% 7. TAP
subplot(3,3,7);
zoom on; hold on ;grid on;

hX = stairs(0:nt+1, tap_log');%Measured
set(hX, 'LineWidth', 1.5)

legendHandles = hX(:)';       
legendStrings = [compose("transformer %d-%d (branch %d)", [mpc.branch(transformer_indices,1:2),transformer_indices])];
legend(legendHandles, legendStrings, 'Location', 'northeast');

title('Tap')
xlabel('k')
ylabel('Tap')

% 8. Criterion on V
subplot(3,3,8)
zoom on; hold on ;grid on;

hX = stairs(0:nt, J_log(3,:)');
set(hX, 'LineWidth', 1.5);

x_last = nt;
y_last = J_log(3,nt);

plot(x_last, y_last, 'ro','MarkerFaceColor','r');
text(x_last, y_last, sprintf(' %.3e', y_last), ...
    'VerticalAlignment','bottom','HorizontalAlignment','left');

title('Criterion on V :');
xlabel('k'); 
ylabel('sum (V-V_{target})^2')

% 9. Capacitor
if add_cap_bool
    subplot(3,3,9);
    zoom on; hold on ;grid on;

    hX=stairs(0:nt,cap_log(target_idx,:)');
    set(hX, 'LineWidth', 1.5)

    legendStrings = [compose("Load %d", target_idx')];
    legend(legendStrings{:});

    title('Capacitors');
    xlabel('k');
    ylabel('Capacitors')
end

% %% DISPLAY REPORT — Axis 1: controller comparison
% 
% fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures', 'ComparaisonVersion');
% if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
% prefix = 'alt';
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