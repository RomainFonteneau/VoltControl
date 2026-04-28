clear all

%% Parameters definition
seed = 1234;
nt = 50;   % Number of simulation steps
N  = 3;    % MPC prediction horizon (steps)

% MPC cost weights (presentation notation):
%   min sum_i [ ||V_load - V_ref||^2 + alpha*||Q_gen||^2 + beta*||u_OLTC||^2 ]
alpha = 0.001;  % Weight on Q_gen: voltage tracking dominates (V term has coefficient 1)
beta  = 1;      % Weight on OLTC magnitude: tap used only when generators insufficient

u_max      = 0.2;   % Maximum generator voltage increment per step (p.u.)


% Gaussian noise standard deviations
sigma_rate = 0;  % Rate used to compute the std deviation of noise around the nominal

% Jacobian recalculation period: 0 = every step, n = every n steps
t_jac = 0;

%----------Capacitor (lowest priority actuator)------------
add_cap_bool  = true;
Q_lim         = 0.2;    % Soft reactive power limit (p.u.)
t_cap          = 5;     % Minimum number of steps between two switches on the same bus (cooldown)
var_v_cap_last = 0.01;  % Max allowed voltage variation (p.u.) at the target bus over the past t_cap_last steps
t_cap_last     = 3;     % Number of past steps used for the voltage stability check
var_v_cap_next = 0.01;  % Max allowed reference variation (p.u.) over the next t_cap_next steps
t_cap_next     = 3;     % Number of future steps used for the reference stability check
cap            = 5;     % Capacitor/coil size (Mvar); adding a coil is equivalent to cap = -cap

%---------OLTC (integrated as MPC input)------------
% run_pf_oltc_step parameters — must match the function internals
step_size = 0.00625;  % Tap step size (p.u.)
tol       = 0.005;    % Voltage deadband for target_v encoding in run_pf_oltc_step
tap_min   = 0.9;      % Physical tap lower bound
tap_max   = 1.1;      % Physical tap upper bound

% Threshold for converting the continuous QP output u_OLTC to a discrete tap request.
% If |u_OLTC(j)| > epsilon_oltc, a tap movement is requested on transformer j.
epsilon_oltc = step_size / 2;

% Transient protection: tap commands are blocked by hard QP constraints when the
% secondary bus voltage has been unstable recently or the reference will change soon.
t_oltc          = 0;     % Cooldown steps after each tap action
var_v_oltc_last = 0.01;  % Max voltage deviation (p.u.) at secondary bus over past t_oltc_last steps
t_oltc_last     = 3;     % Number of past steps checked for voltage stability
var_v_oltc_next = 0.01;  % Max reference variation (p.u.) over next t_oltc_next steps
t_oltc_next     = 3;     % Number of future steps checked for reference stability

%% Preparation
rng(seed);

mpc   = case9xx_Bsh();
mpopt = mpoption('verbose', 0, 'out.all', 0, 'pf.enforce_q_lims', 1);
define_constants;
baseMVA = mpc.baseMVA;
fullJac = true;

% Branch indices of the two 400/63 kV OLTC transformers (bus 10->5 and bus 11->6)
transformer_indices = [16; 17];
n_oltc = length(transformer_indices);

% Secondary (LV) buses of OLTC transformers
reg_buses = mpc.branch(transformer_indices, T_BUS);

% Bus index sets
pv_idx    = find(mpc.bus(:, BUS_TYPE) == PV);
slack_idx = find(mpc.bus(:, BUS_TYPE) == REF);
load_idx  = find(mpc.bus(:, BUS_TYPE) == PQ);
target_idx = find(mpc.bus(:, PD) ~= 0);

n_load   = size(load_idx,  1);
n_gen    = size(pv_idx,    1);
n_slack  = size(slack_idx, 1);
n_bus    = n_load + n_gen + n_slack;
n_target = size(target_idx, 1);
n_state  = n_load + n_gen + n_gen;   % original (non-augmented) state size
n_input  = n_gen + n_oltc;           % extended input size

% gen_number(i) = row index in mpc.gen of the generator at bus i (0 if none)
gen_number = zeros(n_bus, 1);
for i = 1:n_gen + n_slack
    gen_number(mpc.gen(i,1)) = i;
end

% Hard limits
Q_gen_max  = mpc.gen(pv_idx, QMAX) / baseMVA;
Q_gen_min  = mpc.gen(pv_idx, QMIN) / baseMVA;
V_gen_max  = mpc.bus(pv_idx, VMAX);
V_gen_min  = mpc.bus(pv_idx, VMIN);
V_load_max = mpc.bus(load_idx, VMAX);
V_load_min = mpc.bus(load_idx, VMIN);

% Nominal load and generation values (noise is re-drawn from these each step)
Pd_nom = mpc.bus(target_idx, PD);
Qd_nom = mpc.bus(target_idx, QD);
Pg_nom = mpc.gen(gen_number(pv_idx), PG);

sigma_load_P = Pd_nom*sigma_rate;   % Std of active load noise (MW)
sigma_load_Q = Qd_nom*sigma_rate;   % Std of reactive load noise (Mvar)
sigma_gen_P  = Pg_nom*sigma_rate;   % Std of active generation noise (MW, non-slack generators only)

% Position of OLTC secondary buses within load_idx (for transient detection)
[~, reg_pos] = ismember(reg_buses, load_idx);

% For each transformer, the index of its secondary bus within target_idx.
% Used to check future reference stability (0 if the bus has no voltage reference).
oltc_in_target = zeros(n_oltc, 1);
for j = 1:n_oltc
    idx = find(target_idx == reg_buses(j));
    if ~isempty(idx)
        oltc_in_target(j) = idx;
    end
end

% Cooldown counters: oltc_availability(j) = steps remaining before transformer j
% can request a tap movement. Implemented as hard QP constraints on u_OLTC(j).
oltc_availability = zeros(n_oltc, 1);
cap_availability  = zeros(n_bus,  1);

% Previous OLTC command, stored in the augmented state to model the two-step tap delay.
% Initialised to zero: no tap command was in flight before the simulation starts.
u_OLTC_prev = zeros(n_oltc, 1);

% Tap limit flags — ONLY used as booleans; the numeric value of final_taps is
% never used in the control logic beyond this comparison.
tap_at_min = mpc.branch(transformer_indices, TAP) <= tap_min;
tap_at_max = mpc.branch(transformer_indices, TAP) >= tap_max;

% Selection matrices for dynamic tap constraints.
% sel_oltc{j} (N x N*n_input): selects u_OLTC(j) from every horizon step of U.
sel_oltc = cell(n_oltc, 1);
for j = 1:n_oltc
    e_j = zeros(n_input, 1);
    e_j(n_gen + j) = 1;
    sel_oltc{j} = kron(eye(N), e_j');   % (N x N*n_input)
end

%% Objective
V_load_min_target = mpc.bus(target_idx, VMIN);
V_load_max_target = mpc.bus(target_idx, VMAX);

% Voltage reference for target buses (constant; structure supports time-varying references)
V_target=[1.03;0.97;0.98];% Fixed reference for buses 5, 6 and 8
%V_target=V_load_min_target+diag(rand(size(target_idx)))*(V_load_max_target-V_load_min_target);

Yc_V = repmat(V_target, 1, nt + max(N, t_cap_next) + 1);
Yc_Q = zeros(n_gen, nt + max(N, t_cap_next) + 1);
Yc   = repmat([V_target; zeros(n_gen,1)], nt + max(N, t_cap_next) + 1, 1);

%% Constant MPC matrices (computed once — now include OLTC dimensions)
[C, PSI, C_tilde, Hu, Hx, b1, H_gamma] = ...
    cal_mat_cst(mpc, N, u_max, alpha, beta, n_oltc, step_size);


% C_orig: output matrix for the non-augmented state x, used in Cal_criterion.
% C has zero columns for the u_OLTC_prev part of x_tilde, so C_orig*x = C*x_tilde.
C_orig = C(:, 1:n_state);   % (n_output x n_state)

%% Initial power flow
results = runpf(mpc, mpopt);

% State vector: x = [V_load (n_load x 1); Q_gen (n_gen x 1); V_gen (n_gen x 1)]
x = [results.bus(load_idx, VM); ...
     results.gen(gen_number(pv_idx), QG) / baseMVA; ...
     results.bus(pv_idx, VM)];

% Initial secondary bus voltages (needed to encode target_v at the first step)
V_reg = results.bus(reg_buses, VM);

% Initial Jacobian
JAC = makeJac(results, fullJac);

% Logging initialisation
X_log      = zeros(n_state,  nt+1);   X_log(:,1) = x;
u_log      = zeros(n_gen,    nt);     % generator voltage increments (n_gen x nt)
u_oltc_log = zeros(n_oltc,   nt);     % OLTC commands from QP        (n_oltc x nt)
J_MPC      = 0;                       % MPC cost placeholder for k=0 (no QP at initial step)
[J_Q, J_V, J_VQ] = Cal_criterion(results, C_orig*x, Yc_V, Yc_Q, 1, alpha);
J_log      = zeros(4, nt+1);          J_log(:,1) = [J_MPC; J_Q; J_V; J_VQ];
tap_log    = zeros(n_oltc, nt+2);     tap_log(:,1:2) = repmat(mpc.branch(transformer_indices, TAP), 1, 2);
cap_log    = zeros(n_bus, 1);         % cap_log remains dynamic (shared interface with add_cap.m)
Pd_log     = zeros(n_target, nt);
Qd_log     = zeros(n_target, nt);
Pg_log     = zeros(n_gen,    nt);
V_sec_log  = zeros(n_oltc,    nt);

%% Control loop
for k = 0:nt-1

    % State-dependent MPC matrices from the current Jacobian
    [Aineq, b2, H, f1, f2, B_tilde, A_tilde, J1] = ...
    Cal_mat_var(mpc, x, N, C_tilde, Hu, Hx, PSI, H_gamma, JAC, transformer_indices);

    % Augmented state: append the in-flight OLTC command (decided at k-1).
    % The augmented model A_tilde propagates u_OLTC_prev into x at this step,
    % correctly accounting for the two-step tap delay without any explicit correction.
    x_tilde = [x; u_OLTC_prev];

    % --- Decrement cooldown counters ---
    oltc_availability = max(oltc_availability - 1, zeros(n_oltc, 1));

    % --- Build dynamic constraints: tap limits + transient protection ---
    Aineq_dyn = [];
    b_dyn     = [];
    for j = 1:n_oltc
        % Past voltage stability at secondary bus
        if k >= t_oltc_last
            recent_V    = X_log(reg_pos(j), k-t_oltc_last+2:k+1);
            past_stable = all(abs(recent_V - x(reg_pos(j))) < var_v_oltc_last);
        else
            past_stable = false;
        end

        % Future reference stability
        if oltc_in_target(j) > 0
            future_refs   = Yc_V(oltc_in_target(j), k+2 : k+t_oltc_next+1);
            future_stable = all(abs(future_refs - Yc_V(oltc_in_target(j), k+2)) < var_v_oltc_next);
        else
            future_stable = true;
        end

        tap_blocked = (oltc_availability(j) > 0) || ~past_stable || ~future_stable;

        if tap_blocked
            % Force u_OLTC(j) = 0 for all N horizon steps
            Aineq_dyn = [Aineq_dyn;  sel_oltc{j}; -sel_oltc{j}];
            b_dyn     = [b_dyn;      zeros(2*N, 1)];
        elseif tap_at_max(j)
            % Tap at upper limit: u_OLTC(j) <= 0
            Aineq_dyn = [Aineq_dyn;  sel_oltc{j}];
            b_dyn     = [b_dyn;      zeros(N, 1)];
        elseif tap_at_min(j)
            % Tap at lower limit: u_OLTC(j) >= 0  (encoded as -u_OLTC(j) <= 0)
            Aineq_dyn = [Aineq_dyn; -sel_oltc{j}];
            b_dyn     = [b_dyn;      zeros(N, 1)];
        end
    end

    Aineq_full = [Aineq;           Aineq_dyn];
    b_full     = [b1 - b2*x_tilde; b_dyn];

    % --- Priority 1: MPC (generator voltage + OLTC) ---
    [U, J_qp, exitflag] = quadprog(H, f1*x_tilde - f2*Yc(k*(n_target+n_gen)+1:(k+N)*(n_target+n_gen), 1), ...
        Aineq_full, b_full, [], [], [], [], [], optimset('Display', 'off'));
    if exitflag <= 0
        warning('QP failed at step %d — applying zero input', k);
        U = zeros(n_input * N, 1);
    end

    % Receding horizon: apply only the first step of the optimal sequence
    u_k         = U(1:n_gen);
    u_OLTC_star_k = U(n_gen+1:n_input);
    u_k1 = U(n_input+1:n_input+n_gen);
    u_OLTC_star_k1=U(n_input+n_gen+1:2*n_input);

    u_log(:, k+1)      = u_k;
    u_oltc_log(:, k+1) = u_OLTC_star_k;

    % Reconstruct full MPC cost for logging
    % J_full = e_free'*PSI*e_free + 2*J_qp
    % where e_free = J1*x_tilde - Yc is the free-response tracking error (U=0)
    % and the factor 2 corrects for quadprog's 0.5 convention (J_qp = 0.5*U'HU + f'U).
    e_free = J1*x_tilde - Yc(k*(n_target+n_gen)+1:(k+N)*(n_target+n_gen), 1);
    JMPC = e_free' * PSI * e_free + 2*J_qp;

    % Convert continuous u_OLTC_star_k to discrete value, compute target_v and reset
    % cooldown
    target_v = V_reg;
    for j = 1:n_oltc
        if u_OLTC_star_k(j) > epsilon_oltc
            target_v(j) = V_reg(j) - 2*tol;   % tap increase requested
            u_OLTC_star_k(j)=step_size;
            oltc_availability(j) = t_oltc; %Reset cooldown
        elseif u_OLTC_star_k(j) < -epsilon_oltc
            target_v(j) = V_reg(j) + 2*tol;   % tap decrease requested
            u_OLTC_star_k(j)=-step_size;
            oltc_availability(j) = t_oltc; %Reset cooldown
        end
    end

    % Convert continuous u_OLTC_star_k1 to discrete value
     for j = 1:n_oltc
        if u_OLTC_star_k1(j) > epsilon_oltc
            u_OLTC_star_k1(j)=step_size;
        elseif u_OLTC_star_k1(j) < -epsilon_oltc
            u_OLTC_star_k1(j)=-step_size;
        end
    end
    V_sec_log(:,k+1)=target_v;


    % State predictions for add_cap
    u_full_k   = [u_k; u_OLTC_star_k];
    x_tilde_k1 = A_tilde * x_tilde + B_tilde * u_full_k;
    x_estimated_k1 = x_tilde_k1(1:n_state);

    u_full_k1  = [u_k1; u_OLTC_star_k1];
    x_tilde_k2 = A_tilde * x_tilde_k1 + B_tilde * u_full_k1;
    x_estimated_k2 = x_tilde_k2(1:n_state);

    % --- Priority 2: Capacitor/coil decision ---
    % use_add_cap is a pre-filter: skips the call (and the costly Cal_Cv_capCq_cap inside)
    % when all generators are already within [-Q_lim, Q_lim]. add_cap.m would reach the
    % same conclusion internally, but this avoids the Jacobian-based sensitivity computation.
    use_add_cap = ~(all(x_estimated_k2(n_load+1:n_load+n_gen) < Q_lim) && ...
                    all(x_estimated_k2(n_load+1:n_load+n_gen) > -Q_lim));
    if k < t_cap_last || ~use_add_cap || ~add_cap_bool
        cap_log = [cap_log, cap_log(:, end)];
        cap_availability = max(cap_availability - ones(n_bus,1), zeros(n_bus,1));
    else
        [mpc, cap_availability, cap_log] = add_cap(mpc, JAC, x_estimated_k1, x_estimated_k2, ...
            cap_availability, cap_log, t_cap, var_v_cap_last, t_cap_last, var_v_cap_next, t_cap_next, ...
            cap, k, X_log, Yc_V, Q_lim);
    end

    % Apply MPC command: increment generator voltage setpoints, clamped to physical limits
    mpc.gen(gen_number(pv_idx), VG) = min(max(mpc.gen(gen_number(pv_idx), VG) + u_k, V_gen_min), V_gen_max);

    % Inject Gaussian noise on loads and non-slack generation
    mpc.bus(target_idx, PD) = Pd_nom + sigma_load_P .* randn(length(target_idx), 1);
    mpc.bus(target_idx, QD) = Qd_nom + sigma_load_Q .* randn(length(target_idx), 1);
    mpc.gen(gen_number(pv_idx), PG) = Pg_nom + sigma_gen_P .* randn(n_gen, 1);

    % Run power flow and advance OLTC tap one step toward target_v
    [mpc, results, final_taps] = run_pf_oltc_step(mpc, target_v, transformer_indices);

    % Read actual state
    x = [results.bus(load_idx, VM); ...
         results.gen(gen_number(pv_idx), QG) / baseMVA; ...
         results.bus(pv_idx, VM)];

    % Update V_reg and tap limit flags
    V_reg      = results.bus(reg_buses, VM);
    tap_at_min = final_taps <= tap_min;
    tap_at_max = final_taps >= tap_max;

    % Store command as in-flight tap for the next step's augmented state
    u_OLTC_prev = u_OLTC_star_k;

    X_log(:, k+2)   = x;
    tap_log(:, k+3) = final_taps;    % offset +3: tap_log has 2 pre-filled initial columns
    Pd_log(:, k+1)  = results.bus(target_idx, PD);
    Qd_log(:, k+1)  = results.bus(target_idx, QD);
    Pg_log(:, k+1)  = results.gen(gen_number(pv_idx), PG);

    [J_Q, J_V, J_VQ] = Cal_criterion(results, C_orig*x, Yc_V, Yc_Q, k+1, alpha);
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
% hold on; grid on;
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
% hold on; grid on;
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
% hold on; grid on;
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
%     hY=stairs(0:nt,[-ones(n_gen,nt+1)*Q_lim;ones(n_gen,nt+1)*Q_lim]');
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
%     legendHandles = [hX(:)',hY(1),hZ(1)];       
%     legendStrings = [compose("Gen %d", pv_idx'),"Targeted limits","Constraints"];
% else
%     legendHandles = [hX(:)',hZ(1)];       
%     legendStrings = [compose("Gen %d", pv_idx'),"Constraints"];
% end
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% title('Reactive power in generators');
% xlabel('k');
% ylabel('Reactive power (p.u)')
% 
% %% Voltage in generators
% figure;
% hold on ;grid on;
% hX = stairs(0:nt,X_log(n_load+n_gen+1:end,:)');
% set(hX, 'LineWidth', 1.5)
% hZ = stairs(0:nt,[repmat(V_gen_max,1,nt+1);repmat(V_gen_min,1,nt+1)]');
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
% legendHandles = [hX(:)',hZ(1)];       
% legendStrings = [compose("Gen %d", pv_idx'),"Constraints" ];
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% title('Voltage in generators');
% xlabel('k');
% ylabel('Voltage (p.u)')
% 
% %% Voltage input (generators)
% figure;
% hold on ;grid on;
% hX=stairs(0:(nt-1),u_log');
% set(hX, 'LineWidth', 1.5)
% hZ=stairs(0:(nt-1),[-u_max*ones(N,nt);u_max*ones(N,nt)]');
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
% legendHandles = [hX(:)',hZ(1)];       
% legendStrings = [compose("Gen %d", pv_idx'),"Constraints" ];
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% title('Voltage input');
% xlabel('k');
% ylabel('Voltage input (p.u)')
% 
% %% OLTC commands from QP
% figure;
% hold on; grid on;
% hX = stairs(0:(nt-1), u_oltc_log');
% set(hX, 'LineWidth', 1.5)
% plot([0 nt-1], [epsilon_oltc epsilon_oltc],  'k--', 'LineWidth', 1);
% plot([0 nt-1], [-epsilon_oltc -epsilon_oltc], 'k--', 'LineWidth', 1);
% legendHandles = hX(:)';
% legendStrings = [compose("OLTC %d (branch %d)", [(1:n_oltc)', transformer_indices])];
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% title('OLTC commands from QP');
% xlabel('k'); ylabel('u_{OLTC} (p.u.)')
% 
% %% Capacitor
% if add_cap_bool
%     figure;
%     hold on ;grid on;
%     hX=stairs(0:nt,cap_log(target_idx,:)');
%     set(hX, 'LineWidth', 1.5)
%     legendStrings = [compose("Load %d", target_idx')];
%     legend(legendStrings{:});
%     title('Capacitors');
%     xlabel('k');
%     ylabel('Capacitors')
% end
% 
% %% Active and reactive powers
% figure;
% 
% % 1. Active load power
% subplot(3,1,1);
% hold on; grid on;
% hX = stairs(1:nt, Pd_log');
% set(hX, 'LineWidth', 1.5);
% hN = plot([1 nt], repmat(Pd_nom', 2, 1), '--', 'LineWidth', 1);
% for i = 1:n_target
%     hN(i).Color = hX(i).Color;
% end
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
% for i = 1:n_target
%     hN(i).Color = hX(i).Color;
% end
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
% for i = 1:n_gen
%     hN(i).Color = hX(i).Color;
% end
% legend(hX(:)', compose("Gen %d", pv_idx'), 'Location', 'northeast');
% title('Active generation power (non-slack)');
% xlabel('k'); ylabel('P_{gen} (MW)');

%% Summary
figure;

% 1. Voltage in targeted loads
subplot(3,3,1);
hold on; grid on;
idx = find(ismember(load_idx, target_idx));
hX = stairs(0:nt, X_log(idx,:)');
set(hX, 'LineWidth', 1.5)
hY = stairs(1:nt, Yc_V(1:n_target,1:nt)', 'x');
for i =1:n_target
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
legendStrings = [compose("Load %d", target_idx'),"Target","Constraints" ];
legend(legendHandles, legendStrings, 'Location', 'northeast');
title('Voltage in targeted loads');
xlabel('k');
ylabel('Voltage (p.u)')

% 2. Reactive in generators
subplot(3,3,2);
hold on; grid on;
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
    hY=stairs(0:nt,[-ones(n_gen,nt+1)*Q_lim;ones(n_gen,nt+1)*Q_lim]');
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
hold on; grid on;
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
legendStrings = [compose("Gen %d", pv_idx'),"Constraints" ];
legend(legendHandles, legendStrings, 'Location', 'northeast');
title('Voltage in generators');
xlabel('k');
ylabel('Voltage (p.u)')

% 4. Voltage in loads or capacitors
subplot(3,3,4)
hold on; grid on;

if not(add_cap_bool)
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
    legendStrings = [compose("Load %d", load_idx'),"Constraints" ];
    legend(legendHandles, legendStrings, 'Location', 'northeast');
    title('Voltage in loads');
    xlabel('k');
    ylabel('Voltage (p.u)')
else 
    hX=stairs(0:nt,cap_log(target_idx,:)');
    set(hX, 'LineWidth', 1.5)
    legendStrings = [compose("Load %d", target_idx')];
    legend(legendStrings{:});
    title('Capacitors');
    xlabel('k');
    ylabel('Capacitors')

end
% 5. Case
subplot(3,3,5)
img = imread('case.png');
imshow(img);
title('Case considered');

% 6. Voltage input (generators)
subplot(3,3,6)
hold on ;grid on;
hX=stairs(0:(nt-1),u_log');
set(hX, 'LineWidth', 1.5)
hZ=stairs(0:(nt-1),[-u_max*ones(N,nt);u_max*ones(N,nt)]');
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

% 7. Tap
subplot(3,3,7);
hold on ;grid on;
hX = stairs(0:nt+1, tap_log');
set(hX, 'LineWidth', 1.5)
legendHandles = hX(:)';       
legendStrings = [compose("transformer %d-%d (branch %d)", [mpc.branch(transformer_indices,1:2),transformer_indices])];
legend(legendHandles, legendStrings, 'Location', 'northeast');
title('Tap')
xlabel('k')
ylabel('Tap')

% 8. Criterion on V
subplot(3,3,8)
hold on ;grid on;
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

% 9. OLTC commands from QP
subplot(3,3,9);
zoom on; hold on; grid on;
hX = stairs(0:(nt-1), u_oltc_log');
set(hX, 'LineWidth', 1.5)
hEps_pos  = plot([0 nt-1], [ epsilon_oltc  epsilon_oltc], 'k--', 'LineWidth', 1);
hEps_neg  = plot([0 nt-1], [-epsilon_oltc -epsilon_oltc], 'k--', 'LineWidth', 1);
hStep_pos = plot([0 nt-1], [ step_size  step_size], 'k:', 'LineWidth', 1);
hStep_neg = plot([0 nt-1], [-step_size -step_size], 'k:', 'LineWidth', 1);
legendHandles = [hX(:)', hEps_pos, hStep_pos];
oltc_labels = cellstr(compose("OLTC %d (branch %d)", [(1:n_oltc)', transformer_indices]));
legendStrings = [oltc_labels; ...
                 {sprintf('epsilon = %.4f (deadband)', epsilon_oltc)}; ...
                 {sprintf('step size = %.5f (bound)', step_size)}]';
legend(legendHandles, legendStrings, 'Location', 'northeast');
title('OLTC commands from QP');
xlabel('k');
ylabel('u_{OLTC} (p.u.)')


% %% Soutenance finale
% figure;
% 
% % 1. Voltage in targeted loads
% subplot(2,2,1);
% hold on; grid on;
% idx = find(ismember(load_idx, target_idx));
% hX = stairs(0:nt, X_log(idx,:)');
% set(hX, 'LineWidth', 1.5)
% hY = stairs(1:nt, Yc_V(1:n_target,1:nt)', 'x');
% for i =1:n_target
%     hY(i).Color = hX(i).Color;
%     hY(i).LineStyle = 'none';   
%     hY(i).Marker = 'x';        
%     hY(i).LineWidth = 0.8;
% end
% hZ = stairs(1:nt,[repmat(V_load_max(idx),1,nt);repmat(V_load_min(idx),1,nt)]');
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
% legendHandles = [hX(:)',hY(1),hZ(1)];       
% legendStrings = [compose("Load %d", target_idx'),"Target","Constraints" ];
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% title('Voltage in targeted loads');
% xlabel('k');
% ylabel('Voltage (p.u)')
% 
% % 2. Reactive in generators
% subplot(2,2,2);
% hold on; grid on;
% hX = stairs(0:nt,X_log(n_load+1:n_load+n_gen,:)');
% set(hX, 'LineWidth', 1.5)
% hZ = stairs(0:nt,[repmat(Q_gen_max,1,nt+1);repmat(Q_gen_min,1,nt+1)]');
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
% if add_cap_bool
%     hY=stairs(0:nt,[-ones(n_gen,nt+1)*Q_lim;ones(n_gen,nt+1)*Q_lim]');
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
%     legendHandles = [hX(:)',hY(1),hZ(1)];       
%     legendStrings = [compose("Gen %d", pv_idx'),"Targeted limits","Constraints"];
% else
%     legendHandles = [hX(:)',hZ(1)];       
%     legendStrings = [compose("Gen %d", pv_idx'),"Constraints"];
% end
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% title('Reactive power in generators');
% xlabel('k');
% ylabel('Reactive power (p.u)')
% 
% % 4. Voltage in loads or capacitors
% subplot(2,2,4)
% hold on; grid on;
% 
% if not(add_cap_bool)
%     hX = stairs(0:nt, X_log(1:n_load,:)');
%     set(hX, 'LineWidth', 1.5)
%     hZ = stairs(1:nt,[repmat(V_load_max,1,nt);repmat(V_load_min,1,nt)]');
%     for i = 1:n_load
%         hZ(i).Color = hX(i).Color;
%         hZ(i+n_load).Color = hX(i).Color;
%         hZ(i).LineStyle = 'none';  
%         hZ(i+n_load).LineStyle = 'none';   
%         hZ(i).Marker = 'diamond'; 
%         hZ(i+n_load).Marker = 'diamond';        
%         hZ(i).LineWidth = 0.8;
%         hZ(i+n_load).LineWidth = 0.8;
%     end
%     legendHandles = [hX(:)',hZ(1)];       
%     legendStrings = [compose("Load %d", load_idx'),"Constraints" ];
%     legend(legendHandles, legendStrings, 'Location', 'northeast');
%     title('Voltage in loads');
%     xlabel('k');
%     ylabel('Voltage (p.u)')
% else 
%     hX=stairs(0:nt,cap_log(target_idx,:)');
%     set(hX, 'LineWidth', 1.5)
%     legendStrings = [compose("Load %d", target_idx')];
%     legend(legendStrings{:});
%     title('Capacitors');
%     xlabel('k');
%     ylabel('Capacitors')
% 
% end
% 
% % 7. Tap
% subplot(2,2,3);
% hold on ;grid on;
% hX = stairs(0:nt+1, tap_log');
% set(hX, 'LineWidth', 1.5)
% legendHandles = hX(:)';       
% legendStrings = [compose("transformer %d-%d (branch %d)", [mpc.branch(transformer_indices,1:2),transformer_indices])];
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% title('Tap')
% xlabel('k')
% ylabel('Tap')
% 
% %%%%%%%%%%%
% figure;
% 
% % 8. Criterion on V
% subplot(3,1,3)
% hold on ;grid on;
% hX = stairs(0:nt, J_log(3,:)');
% set(hX, 'LineWidth', 1.5);
% x_last = nt;
% y_last = J_log(3,nt);
% plot(x_last, y_last, 'ro','MarkerFaceColor','r');
% text(x_last, y_last, sprintf(' %.3e', y_last), ...
%     'VerticalAlignment','bottom','HorizontalAlignment','left');
% title('Criterion on V');
% xlabel('k'); 
% ylabel('sum (V-V_{target})^2')
% 
% % 8. V_target_secondary
% subplot(3,1,2)
% hold on ;grid on;
% hX = stairs(1:nt, [V_sec_log;X_log(reg_pos,1:end-1)]');
% set(hX, 'LineWidth', 1.5);
% legendHandles = hX(:)';       
% legendStrings = ["Voltage load 5","Voltage load 6","Voltage given to OLTC 10-5","Voltage given to OLTC 11-6"];
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% title('OLTC voltage input');
% xlabel('k'); 
% ylabel('V_target_secondary')
% 
% % 9. OLTC commands from QP
% subplot(3,1,1);
% hold on; grid on;
% hX = stairs(0:(nt-1), u_oltc_log');
% set(hX, 'LineWidth', 1.5)
% hEps_pos  = plot([0 nt-1], [ epsilon_oltc  epsilon_oltc], 'k--', 'LineWidth', 1);
% hEps_neg  = plot([0 nt-1], [-epsilon_oltc -epsilon_oltc], 'k--', 'LineWidth', 1);
% hStep_pos = plot([0 nt-1], [ step_size  step_size], 'k:', 'LineWidth', 1);
% hStep_neg = plot([0 nt-1], [-step_size -step_size], 'k:', 'LineWidth', 1);
% legendHandles = [hX(:)', hEps_pos, hStep_pos];
% oltc_labels = cellstr([compose("transformer %d-%d (branch %d)", [mpc.branch(transformer_indices,1:2),transformer_indices])]);
% legendStrings = [oltc_labels; ...
%                  {sprintf('epsilon = %.4f (deadband)', epsilon_oltc)}; ...
%                  {sprintf('step size = %.5f (bound)', step_size)}]';
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% title('OLTC commands from QP');
% xlabel('k');
% ylabel('u_{OLTC} (p.u.)')

% %% DISPLAY REPORT
% 
% fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'figures', 'JAC');
% if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
% prefix = 'main';
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
% % ---- Figure 3: Tap / U_oltc / Target voltage given to OLTC ----
% fig3 = figure('Units','centimeters','Position',[2 2 16 12]);
% 
% subplot(3,1,1); hold on; grid on;
% hT = stairs(0:nt+1, tap_log');
% set(hT, 'LineWidth', 1.5);
% legend(hT(:)', compose("Transfo %d to %d", mpc.branch(transformer_indices,1:2)), 'Location','northeast');
% title('OLTC tap position'); xlabel('k'); ylabel('Tap ratio (p.u.)');
% 
% subplot(3,1,2);hold on; grid on;
% hX = stairs(0:(nt-1), u_oltc_log');
% set(hX, 'LineWidth', 1.5)
% hEps_pos  = plot([0 nt-1], [ epsilon_oltc  epsilon_oltc], 'k--', 'LineWidth', 1);
% hEps_neg  = plot([0 nt-1], [-epsilon_oltc -epsilon_oltc], 'k--', 'LineWidth', 1);
% hStep_pos = plot([0 nt-1], [ step_size  step_size], 'k:', 'LineWidth', 1);
% hStep_neg = plot([0 nt-1], [-step_size -step_size], 'k:', 'LineWidth', 1);
% legendHandles = [hX(:)', hEps_pos, hStep_pos];
% oltc_labels = cellstr([compose("transformer %d-%d (branch %d)", [mpc.branch(transformer_indices,1:2),transformer_indices])]);
% legendStrings = [oltc_labels; ...
%                  {sprintf('epsilon = %.4f (deadband)', epsilon_oltc)}; ...
%                  {sprintf('step size = %.5f (bound)', step_size)}]';
% legend(legendHandles, legendStrings, 'Location', 'northeast');
% title('OLTC commands from QP');
% xlabel('k');
% ylabel('u_{OLTC} (p.u.)')
% 
% subplot(3,1,3); hold on; grid on;
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