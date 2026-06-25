clear all
close all
%% Parameters definition
% General parameters ------------------------------------------------------
seed = 1234;
nt = 300;   % Number of simulation steps
t_jac = 1; % Jacobian recalculation period: 0= never 1 = every step, n = every n steps
Q_lim_rate = 0.5;  % Soft reactive power limit is Q_lim_rate*QMIN/MAX (p.u.)

% Physical parameters -----------------------------------------------------
u_pv_max   = 0.025;    % Maximum generator voltage increment per step (p.u.)
tap_min   = 0.9;      % Physical tap lower bound
tap_max   = 1.1;      % Physical tap upper bound
tap_step_size = 0.00625;  % Tap step size (p.u.)
tol_tap       = 0.005;    % Voltage deadband for changing tap ratio
t_unavailability_cap = 5;     % Minimum number of steps between two switches on the same bus (cooldown)
t_min_out_Q_lim = 5; % Number of steps our of reactive soft limit required to add a capacitor 
cap           = 5;     % Capacitor/coil size (Mvar); adding a coil is equivalent to cap = -cap
first_delay_oltc=3;   % Number of consecutive out-of-tolerance steps required before the first tap action fires.
                      % Once the threshold is reached, subsequent tap actions fire every step as long as
                      % the secondary voltage remains out of tolerance on the same side with an unchanged setpoint.
                      % The full delay is required again if the voltage returns within tolerance,
                      % crosses to the other side, or the setpoint changes.
% MPC parameters ----------------------------------------------------------
N  = 4;    % MPC prediction horizon (steps)
alpha = 0;  % Weight on Q_gen: voltage tracking dominates (V term has coefficient 1)

% Threshold for converting the continuous QP output u_OLTC to a discrete tap request.
% If |u_OLTC(j)| > epsilon_oltc, a tap movement is requested on transformer j.
epsilon_oltc = 0.5*tap_step_size;

% --- Parameters given to custom_runpf
parameters =[tap_step_size,tol_tap,tap_min,tap_max,cap,first_delay_oltc];

% Scenario parameters -----------------------------------------------------
sigma_rate = 0;  % Rate used to compute the std deviation of noise around the nominal for production and consumption
V_target = [0.98;1;0.99;1.02]; %Voltage demanded by loads

%% Preparation
rng(seed);

% Loading case ------------------------------------------------------------
mpc   = Extented_case();
mpc   = scenario(mpc,0);
mpopt = mpoption('verbose', 0, 'out.all', 0, 'pf.enforce_q_lims', 1);
define_constants;
baseMVA = mpc.baseMVA;

% Index sets --------------------------------------------------------------
pv_idx    = find(mpc.bus(:, BUS_TYPE) == PV);
pq_idx    = find(mpc.bus(:, BUS_TYPE) == PQ);
load_idx  = find(mpc.bus(:, PD) ~= 0);
cap_idx   = [4;6;8;11;12;13;14];
oltc_idx  = [5;6;20;21;22;23]; %branch idx
secondary_oltc_idx = mpc.branch(oltc_idx,T_BUS);
primary_oltc_idx = mpc.branch(oltc_idx,F_BUS);

% Size of index sets ------------------------------------------------------
n_pv=size(pv_idx,1);
n_pq=size(pq_idx,1);
n_load   = size(load_idx,  1);
n_gen = size(mpc.gen,1);
n_cap  = size(cap_idx,1);
n_oltc = size(oltc_idx,1);
n_bus    = size(mpc.bus,1);
n_state = n_pq+n_pv+n_pv+n_oltc;
n_input =n_pv+n_oltc;
n_branch=size(mpc.branch,1);

% Local index in x
load_idx_in_x = zeros(n_load, 1);
for i = 1:n_load
    bus = load_idx(i);
    pq_pos = find(pq_idx == bus);
    pv_pos = find(pv_idx == bus);
    if ~isempty(pq_pos)
        load_idx_in_x(i) = pq_pos;
    else
        load_idx_in_x(i) = n_pq + n_pv + pv_pos;
    end
end
secondary_oltc_idx_in_x = zeros(n_oltc, 1);
for j = 1:n_oltc
    bus = secondary_oltc_idx(j);
    pq_pos = find(pq_idx == bus);
    pv_pos = find(pv_idx == bus);
    if ~isempty(pq_pos)
        secondary_oltc_idx_in_x(j) = pq_pos;
    else
        secondary_oltc_idx_in_x(j) = n_pq + n_pv + pv_pos;
    end
end

% gen_number(i) = row index in mpc.gen of the generator at bus i (0 if none)
gen_number = zeros(n_bus, 1);
for i = 1:n_gen
    gen_number(mpc.gen(i,1)) = i;
end

%Reactive soft limits
Q_lim_up=Q_lim_rate*mpc.gen(gen_number(pv_idx),QMAX)/baseMVA;
Q_lim_down=Q_lim_rate*mpc.gen(gen_number(pv_idx),QMIN)/baseMVA;

% Nominal load and generation values (noise is re-drawn from these each step)
Pd_nom = mpc.bus(:, PD);
Qd_nom = mpc.bus(:, QD);
Pg_nom = mpc.gen(:, PG);

sigma_load_P = Pd_nom*sigma_rate;   % Std of active load noise (MW)
sigma_load_Q = Qd_nom*sigma_rate;   % Std of reactive load noise (Mvar)
sigma_gen_P  = Pg_nom*sigma_rate;   % Std of active generation noise (MW, non-slack generators only)

cap_availability  = zeros(n_bus,  1); %Cooldown counters: cap_availability(j) = 
% steps remaining before capacitor j can request an activation

% Previous OLTC command, stored in the state to model the two-step tap delay.
% Initialized to zero: no tap command was in flight before the simulation starts.
u_oltc_tap_in_flight = zeros(n_oltc, 1);

% Estimated replica of the internal persistence counter of custom_runpf (oltc_counter).
% oltc_counter_est(j) = number of consecutive steps OLTC j has been out of tolerance
% on the same side of its setpoint: positive if voltage is above setpoint, negative if below, 0 if within tolerance.
% It is updated each step using the estimated secondary voltage x_estimated_k1 and the
% setpoint u_oltc_V_k, mirroring exactly the logic in custom_runpf.
% Used to build the delay constraints for the QP and to decide anticipation.
oltc_counter_est=zeros(n_oltc,1);

% Setpoint sent to the OLTC at the previous step (u_oltc_V_k from step k-1).
% Used to detect setpoint changes when updating oltc_counter_est, mirroring
% the same detection done in custom_runpf with its own prev_u_oltc_V.
% Initialised to 1 (in-tolerance neutral value) so that no reset is triggered at k=0.
prev_u_oltc_V=ones(n_oltc,1);

% Voltage value accepted by OLTC
V_max_oltc = mpc.bus(secondary_oltc_idx,VMAX);
V_min_oltc = mpc.bus(secondary_oltc_idx,VMIN);

%% Initial power flow

results = runpf(mpc, mpopt);

% State vector: x = [V_pq (n_pq x 1); Q_pv (n_pv x 1); V_pv (n_pv x 1); U_oltc_tap_prev (n_oltc x 1)]
x = [results.bus(pq_idx, VM); ...
    results.gen(gen_number(pv_idx), QG) / baseMVA; ...
    results.bus(pv_idx, VM);
    u_oltc_tap_in_flight];

% Logging initialization
X_log      = zeros(n_state,  nt+1);   
X_log(:,1) = x;

X_est_log = zeros(n_state,nt+1);
X_est_log(:,1)=x;

U_log  = zeros(n_pv+n_oltc+n_oltc+n_oltc+n_oltc, nt); %u_pv,u_oltc computed by qp,u_oltc discrete,voltage given to oltc,true tap modification

tap_log    = zeros(n_branch, nt+1);     
tap_log(oltc_idx,1) = mpc.branch(oltc_idx, TAP);

cap_log    = zeros(n_bus, nt+1);         

Pd_log     = zeros(n_bus, nt+1);    
Pd_log(load_idx,1)=results.bus(load_idx,PD);

Qd_log     = zeros(n_bus, nt+1);    
Qd_log(load_idx,1)=results.bus(load_idx,QD);

Pg_log     = zeros(n_bus,nt+1);
Pg_log(pv_idx,1)=results.gen(gen_number(pv_idx),PG);

TAP_estimation_error_log  = zeros(n_branch, nt+1);
tap_est_discrete=ones(n_branch,1);

Crit_log=zeros(1,nt+1);                                                                   
Crit_log(1,1)=sum((x(load_idx_in_x)-V_target).^2);   

%% Control loop
for k = 0:nt-1

    %--- Scenario ---
    mpc = scenario(mpc, k);

    %--- Update jacobian and matrices ---
    if mod(k,t_jac)==0 
        % Compute jacobian
        JAC = makeJac(results, true);

        % Extract voltage
        voltage=results.bus(:,VM);
    
        % Compute sensitivity matrices
        [Cv,Cq] = Cal_CvCq(mpc,JAC);
        [Cv_tap,Cq_tap] = Cal_Cv_tapCq_tap(mpc,JAC,voltage,tap_est_discrete);
        [Cv_cap,Cq_cap] = Cal_Cv_capCq_cap(mpc,JAC,voltage);
    
        Cv=Cv(pq_idx,pv_idx);
        Cq=Cq(pv_idx,pv_idx);
        Cv_tap=Cv_tap(pq_idx,oltc_idx);
        Cq_tap=Cq_tap(pv_idx,oltc_idx);

        %Build MPC matrices
        [A,B,H,f1,f2,Yref,A_ineq,b_ineq_x0,b_ineq] = Cal_MPC_matrices(mpc,Cv,Cq,Cv_tap,Cq_tap,N,alpha,u_pv_max,tap_step_size,oltc_idx,V_target);
    end

    % --- Dynamic tap bound constraints ---
    % One row per OLTC per direction: sum of u_oltc(j) over the horizon
    % is bounded by the remaining margin to each physical limit.

    margin_up   = tap_max - tap_est_discrete(oltc_idx);   % (n_oltc x 1)
    margin_down = tap_est_discrete(oltc_idx) - tap_min;   % (n_oltc x 1)
    
    % S_oltc: (n_oltc x n_input*N), entry (j, (i-1)*n_input + n_pv + j) = 1
    % Selects and sums u_oltc(j) across all N horizon steps
    S_oltc = repmat([zeros(n_oltc, n_pv), eye(n_oltc)], 1, N);
    
    A_ineq_dyn = [A_ineq;  S_oltc; -S_oltc];
    b_ineq_dyn = [ b_ineq_x0 * x + b_ineq; margin_up; margin_down];

    % --- Dynamic OLTC delay constraints ---
    % The physical OLTC will not act until its internal counter reaches first_delay_oltc.
    % The MPC must not plan a tap action before that threshold is reachable, otherwise
    % it would schedule a movement that custom_runpf will silently ignore.
    %
    % steps_to_fire(j): number of additional out-of-tolerance steps still needed for
    % OLTC j before its counter reaches first_delay_oltc in the current direction.
    % During these steps the corresponding u_oltc(j) is blocked to zero.
    % Once the counter has reached the threshold (steps_to_fire == 0), the normal
    % tap amplitude constraints from Cal_MPC_matrices apply and the OLTC can act
    % every step as long as the voltage remains out of tolerance on the same side.

    steps_to_fire = max(0, first_delay_oltc - abs(oltc_counter_est));  % (n_oltc x 1)

    % Build one blocking row per OLTC per blocked step:
    A_ineq_oltc_delay=[];
    b_ineq_oltc_delay=[];
    for j=1:n_oltc
        if oltc_counter_est(j) > 0
            % Counter is positive: voltage is currently above setpoint, counter is accumulating
            % toward a tap increase. steps_to_fire(j) steps are still needed before the
            % counter reaches first_delay_oltc, so tap increases are blocked for that many steps.
            % Tap decreases are blocked for the full first_delay_oltc steps: the counter
            % would have to reset to 0 first and then re-accumulate in the other direction.

            % Block tap increase for the remaining steps before threshold
            A_rows_under_construction=zeros(steps_to_fire(j)-1,N*n_input);
            cols_selected = n_input*(1:steps_to_fire(j)-1) - n_oltc + j;
            A_rows_under_construction(:,cols_selected)=eye(steps_to_fire(j)-1);
            A_ineq_oltc_delay=[A_ineq_oltc_delay;A_rows_under_construction];
            b_ineq_oltc_delay=[b_ineq_oltc_delay;zeros(steps_to_fire(j)-1,1)];
            % Block tap decrease for the full delay (counter must reset and re-accumulate)
            A_rows_under_construction=zeros(first_delay_oltc-1,N*n_input);
            cols_selected = n_input*(1:first_delay_oltc-1) - n_oltc + j;
            A_rows_under_construction(:,cols_selected)=-eye(first_delay_oltc-1);
            A_ineq_oltc_delay=[A_ineq_oltc_delay;A_rows_under_construction];
            b_ineq_oltc_delay=[b_ineq_oltc_delay;zeros(first_delay_oltc-1,1)];
        elseif oltc_counter_est(j) < 0
            % Counter is negative: voltage is currently below setpoint, counter is accumulating
            % toward a tap decrease. Symmetric to the positive case.

            % Block tap increase for the full delay (counter must reset and re-accumulate)
            A_rows_under_construction=zeros(first_delay_oltc-1,N*n_input);
            cols_selected = n_input*(1:first_delay_oltc-1) - n_oltc + j;
            A_rows_under_construction(:,cols_selected)=eye(first_delay_oltc-1);
            A_ineq_oltc_delay=[A_ineq_oltc_delay;A_rows_under_construction];
            b_ineq_oltc_delay=[b_ineq_oltc_delay;zeros(first_delay_oltc-1,1)];
            % Block tap decrease for the remaining steps before threshold
            A_rows_under_construction=zeros(steps_to_fire(j)-1,N*n_input);
            cols_selected = n_input*(1:steps_to_fire(j)-1) - n_oltc + j;
            A_rows_under_construction(:,cols_selected)=-eye(steps_to_fire(j)-1);
            A_ineq_oltc_delay=[A_ineq_oltc_delay;A_rows_under_construction];
            b_ineq_oltc_delay=[b_ineq_oltc_delay;zeros(steps_to_fire(j)-1,1)];
        else
            % Counter is zero: voltage is within tolerance (or just reset).
            % The full delay first_delay_oltc must elapse before any tap action in either direction.

            % Block tap increase for the full delay
            A_rows_under_construction=zeros(first_delay_oltc-1,N*n_input);
            cols_selected = n_input*(1:first_delay_oltc-1) - n_oltc + j;
            A_rows_under_construction(:,cols_selected)=eye(first_delay_oltc-1);
            A_ineq_oltc_delay=[A_ineq_oltc_delay;A_rows_under_construction];
            b_ineq_oltc_delay=[b_ineq_oltc_delay;zeros(first_delay_oltc-1,1)];
            % Block tap decrease for the full delay
            A_rows_under_construction=zeros(first_delay_oltc-1,N*n_input);
            cols_selected = n_input*(1:first_delay_oltc-1) - n_oltc + j;
            A_rows_under_construction(:,cols_selected)=-eye(first_delay_oltc-1);
            A_ineq_oltc_delay=[A_ineq_oltc_delay;A_rows_under_construction];
            b_ineq_oltc_delay=[b_ineq_oltc_delay;zeros(first_delay_oltc-1,1)];
        end
    end
    A_ineq_dyn = [A_ineq_dyn; A_ineq_oltc_delay];
    b_ineq_dyn = [b_ineq_dyn; b_ineq_oltc_delay];

    % --- Priority 1: MPC (generator voltage + OLTC) ---
    [U, ~,exitflag] = quadprog(H, f1*x - f2*Yref, A_ineq_dyn,b_ineq_dyn, [], [], [], [], [], optimset('Display', 'off'));
    if exitflag <= 0
        warning('QP failed at step %d with flag %d — applying zero input', k,exitflag);
        U = zeros(n_input * N, 1);
    end

    % Receding horizon: apply only the first step of the optimal sequence
    u_pv_k        = U(1:n_pv);
    u_oltc_tap_k  = U(n_pv+1:n_pv+n_oltc);

    U_log(1:n_input,k+1)=[u_pv_k;u_oltc_tap_k];

    % Convert the continuous QP output u_oltc_tap_k to a discrete tap request.
    % Values within the deadband [-epsilon_oltc, +epsilon_oltc] are treated as no action.
    % Values outside the deadband are mapped to exactly one tap step (+/-tap_step_size).
    u_oltc_tap_discrete_k=zeros(n_oltc,1);
    for j = 1:n_oltc
        if u_oltc_tap_k(j) > epsilon_oltc
            u_oltc_tap_discrete_k(j)=tap_step_size;
        elseif u_oltc_tap_k(j) < -epsilon_oltc
            u_oltc_tap_discrete_k(j)=-tap_step_size;
        end
    end
    
    U_log(n_pv+n_oltc+1:n_pv+n_oltc+n_oltc,k+1)=u_oltc_tap_discrete_k;

    % --- Priority 2: Capacitor/coil decision ---

    % State estimation
    u_k           = [u_pv_k;u_oltc_tap_discrete_k];
    x_estimated_k1 = A*x + B*u_k;
    
    % Capacitor logic
    if k>=t_min_out_Q_lim 
        last_Q = [X_log(n_pq+1:n_pq+n_pv, k+2-t_min_out_Q_lim+1:k+1), x_estimated_k1(n_pq+1:n_pq+n_pv,:)]; 
        is_out_up = (sum(last_Q>Q_lim_up,2)==t_min_out_Q_lim);
        is_out_down = (sum(last_Q<Q_lim_down,2)==t_min_out_Q_lim);
        if any(is_out_down) || any(is_out_up)
            disp(['k=', num2str(k), ' | any_out_up=', num2str(any(is_out_up)), ' | any_out_down=', num2str(any(is_out_down))]);
            u_cap_k=add_cap(mpc,x_estimated_k1,Cv_cap,Cq_cap,cap_availability,cap_idx,cap,Q_lim_rate,cap_log(:,k+1));
        else
            u_cap_k=zeros(n_cap,1);
        end    
    else
        u_cap_k=zeros(n_cap,1);
    end

    cap_availability(cap_idx)=max(0,cap_availability(cap_idx)+abs(u_cap_k)*(t_unavailability_cap+1)-1);
    x_estimated_k1(1:n_pq)=x_estimated_k1(1:n_pq)+Cv_cap(pq_idx,cap_idx)*u_cap_k*cap/baseMVA;%Update x estimated
    x_estimated_k1(n_pq+1:n_pq+n_pv)=x_estimated_k1(n_pq+1:n_pq+n_pv)+Cq_cap(pv_idx,cap_idx)*u_cap_k*cap/baseMVA;
    
    % --- Compute voltage setpoint sent to the OLTC (u_oltc_V_k) ---
    %
    % custom_runpf does not receive a direct tap command. Instead it receives a
    % voltage setpoint and decides internally whether to move the tap based on
    % how long the secondary voltage has been out of tolerance on the same side.
    %
    % Mapping from discrete tap request to voltage setpoint:
    %   tap increase requested  (+tap_step_size) -> setpoint = V_min_oltc
    %                                               (pushes secondary above tolerance, counter accumulates upward)
    %   tap decrease requested  (-tap_step_size) -> setpoint = V_max_oltc
    %                                               (pushes secondary below tolerance, counter accumulates downward)
    %   no tap action requested (0)              -> setpoint = estimated secondary voltage
    %                                               (secondary appears in tolerance, counter resets to 0)
    u_oltc_V_k = x_estimated_k1(secondary_oltc_idx_in_x);
    for j = 1:n_oltc
        if u_oltc_tap_discrete_k(j) == tap_step_size %Increase tap
            u_oltc_V_k(j) = V_min_oltc(j); %Decrease voltage
        elseif u_oltc_tap_discrete_k(j) == -tap_step_size %Decrease tap
            u_oltc_V_k(j) = V_max_oltc(j); %Increase voltage
        end
    end
    
    % --- Anticipation: start accumulating the OLTC counter before the planned tap step ---
    %
    % The QP may plan a tap action at horizon step t_fire_j > 1, but the physical OLTC
    % will only act once its internal counter has been out of tolerance for first_delay_oltc
    % consecutive steps. If the counter has not yet accumulated enough steps, the tap action
    % planned by the QP will not fire on time.
    %
    % Anticipation: if the number of out-of-tolerance steps still needed (steps_needed_j)
    % is >= t_fire_j, push the OLTC setpoint outside tolerance NOW (even though no immediate
    % tap action was requested) so the counter reaches the threshold exactly when t_fire_j arrives.
    % This may override the neutral setpoint set in the block above for OLTCs with no immediate action.
    %
    % Credit is given for counter steps already accumulated in the correct direction,
    % so anticipation is skipped when the counter is already on track.
    for j = 1:n_oltc
        % Extract planned tap commands for OLTC j over the full horizon
        u_oltc_j_horizon = U((0:N-1)*n_input + n_pv + j);        % Zero out values below the deadband to avoid triggering on QP noise
        u_oltc_j_horizon(abs(u_oltc_j_horizon) < epsilon_oltc) = 0;
        % First horizon step at which a tap action is planned
        t_fire_j = find(u_oltc_j_horizon, 1, 'first');
        if ~isempty(t_fire_j)
            dir_j = sign(u_oltc_j_horizon(t_fire_j));   % +1: tap increase, -1: tap decrease
            % Steps still needed to reach the threshold in direction dir_j.
            % max(0, ...) gives credit for counter steps already in the right direction.
            steps_needed_j = first_delay_oltc - max(0, dir_j * oltc_counter_est(j));
            if steps_needed_j >= t_fire_j
                %In matpower tap=HV/LV
                %Increase LV => Decrease tap
                if dir_j == 1
                    u_oltc_V_k(j) = V_min_oltc(j);    % anticipate tap increase
                elseif dir_j == -1
                    u_oltc_V_k(j) = V_max_oltc(j);    % anticipate tap decrease
                end
            end
        end
    end

    % Update oltc_counter_est to mirror the internal persistence counter of custom_runpf.
    % Uses the estimated secondary voltages at k+1 (x_estimated_k1) and the setpoint
    % u_oltc_V_k just computed. The update rules are identical to those in custom_runpf:
    %   - Setpoint change         -> reset to 0 (full delay required again)
    %   - Out of tolerance, same side -> increment counter in that direction
    %   - Out of tolerance, side change -> restart from +/-1 (full delay required in new direction)
    %   - Within tolerance        -> reset to 0
    v_secondary_est = x_estimated_k1(secondary_oltc_idx_in_x);
    diffs_est = v_secondary_est-u_oltc_V_k ;
    %In matpower TAP=HV/LV 
    %Boost LV => decrease tap

    for j = 1:n_oltc
        % Setpoint change: reset counter so the full delay is required again.
        % If the voltage is simultaneously out of tolerance, the increment below
        % will run in the same iteration and bring the counter to +/-1.
        % This matches exactly the behaviour of custom_runpf (same two-block structure).
        if u_oltc_V_k(j) ~= prev_u_oltc_V(j)
            oltc_counter_est(j) = 0;
        end
        if abs(diffs_est(j)) > tol_tap
            new_sign = sign(diffs_est(j));
            if oltc_counter_est(j) == 0 || sign(oltc_counter_est(j)) == new_sign
                % Voltage out of tolerance on the same side: increment counter.
                % Or voltage out of tolerance for the first time, intialize
                % counter
                oltc_counter_est(j) = oltc_counter_est(j) + new_sign;
            else
                % Voltage switched side: restart count from +/-1
                oltc_counter_est(j) = new_sign;
            end
        else
            % Voltage back within tolerance: reset
            oltc_counter_est(j) = 0;
        end
    end

    prev_u_oltc_V = u_oltc_V_k;
    
    U_log(n_pv+n_oltc+n_oltc+1:n_pv+n_oltc+n_oltc+n_oltc,k+1)=u_oltc_V_k;

    % --- Apply input and run powerflow ---

    [mpc,results] = custom_runpf(mpc,u_pv_k,u_oltc_V_k,u_cap_k,parameters,oltc_idx,cap_idx);
    %mpc: mpc used next step with tap modified
    %results: power flow of mpc from the current step with u_pv_k, u_cap_k
    %applied and before modifying the tap (delayed)
    
    %--- Read measurements ---
    V_pq=results.bus(pq_idx, VM);
    Q_pv=results.gen(gen_number(pv_idx), QG) / baseMVA;
    V_pv=results.bus(pv_idx, VM);
    V_secondary_oltc_idx=results.bus(secondary_oltc_idx,VM);

    % Reconstruct the tap movement that custom_runpf actually applied at this step,
    % using the measured secondary voltage and the setpoint sent (u_oltc_V_k).
    % This may differ from u_oltc_tap_discrete_k when x_estimated and the measured
    % state diverge: the counter in custom_runpf is driven by the true secondary
    % voltage, not the estimated one, so the physical tap may fire at a different
    % step than predicted.
    % The measured counter state is approximated by oltc_counter_est: a tap action
    % is confirmed only if the counter has reached the threshold in the right direction.
    u_oltc_tap_in_flight=zeros(n_oltc,1);
    for j=1:n_oltc
        if V_secondary_oltc_idx(j)-u_oltc_V_k(j)>tol_tap && oltc_counter_est(j)>=first_delay_oltc
            u_oltc_tap_in_flight(j) = tap_step_size;
        elseif V_secondary_oltc_idx(j)-u_oltc_V_k(j)<-tol_tap && oltc_counter_est(j)<=-first_delay_oltc
            u_oltc_tap_in_flight(j) = -tap_step_size;
        end
    end
    
    U_log(n_pv+n_oltc+n_oltc+n_oltc+1:end,k+1)=u_oltc_tap_in_flight;

    % Read new state
    x = [results.bus(pq_idx, VM); ...
        results.gen(gen_number(pv_idx), QG) / baseMVA; ...
        results.bus(pv_idx, VM);
        u_oltc_tap_in_flight];

    %--- Tap estimation ---
    R_oltc  = results.branch(oltc_idx, BR_R);
    X_oltc  = results.branch(oltc_idx, BR_X);
    B_oltc  = results.branch(oltc_idx, BR_B);
    V_F = results.bus(primary_oltc_idx,   VM);
    V_T = results.bus(secondary_oltc_idx, VM);
    PF_pu   = results.branch(oltc_idx, PF) / baseMVA;
    QF_pu   = results.branch(oltc_idx, QF) / baseMVA;

    % Remove the shunt magnetizing current contribution from the measured flows
    PF_series = PF_pu;                             % B has no effect on active power (pure susceptance)
    QF_series = QF_pu - (B_oltc/2) .* V_F.^2; % Reactive flow through the series impedance only

    tap_estimated = V_F ./ (V_T + (R_oltc .* PF_series + X_oltc .* QF_series) ./ V_T);

    tap_est_discrete(oltc_idx) = round((tap_estimated - tap_min) / tap_step_size) * tap_step_size + tap_min;
    tap_est_discrete = min(max(tap_est_discrete, tap_min), tap_max);

    tap_measured = results.branch(oltc_idx, TAP);

    % --- Logging ---
    X_log(:,k+2)=x;
    X_est_log(:,k+2)=x_estimated_k1;
    tap_log(oltc_idx,k+2)    = results.branch(oltc_idx, TAP);
    cap_log(cap_idx ,k+2)    = cap_log(cap_idx,k+1) + u_cap_k;      
    Pd_log(load_idx,k+2)     = results.bus(load_idx,PD);
    Qd_log(load_idx,k+2)     = results.bus(load_idx,QD);
    Pg_log(pv_idx,k+2)       = results.gen(gen_number(pv_idx),PG);
    TAP_estimation_error_log(oltc_idx,k+2)=tap_measured-tap_est_discrete(oltc_idx);
    Crit_log(1,k+2)=sum((x(load_idx_in_x)-V_target).^2);

       
end

%% DISPLAY

close all

%% Summary 
figure;
% 1. Voltage in loads
subplot(3,3,1);
hold on; grid on;
hX = stairs(0:nt, X_log(load_idx_in_x,:)');
set(hX, 'LineWidth', 1.5)
hZ = stairs(0:nt,[repmat(mpc.bus(load_idx,VMAX),1,nt+1);repmat(mpc.bus(load_idx,VMIN),1,nt+1)]');
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
hT = stairs(0:nt, repmat(V_target, 1, nt+1)');
for i = 1:n_load
    hT(i).Color = hX(i).Color;
    hT(i).LineStyle = '--';
    hT(i).LineWidth = 0.8;
end
legendHandles = [hX(:)', hT(1), hZ(1)];
legendStrings = [compose("Load %d", load_idx'), "V target", "Constraints"];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Voltage in loads');
xlabel('k');
ylabel('Voltage (p.u)')

% 2. Reactive in pv generators
subplot(3,3,2);
hold on; grid on;
hX = stairs(0:nt,X_log(n_pq+1:n_pq+n_pv,:)');
set(hX, 'LineWidth', 1.5)
hZ = stairs(0:nt,[repmat(mpc.gen(gen_number(pv_idx),QMAX),1,nt+1)/baseMVA;repmat(mpc.gen(gen_number(pv_idx),QMIN),1,nt+1)/baseMVA]');
for i = 1:n_pv
    hZ(i).Color = hX(i).Color;
    hZ(i+n_pv).Color = hX(i).Color;
    hZ(i).LineStyle = 'none';  
    hZ(i+n_pv).LineStyle = 'none';   
    hZ(i).Marker = 'diamond'; 
    hZ(i+n_pv).Marker = 'diamond';        
    hZ(i).LineWidth = 0.8;
    hZ(i+n_pv).LineWidth = 0.8;
end

hY = stairs(0:nt,repmat([Q_lim_up;Q_lim_down],1,nt+1)');
for i = 1:n_pv
    hY(i).Color = hX(i).Color;
    hY(i+n_pv).Color = hX(i).Color;
    hY(i).LineStyle = 'none';  
    hY(i+n_pv).LineStyle = 'none';   
    hY(i).Marker = 'x'; 
    hY(i+n_pv).Marker = 'x';        
end

legendHandles = [hX(:)',hY(1),hZ(1)];       
legendStrings = [compose("Gen %d", pv_idx'),"Targeted limits","Constraints"];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Reactive power in generators');
xlabel('k');
ylabel('Reactive power (p.u)')

% 3. Voltage in generators
subplot(3,3,3);
hold on; grid on;
hX = stairs(0:nt,X_log(n_pq+n_pv+1:n_pq+n_pv+n_pv,:)');
set(hX, 'LineWidth', 1.5)
hZ = stairs(0:nt,[repmat(mpc.bus(pv_idx,VMAX),1,nt+1);repmat(mpc.bus(pv_idx,VMIN),1,nt+1)]');
for i = 1:n_pv
    hZ(i).Color = hX(i).Color;
    hZ(i+n_pv).Color = hX(i).Color;
    hZ(i).LineStyle = 'none';  
    hZ(i+n_pv).LineStyle = 'none';   
    hZ(i).Marker = 'diamond'; 
    hZ(i+n_pv).Marker = 'diamond';        
    hZ(i).LineWidth = 0.8;
    hZ(i+n_pv).LineWidth = 0.8;
end
legendHandles = [hX(:)',hZ(1)];       
legendStrings = [compose("Gen %d", pv_idx'),"Constraints" ];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Voltage in generators');
xlabel('k');
ylabel('Voltage (p.u)')

% 4. Capacitors
subplot(3,3,4)
hold on; grid on;

hX=stairs(0:nt,cap_log(cap_idx,:)');
set(hX, 'LineWidth', 1.5)
legendStrings = [compose("Bus %d", cap_idx')];
lg = legend(legendStrings{:});
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Capacitors');
xlabel('k');
ylabel('Capacitors')

% 5. Criterion
subplot(3,3,5)

hold on; grid on;
h = stairs(0:nt, Crit_log');
set(h, 'LineWidth', 1.5)

legendHandles = h(:)';
legendStrings = {'Criterion'};
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Criterion')
xlabel('k')
ylabel('Criterion')

% 6. Voltage input (generators)
subplot(3,3,6)
hold on ;grid on;
hX=stairs(0:(nt-1),U_log(1:n_pv,:)');
set(hX, 'LineWidth', 1.5)

hZ=stairs(0:(nt-1),[-u_pv_max*ones(1,nt);u_pv_max*ones(1,nt)]');
hZ(1).Marker = 'diamond'; 
hZ(2).Marker = 'diamond';        
hZ(1).LineWidth = 0.8;
hZ(2).LineWidth = 0.8;
hZ(2).Color=hZ(1).Color;
hZ(1).LineStyle = 'none';  
hZ(2).LineStyle = 'none'; 

legendHandles = [hX(:)',hZ(1)];       
legendStrings = [compose("Gen %d", pv_idx'),"Constraints" ];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Voltage input');
xlabel('k');
ylabel('Voltage input (p.u)')

% 7. Tap
subplot(3,3,7);
hold on ;grid on;
hX = stairs(0:nt, tap_log(oltc_idx,:)');
set(hX, 'LineWidth', 1.5)

hZ=stairs(0:nt,[tap_min*ones(1,nt+1);tap_max*ones(1,nt+1)]');
hZ(1).Marker = 'diamond'; 
hZ(2).Marker = 'diamond';        
hZ(1).LineWidth = 2;
hZ(2).LineWidth = 2;
hZ(2).Color=hZ(1).Color;
hZ(1).LineStyle = 'none';  
hZ(2).LineStyle = 'none'; 

legendHandles = [hX(:)',hZ(1)];       
legendStrings = [compose("OLTC %d-%d (branch %d)", [mpc.branch(oltc_idx,F_BUS:T_BUS),oltc_idx])',"Limits"];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Tap')
xlabel('k')
ylabel('Tap')

% 8. OLTC commands from QP
subplot(3,3,8);
hold on; grid on;
hX = stairs(0:(nt-1), U_log(n_pv+1:n_pv+n_oltc,:)');
set(hX, 'LineWidth', 1.5)
hEps_pos  = plot([0 nt-1], [ epsilon_oltc  epsilon_oltc], 'k--', 'LineWidth', 1);
hEps_neg  = plot([0 nt-1], [-epsilon_oltc -epsilon_oltc], 'k--', 'LineWidth', 1);
hStep_pos = plot([0 nt-1], [ tap_step_size  tap_step_size], 'k:', 'LineWidth', 1);
hStep_neg = plot([0 nt-1], [-tap_step_size -tap_step_size], 'k:', 'LineWidth', 1);
legendHandles = [hX(:)', hEps_pos, hStep_pos];
legendStrings = [compose("OLTC %d-%d (branch %d)", [mpc.branch(oltc_idx,F_BUS:T_BUS), oltc_idx])', ...
    sprintf('epsilon = %.4f (deadband)', epsilon_oltc), ...
    sprintf('step size = %.5f (bound)', tap_step_size)];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('OLTC commands from QP');
xlabel('k');
ylabel('u_{OLTC} (p.u.)')

% 9. Voltage given to OLTC
subplot(3,3,9);
hold on; grid on;
hX = stairs(0:(nt-1), U_log(n_pv+n_oltc+n_oltc+1:n_pv+n_oltc+n_oltc+n_oltc,:)');
set(hX, 'LineWidth', 1.5)

legendHandles = hX(:)';
legendStrings = compose("OLTC %d-%d (branch %d)", [mpc.branch(oltc_idx,F_BUS:T_BUS),oltc_idx])';
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Voltage given to OLTC')
xlabel('k')
ylabel('Voltage (p.u)')
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10)

%% ESTIMATION
figure;

% TAP estimation
subplot(3,1,1)
hold on; grid on;

hErr = stairs(0:nt, TAP_estimation_error_log(oltc_idx,:)');
set(hErr, 'LineWidth', 1.5);

% Horizontal lines for all possible tap errors (difference between any two tap values)
tap_values = tap_min:tap_step_size:tap_max;
for tv = tap_values
    yline(tv - tap_values, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end
yline(0, 'k-', 'LineWidth', 1.2);

colors = get(gca, 'ColorOrder');
for j = 1:n_oltc
    c = colors(mod(j-1, size(colors,1)) + 1, :);
    hErr(j).Color = c;
end

legendHandles = hErr(:)';
legendStrings = compose("OLTC %d-%d (branch %d)", [mpc.branch(oltc_idx,F_BUS:T_BUS), oltc_idx])';
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);

title('Tap estimation error');
xlabel('k');
ylabel('tap_{est} - tap_{measured} (p.u.)');

% Difference between TAP modification requested and real
subplot(3,1,2);
hold on; grid on;
h = stairs(0:nt-1, ...
    U_log(n_pv+n_oltc+1:n_pv+n_oltc+n_oltc,:)' ...
    - U_log(n_pv+n_oltc+n_oltc+n_oltc+1:end, :)');
set(h, 'LineWidth', 1.5)

legendHandles = h(:)';
legendStrings = compose("OLTC %d-%d (branch %d)", [mpc.branch(oltc_idx,F_BUS:T_BUS), oltc_idx])';
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Difference between TAP modification requested and real')
xlabel('k')
ylabel('Tap modification')

% OLTC secondary voltage estimated vs measured
subplot(3,1,3);
hold on; grid on;
hMeas = stairs(0:nt, X_log(secondary_oltc_idx_in_x,:)');
set(hMeas, 'LineWidth', 1.5)

% Estimated secondary voltage
hEst = stairs(0:nt, X_est_log(secondary_oltc_idx_in_x,:)');
for j = 1:n_oltc
    hEst(j).Color     = hMeas(j).Color;
    hEst(j).LineStyle = '--';
    hEst(j).LineWidth = 1.2;
end

legendHandles = [hMeas(:)', hEst(1)];
legendStrings = [compose("OLTC %d-%d (branch %d)", [mpc.branch(oltc_idx,F_BUS:T_BUS), oltc_idx])',...
    sprintf('V estimated (secondary)')];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Estimated secondary voltage vs measured')
xlabel('k')
ylabel('Voltage (p.u)')
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10)

%% Information relative to tap

figure;

% OLTC commands from QP
subplot(3,1,1);
hold on; grid on;
hX = stairs(0:(nt-1), U_log(n_pv+1:n_pv+n_oltc,:)');
set(hX, 'LineWidth', 1.5)
hEps_pos  = plot([0 nt-1], [ epsilon_oltc  epsilon_oltc], 'k--', 'LineWidth', 1);
hEps_neg  = plot([0 nt-1], [-epsilon_oltc -epsilon_oltc], 'k--', 'LineWidth', 1);
hStep_pos = plot([0 nt-1], [ tap_step_size  tap_step_size], 'k:', 'LineWidth', 1);
hStep_neg = plot([0 nt-1], [-tap_step_size -tap_step_size], 'k:', 'LineWidth', 1);
legendHandles = [hX(:)', hEps_pos, hStep_pos];
legendStrings = [compose("OLTC %d-%d (branch %d)", [mpc.branch(oltc_idx,F_BUS:T_BUS), oltc_idx])', ...
    sprintf('epsilon = %.4f (deadband)', epsilon_oltc), ...
    sprintf('step size = %.5f (bound)', tap_step_size)];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('OLTC commands from QP');
xlabel('k');
ylabel('u_{OLTC} (p.u.)')

% Voltage given to OLTC vs estimated
subplot(3,1,2);
hold on; grid on;
hX = stairs(0:(nt-1), U_log(n_pv+n_oltc+n_oltc+1:n_pv+n_oltc+n_oltc+n_oltc,:)');
set(hX, 'LineWidth', 1.5)

% Estimated secondary voltage (used to build u_oltc_V_k, logged one step ahead)
hEst = stairs(0:(nt-1), X_est_log(secondary_oltc_idx_in_x,2:nt+1)');
for j = 1:n_oltc
    hEst(j).Color     = hX(j).Color;
    hEst(j).LineStyle = '--';
    hEst(j).LineWidth = 1.2;
end

legendHandles = [hX(:)', hEst(1)];
legendStrings = [compose("OLTC %d-%d (branch %d)", [mpc.branch(oltc_idx,F_BUS:T_BUS), oltc_idx])',...
    sprintf('V estimated (secondary)')];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Voltage given to OLTC vs estimated secondary voltage')
xlabel('k')
ylabel('Voltage (p.u)')

% Tap
subplot(3,1,3);
hold on ;grid on;
hX = stairs(0:nt, tap_log(oltc_idx,:)');
set(hX, 'LineWidth', 1.5)

hZ=stairs(0:nt,[tap_min*ones(1,nt+1);tap_max*ones(1,nt+1)]');
hZ(1).Marker = 'diamond'; 
hZ(2).Marker = 'diamond';        
hZ(1).LineWidth = 2;
hZ(2).LineWidth = 2;
hZ(2).Color=hZ(1).Color;
hZ(1).LineStyle = 'none';  
hZ(2).LineStyle = 'none'; 

legendHandles = [hX(:)',hZ(1)];       
legendStrings = [compose("OLTC %d-%d (branch %d)", [mpc.branch(oltc_idx,F_BUS:T_BUS),oltc_idx])',"Limits"];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Tap')
xlabel('k')
ylabel('Tap')
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10)

%% Voltage and reactive power

figure;
%Voltage in loads
subplot(2,2,1);
hold on; grid on;
hX = stairs(0:nt, X_log(load_idx_in_x,:)');
set(hX, 'LineWidth', 1.5)
hZ = stairs(0:nt,[repmat(mpc.bus(load_idx,VMAX),1,nt+1);repmat(mpc.bus(load_idx,VMIN),1,nt+1)]');
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
hT = stairs(0:nt, repmat(V_target, 1, nt+1)');
for i = 1:n_load
    hT(i).Color = hX(i).Color;
    hT(i).LineStyle = '--';
    hT(i).LineWidth = 0.8;
end
legendHandles = [hX(:)', hT(1), hZ(1)];
legendStrings = [compose("Load %d", load_idx'), "V target", "Constraints"];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Voltage in loads');
xlabel('k');
ylabel('Voltage (p.u)')

% Reactive in pv generators
subplot(2,2,2);
hold on; grid on;
hX = stairs(0:nt,X_log(n_pq+1:n_pq+n_pv,:)');
set(hX, 'LineWidth', 1.5)
hZ = stairs(0:nt,[repmat(mpc.gen(gen_number(pv_idx),QMAX),1,nt+1)/baseMVA;repmat(mpc.gen(gen_number(pv_idx),QMIN),1,nt+1)/baseMVA]');
for i = 1:n_pv
    hZ(i).Color = hX(i).Color;
    hZ(i+n_pv).Color = hX(i).Color;
    hZ(i).LineStyle = 'none';  
    hZ(i+n_pv).LineStyle = 'none';   
    hZ(i).Marker = 'diamond'; 
    hZ(i+n_pv).Marker = 'diamond';        
    hZ(i).LineWidth = 0.8;
    hZ(i+n_pv).LineWidth = 0.8;
end

hY = stairs(0:nt,repmat([Q_lim_up;Q_lim_down],1,nt+1)');
for i = 1:n_pv
    hY(i).Color = hX(i).Color;
    hY(i+n_pv).Color = hX(i).Color;
    hY(i).LineStyle = 'none';  
    hY(i+n_pv).LineStyle = 'none';   
    hY(i).Marker = 'x'; 
    hY(i+n_pv).Marker = 'x';        
end

legendHandles = [hX(:)',hY(1),hZ(1)];       
legendStrings = [compose("Gen %d", pv_idx'),"Targeted limits","Constraints"];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Reactive power in generators');
xlabel('k');
ylabel('Reactive power (p.u)')

% Voltage input (generators)
subplot(2,2,3)
hold on ;grid on;
hX=stairs(0:(nt-1),U_log(1:n_pv,:)');
set(hX, 'LineWidth', 1.5)

hZ=stairs(0:(nt-1),[-u_pv_max*ones(1,nt);u_pv_max*ones(1,nt)]');
hZ(1).Marker = 'diamond'; 
hZ(2).Marker = 'diamond';        
hZ(1).LineWidth = 0.8;
hZ(2).LineWidth = 0.8;
hZ(2).Color=hZ(1).Color;
hZ(1).LineStyle = 'none';  
hZ(2).LineStyle = 'none'; 

legendHandles = [hX(:)',hZ(1)];       
legendStrings = [compose("Gen %d", pv_idx'),"Constraints" ];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Voltage input');
xlabel('k');
ylabel('Voltage input (p.u)')

% Voltage in generators
subplot(2,2,4);
hold on; grid on;
hX = stairs(0:nt,X_log(n_pq+n_pv+1:n_pq+n_pv+n_pv,:)');
set(hX, 'LineWidth', 1.5)
hZ = stairs(0:nt,[repmat(mpc.bus(pv_idx,VMAX),1,nt+1);repmat(mpc.bus(pv_idx,VMIN),1,nt+1)]');
for i = 1:n_pv
    hZ(i).Color = hX(i).Color;
    hZ(i+n_pv).Color = hX(i).Color;
    hZ(i).LineStyle = 'none';  
    hZ(i+n_pv).LineStyle = 'none';   
    hZ(i).Marker = 'diamond'; 
    hZ(i+n_pv).Marker = 'diamond';        
    hZ(i).LineWidth = 0.8;
    hZ(i+n_pv).LineWidth = 0.8;
end
legendHandles = [hX(:)',hZ(1)];       
legendStrings = [compose("Gen %d", pv_idx'),"Constraints" ];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Voltage in generators');
xlabel('k');
ylabel('Voltage (p.u)')
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10)

%% Scenario

figure;
% 1. Active load per bus
subplot(3,1,1);
hold on; grid on;
hX = stairs(0:nt, Pd_log(load_idx,:)');
set(hX, 'LineWidth', 1.5)
legendStrings = compose("Bus %d", load_idx');
lg = legend(legendStrings{:}, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Active load');
xlabel('k');
ylabel('P (MW)')

% 2. Reactive load per bus
subplot(3,1,2);
hold on; grid on;
hX = stairs(0:nt, Qd_log(load_idx,:)');
set(hX, 'LineWidth', 1.5)
legendStrings = compose("Bus %d", load_idx');
lg = legend(legendStrings{:}, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Reactive load');
xlabel('k');
ylabel('Q (Mvar)')

% 3. Active generation per PV bus
subplot(3,1,3);
hold on; grid on;
hX = stairs(0:nt, Pg_log(pv_idx,:)');
set(hX, 'LineWidth', 1.5)
hZ = stairs(0:nt, [repmat(mpc.gen(gen_number(pv_idx),PMAX),1,nt+1); ...
    repmat(mpc.gen(gen_number(pv_idx),PMIN),1,nt+1)]');
for i = 1:n_pv
    hZ(i).Color        = hX(i).Color;
    hZ(i+n_pv).Color   = hX(i).Color;
    hZ(i).LineStyle    = 'none';
    hZ(i+n_pv).LineStyle = 'none';
    hZ(i).Marker       = 'diamond';
    hZ(i+n_pv).Marker  = 'diamond';
    hZ(i).LineWidth    = 0.8;
    hZ(i+n_pv).LineWidth = 0.8;
end
legendHandles = [hX(:)', hZ(1)];
legendStrings = [compose("Gen %d", pv_idx'), "Constraints"];
lg = legend(legendHandles, legendStrings, 'Location', 'northeast');
lg.ItemHitFcn = @(~,e) set(e.Peer, Visible=~e.Peer.Visible);
title('Active generation (PV buses)');
xlabel('k');
ylabel('P (MW)')
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10)