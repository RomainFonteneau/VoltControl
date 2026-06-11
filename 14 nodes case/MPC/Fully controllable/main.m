clear all

%% Parameters definition
% General parameters ------------------------------------------------------
seed = 1234;
nt = 150;   % Number of simulation steps
t_jac = 0; % Jacobian recalculation period: 0= never 1 = every step, n = every n steps
Q_lim_rate = 0.6;  % Soft reactive power limit is Q_lim_rate*QMIN/MAX (p.u.)

% Physical parameters -----------------------------------------------------
u_pv_max   = 0.05;    % Maximum generator voltage increment per step (p.u.)
tap_min   = 0.9;      % Physical tap lower bound
tap_max   = 1.1;      % Physical tap upper bound
tap_step_size = 0.00625;  % Tap step size (p.u.)
tol_tap       = 0.005;    % Voltage deadband for changing tap ratio
t_unavailability_cap = 5;     % Minimum number of steps between two switches on the same bus (cooldown)
t_min_out_Q_lim = 5; % Number of steps our of reactive soft limit required to add a capacitor 
cap           = 5;     % Capacitor/coil size (Mvar); adding a coil is equivalent to cap = -cap

% MPC parameters ----------------------------------------------------------
N  = 3;    % MPC prediction horizon (steps)
alpha = 0;  % Weight on Q_gen: voltage tracking dominates (V term has coefficient 1)

% Threshold for converting the continuous QP output u_OLTC to a discrete tap request.
% If |u_OLTC(j)| > epsilon_oltc, a tap movement is requested on transformer j.
epsilon_oltc = tap_step_size / 2;

% Scenario parameters -----------------------------------------------------
sigma_rate = 0;  % Rate used to compute the std deviation of noise around the nominal for production and consumption
V_target = [1;1;1;1]; %Voltage demanded by loads

do_open_line=false; 
line_opened_idx = [17;18]; % Idx of opened line
step_opening_line = 25; % Step of opening

do_close_line=false;
step_closing_line = 90; % Step of closing 

do_increase_PG=false; %Increase active production
PG_increased_idx=[10]; 
PG_increase_per_step=2;
step_start_PG_increase=40;
step_end_PG_increase=50;

do_decrease_PG=false; %Decrease active production
PG_decreased_idx=[10]; 
PG_decrease_per_step=3;
step_start_PG_decrease=40;
step_end_PG_decrease=50;

do_increase_PD=false; %Increase active consumption
PD_increased_idx=[14]; 
PD_increase_per_step=3;
step_start_PD_increase=40;
step_end_PD_increase=50;

do_decrease_PD=false; %Decrease active consumption
PD_decreased_idx=[14]; 
PD_decrease_per_step=3;
step_start_PD_decrease=40;
step_end_PD_decrease=50;

do_increase_QD=false; %Increase reactive consumption
QD_increased_idx=[14]; 
QD_increase_per_step=1;
step_start_QD_increase=40;
step_end_QD_increase=50;

do_decrease_QD=false; %Decrease reactive consumption
QD_decreased_idx=[14]; 
QD_decrease_per_step=1;
step_start_QD_decrease=40;
step_end_QD_decrease=50;

% Parameters given to custom_runpf
parameters =[tap_step_size,tol_tap,tap_min,tap_max,cap];

%% Preparation
rng(seed);

% Loading case ------------------------------------------------------------
mpc   = Extented_case();
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
Q_lim_up=Q_lim_rate*mpc.gen(gen_number(pv_idx),QMAX);
Q_lim_down=Q_lim_rate*mpc.gen(gen_number(pv_idx),QMIN);

% Nominal load and generation values (noise is re-drawn from these each step)
Pd_nom = mpc.bus(:, PD);
Qd_nom = mpc.bus(:, QD);
Pg_nom = mpc.gen(:, PG);

sigma_load_P = Pd_nom*sigma_rate;   % Std of active load noise (MW)
sigma_load_Q = Qd_nom*sigma_rate;   % Std of reactive load noise (Mvar)
sigma_gen_P  = Pg_nom*sigma_rate;   % Std of active generation noise (MW, non-slack generators only)

cap_availability  = zeros(n_bus,  1); %Cooldown counters: cap_availability(j) = steps remaining before capacitor j
% can request an activation

% Previous OLTC command, stored in the state to model the two-step tap delay.
% Initialized to zero: no tap command was in flight before the simulation starts.
u_oltc_tap_in_flight = zeros(n_oltc, 1);

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

U_log      = zeros(n_pv+n_oltc+n_oltc, nt); 

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

%% Control loop
for k = 0:nt-1

    %--- Scenario ---
    if do_open_line && k==step_opening_line
        mpc.branch(line_opened_idx,BR_STATUS)=0;
    end

    if do_close_line && k==step_closing_line
        mpc.branch(line_opened_idx,BR_STATUS)=1;
    end

    if do_increase_PG && k>=step_start_PG_increase && k<=step_end_PG_increase
        mpc.gen(gen_number(PG_increased_idx),PG)=min(mpc.gen(gen_number(PG_increased_idx),PMAX),mpc.gen(gen_number(PG_increased_idx),PG)+PG_increase_per_step);
    end

    if do_decrease_PG && k>=step_start_PG_decrease && k<=step_end_PG_decrease
        mpc.gen(gen_number(PG_decreased_idx),PG)=max(mpc.gen(gen_number(PG_decreased_idx),PMIN),mpc.gen(gen_number(PG_decreased_idx),PG)+PG_decrease_per_step);
    end

    if do_increase_PD && k>=step_start_PD_increase && k<=step_end_PD_increase
        mpc.bus(PD_increased_idx,PD)=mpc.bus(PD_increased_idx,PD)+PD_increase_per_step;
    end

    if do_decrease_PD && k>=step_start_PD_decrease && k<=step_end_PD_decrease
        mpc.bus(PD_decreased_idx,PD)=mpc.bus(PD_decreased_idx,PD)+PD_decrease_per_step;
    end

    if do_increase_QD && k>=step_start_QD_increase && k<=step_end_QD_increase
        mpc.bus(QD_increased_idx,QD)=mpc.bus(QD_increased_idx,QD)+QD_increase_per_step;
    end

    if do_decrease_QD && k>=step_start_QD_decrease && k<=step_end_QD_decrease
        mpc.bus(QD_decreased_idx,QD)=mpc.bus(QD_decreased_idx,QD)+QD_decrease_per_step;
    end

    
    %--- Update jacobian and matrices ---
    if mod(k,t_jac)==0 
        % Compute jacobian
        JAC = makeJac(results, true);

        % Extract voltage
        voltage=results.bus(:,VM);
    
        % Compute sensitivity matrices
        [Cv,Cq] = Cal_CvCq(mpc,JAC);
        [Cv_tap,Cq_tap] = Cal_Cv_tapCq_tap(mpc,JAC,voltage);
        [Cv_cap,Cq_cap] = Cal_Cv_capCq_cap(mpc,JAC,voltage);
    
        Cv=Cv(pq_idx,pv_idx);
        Cq=Cq(pv_idx,pv_idx);
        Cv_tap=Cv_tap(pq_idx,oltc_idx);
        Cq_tap=Cq_tap(pv_idx,oltc_idx);

        %Build MPC matrices
        [A,B,H,f1,f2,Yref,A_ineq,b_ineq_x0,b_ineq] = Cal_MPC_matrices(mpc,Cv,Cq,Cv_tap,Cq_tap,N,alpha,u_pv_max,tap_step_size,oltc_idx,V_target);
    end
    
    %--- Tap estimation ---
    R_oltc  = mpc.branch(oltc_idx, BR_R);
    X_oltc  = mpc.branch(oltc_idx, BR_X);
    B_oltc  = mpc.branch(oltc_idx, BR_B);
    V_F = results.bus(primary_oltc_idx,   VM);
    V_T = results.bus(secondary_oltc_idx, VM);
    PF_pu   = results.branch(oltc_idx, PF) / baseMVA;
    QF_pu   = results.branch(oltc_idx, QF) / baseMVA;

    % Remove the shunt magnetizing current contribution from the measured flows
    PF_series = PF_pu;                             % B has no effect on active power (pure susceptance)
    QF_series = QF_pu - (B_oltc/2) .* V_F.^2; % Reactive flow through the series impedance only

    tap_estimated = V_F ./ (V_T + (R_oltc .* PF_series + X_oltc .* QF_series) ./ V_T);

    tap_est_discrete = round((tap_estimated - tap_min) / tap_step_size) * tap_step_size + tap_min;
    tap_est_discrete = min(max(tap_est_discrete, tap_min), tap_max);

    tap_measured = results.branch(oltc_idx, TAP);

    TAP_estimation_error_log(oltc_idx,k+2)=tap_measured-tap_est_discrete;

    % --- Dynamic tap bound constraints ---
    % One row per OLTC per direction: sum of u_oltc(j) over the horizon
    % is bounded by the remaining margin to each physical limit.

    margin_up   = tap_max - tap_est_discrete;   % (n_oltc x 1)
    margin_down = tap_est_discrete - tap_min;   % (n_oltc x 1)
    
    % S_oltc: (n_oltc x n_input*N), entry (j, (i-1)*n_input + n_pv + j) = 1
    % Selects and sums u_oltc(j) across all N horizon steps
    S_oltc = repmat([zeros(n_oltc, n_pv), eye(n_oltc)], 1, N);
    
    A_ineq_dyn = [A_ineq;  S_oltc; -S_oltc];
    b_ineq_dyn = [ b_ineq_x0 * x + b_ineq; margin_up; margin_down];

    % --- Priority 1: MPC (generator voltage + OLTC) ---
    [U, ~,exitflag] = quadprog(H, f1*x - f2*Yref, A_ineq_dyn,b_ineq_dyn, [], [], [], [], [], optimset('Display', 'off'));
    if exitflag <= 0
        warning('QP failed at step %d with flag %d — applying zero input', k,exitflag);
        U = zeros(n_input * N, 1);
    end

    % Receding horizon: apply only the first step of the optimal sequence
    u_pv_k        = U(1:n_pv);
    u_oltc_tap_k  = U(n_pv+1:n_input);

    U_log(1:n_input,k+1)=[u_pv_k;u_oltc_tap_k];

    % Convert continuous u_oltc_tap_k to discrete value
    u_oltc_tap_discrete_k=zeros(n_oltc,1);
    for j = 1:n_oltc
        if u_oltc_tap_k(j) > epsilon_oltc
            u_oltc_tap_discrete_k(j)=tap_step_size;
        elseif u_oltc_tap_k(j) < -epsilon_oltc
            u_oltc_tap_discrete_k(j)=-tap_step_size;
        end
    end

    % --- Priority 2: Capacitor/coil decision ---

    % State estimation
    u_k           = [u_pv_k;u_oltc_tap_discrete_k];
    x_estimated_k1 = A*x + B*u_k;

    %u_cap_k = add_cap();
    if k>t_min_out_Q_lim 
        last_Q=[X_log(n_pq+1:n_pq+n_pv,end-t_min_out_Q_lim+2:end),x_estimated_k1(n_pq+1:n_pq+n_pv,:)];
        is_out_up = (sum(last_Q>Q_lim_up,2)==t_min_out_Q_lim);
        is_out_down = (sum(last_Q<Q_lim_down,2)==t_min_out_Q_lim);
        if any(is_out_down) || any(is_out_up)
            u_cap_k=add_cap(mpc,x_estimated_k1,Cv_cap,Cq_cap,cap_availability,cap_idx,cap,Q_lim_rate,cap_log(:,k+1));
        else
            u_cap_k=zeros(n_cap,1);
        end    
    else
        u_cap_k=zeros(n_cap,1);
    end

    cap_availability(cap_idx)=max(0,cap_availability(cap_idx)+u_cap_k*(t_unavailability_cap+1)-1);
    
    % --- Apply input and run power flow ---

    % Convert discrete u_oltc_tap_k to u_oltc_V_k (voltage given to the OLTC)
    u_oltc_V_k=x_estimated_k1(secondary_oltc_idx_in_x);
    for j = 1:n_oltc
        if u_oltc_tap_discrete_k(j)==tap_step_size
            u_oltc_V_k(j) = u_oltc_V_k(j) - 2*tol_tap;   % tap increase requested
        elseif u_oltc_tap_discrete_k(j)==-tap_step_size
            u_oltc_V_k(j) = u_oltc_V_k(j) + 2*tol_tap;   % tap decrease requested
        end
    end

    U_log(n_input+1:end,k+1)=u_oltc_V_k;

    [mpc,results] = custom_runpf(mpc,u_pv_k,u_oltc_V_k,u_cap_k,parameters,oltc_idx,cap_idx);
    %mpc: mpc used next step with tap modified
    %results: power flow of mpc from the current step with u_pv_k, u_cap_k
    %applied and before modifying the tap (delayed)
    
    %--- Read measurements ---
    V_pq=results.bus(pq_idx, VM);
    Q_pv=results.gen(gen_number(pv_idx), QG) / baseMVA;
    V_pv=results.bus(pv_idx, VM);
    V_secondary_oltc_idx=results.bus(secondary_oltc_idx,VM);

    %Compute OLTC action truly in flight, if x_estimated and x_measured are
    %different the tap movement wanted could be different than the tap
    %movement truly in flight
    u_oltc_tap_in_flight=zeros(n_oltc,1);
    for j=1:n_oltc
        % V2(LV) = V1(HV)/ratio
        if V_secondary_oltc_idx(j)-u_oltc_V_k(j)>tol_tap %Voltage needs to be decreased
            u_oltc_tap_in_flight(j) = tap_step_size; % tap increase in flight
        elseif V_secondary_oltc_idx(j)-u_oltc_V_k(j)<-tol_tap %Voltage needs to be increased
            u_oltc_tap_in_flight(j) = -tap_step_size; % tap decrease in flight
        end
    end

    % Read new state
    x = [results.bus(pq_idx, VM); ...
        results.gen(gen_number(pv_idx), QG) / baseMVA; ...
        results.bus(pv_idx, VM);
        u_oltc_tap_in_flight];

    % --- Logging ---
    X_log(:,k+2)=x;
    tap_log(oltc_idx,k+2)    = results.branch(oltc_idx, TAP);
    cap_log(cap_idx ,k+2)    = cap_log(cap_idx,k+1) + u_cap_k;      
    Pd_log(load_idx,k+2)     = results.bus(load_idx,PD);
    Qd_log(load_idx,k+2)     = results.bus(load_idx,QD);
    Pg_log(pv_idx,k+2)       = results.gen(gen_number(pv_idx),PG);
       
end
%% DISPLAY

close all

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
legendHandles = [hX(:)',hZ(1)];       
legendStrings = [compose("Load %d", load_idx'),"Constraints" ];
legend(legendHandles, legendStrings, 'Location', 'northeast');
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

hY = stairs(0:nt,repmat([mpc.gen(gen_number(pv_idx),QMAX);mpc.gen(gen_number(pv_idx),QMIN)],1,nt+1)'*Q_lim_rate/baseMVA);
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
legend(legendHandles, legendStrings, 'Location', 'northeast');
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
legend(legendHandles, legendStrings, 'Location', 'northeast');
title('Voltage in generators');
xlabel('k');
ylabel('Voltage (p.u)')

% 4. Capacitors
subplot(3,3,4)
hold on; grid on;

hX=stairs(0:nt,cap_log(cap_idx,:)');
set(hX, 'LineWidth', 1.5)
legendStrings = [compose("Bus %d", cap_idx')];
legend(legendStrings{:});
title('Capacitors');
xlabel('k');
ylabel('Capacitors')

% % 5. Case
% ax = subplot(3,3,5);
% img = imread('Extended_case.png');
% imshow(img, 'Parent', ax);
% title(ax, 'Case considered');

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
legend(legendHandles, legendStrings, 'Location', 'northeast');
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
legend(legendHandles, legendStrings, 'Location', 'northeast');
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
legend(legendHandles, legendStrings, 'Location', 'northeast');
title('OLTC commands from QP');
xlabel('k');
ylabel('u_{OLTC} (p.u.)')

% 9. Voltage given to OLTC
subplot(3,3,9);
hold on; grid on;
hX = stairs(0:(nt-1), U_log(n_pv+n_oltc+1:n_pv+n_oltc+n_oltc,:)');
set(hX, 'LineWidth', 1.5)

legendHandles = hX(:)';
legendStrings = compose("OLTC %d-%d (branch %d)", [mpc.branch(oltc_idx,F_BUS:T_BUS),oltc_idx])';
legend(legendHandles, legendStrings, 'Location', 'northeast');
title('Voltage given to OLTC')
xlabel('k')
ylabel('Voltage (p.u)')


figure;
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
legend(legendHandles, legendStrings, 'Location', 'northeast');

title('Tap estimation error');
xlabel('k');
ylabel('tap_{est} - tap_{measured} (p.u.)');