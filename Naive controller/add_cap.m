function [mpc,cap_availability,cap_log] = add_cap(mpc,JAC,x_estimated_k1,x_estimated_k2,cap_availability,cap_log,t_cap,var_v_cap_last,t_cap_last,var_v_cap_next,t_cap_next,cap,k,X_log,Yc_V,Q_lim)
% add_cap  Selects and applies the best single capacitor/coil switching action (lowest priority actuator).
%
% Called after action_oltc (priority 2). Modifies mpc.bus (shunt susceptance) before the
% power flow, so the action takes effect immediately at k+1 (no delay).
%
% A bus is eligible only if:
%   (1) its cooldown has expired (cap_availability == 0)
%   (2) the load voltage at that bus has been stable over the last t_cap_last steps
%   (3) the voltage reference for that bus will not change over the next t_cap_next steps
% These conditions prevent switching during transients (same logic as action_oltc).
%
% Each candidate action (add capacitor or add coil on one eligible bus) is evaluated
% at both k+1 (via x_estimated_k1) and k+2 (via x_estimated_k2, already corrected for
% the MPC and OLTC actions). Hard constraints (VMAX/VMIN, QMAX/QMIN) are checked at
% both horizons. The criterion minimised is sum(Q_gen_out^2), where Q_gen_out are the
% reactive powers of generators currently outside [-Q_lim, Q_lim].
%
% Inputs:
%   mpc              : current MATPOWER case (modified in place if a switch is made)
%   JAC              : full network Jacobian (shared with all sensitivity functions)
%   x_estimated_k1   : predicted state at k+1 [V_load; Q_gen; V_gen]
%   x_estimated_k2   : predicted state at k+2 [V_load; Q_gen; V_gen],
%                      already corrected for MPC and OLTC effects by main.m
%   cap_availability : (n_bus x 1) cooldown counters; bus i is blocked if > 0
%   cap_log          : (n_bus x (k+1)) cumulative capacitor count per bus (for logging)
%   t_cap            : cooldown duration (steps) after each switch on the same bus
%   var_v_cap_last   : max allowed voltage deviation (p.u.) at the target bus over last t_cap_last steps
%   t_cap_last       : number of past steps used for the voltage stability check
%   var_v_cap_next   : max allowed reference variation (p.u.) over next t_cap_next steps
%   t_cap_next       : number of future steps used for the reference stability check
%   cap              : shunt susceptance increment (Mvar); coil = -cap
%   k                : current step (0-based)
%   X_log            : (n_load x (k+1)) logged state history; columns are x(0)..x(k)
%   Yc_V             : (n_target x T) voltage reference trajectory
%   Q_lim            : soft reactive power limit (p.u.); criterion targets Q_gen in [-Q_lim, Q_lim]
%
% Outputs:
%   mpc              : updated MATPOWER case (mpc.bus(:,BS) modified if a switch occurred)
%   cap_availability : updated cooldown counters (all decremented; acting bus reset to t_cap)
%   cap_log          : updated log (new column appended)


%% Preparation
    define_constants;

    baseMVA=mpc.baseMVA;

    pv_idx    = find(mpc.bus(:,BUS_TYPE)==PV);  % controlled generators
    load_idx  = find(mpc.bus(:,BUS_TYPE)==PQ);  % all load buses (including transit buses)
    target_idx= find(mpc.bus(:,PD)~=0);         % buses with active load (switching targets)
    slack_idx = find(mpc.bus(:,BUS_TYPE)==REF);

    n_load = size(load_idx,1);
    n_gen  = size(pv_idx,1);
    n_slack= size(slack_idx,1);
    n_bus  = n_load+n_gen+n_slack;

    gen_number=zeros(n_bus,1); % gen_number(i) = row in mpc.gen for the generator at bus i
    for i =1:n_gen+n_slack
        gen_number(mpc.gen(i,1))=i;
    end

    % Hard constraints (from network data)
    V_load_max=mpc.bus(load_idx,VMAX);
    V_load_min=mpc.bus(load_idx,VMIN);
    Q_gen_max =mpc.gen(pv_idx,QMAX)/baseMVA;
    Q_gen_min =mpc.gen(pv_idx,QMIN)/baseMVA;

    % Sensitivity matrices evaluated at each prediction horizon.
    % Cv_cap(i,j): dV_load(i) / d(one capacitor added on bus j)
    % Cq_cap(i,j): dQ_gen(i)  / d(one capacitor added on bus j)
    % The voltage at the capacitor bus appears in the injection formula Q = B*V^2,
    % so the sensitivities depend on the operating point and differ at k+1 and k+2.
    [Cv_cap_k1,Cq_cap_k1] = Cal_Cv_capCq_cap(mpc,x_estimated_k1,cap,JAC);
    [Cv_cap_k2,Cq_cap_k2] = Cal_Cv_capCq_cap(mpc,x_estimated_k2,cap,JAC);
    
%% Eligibility filter
    % possible_loads(i) = 1 if bus i is a valid switching candidate, 0 otherwise.
    % Buses are eligible only if they are target buses (active load), their cooldown
    % has expired, the voltage has been stable recently, and the reference is stable ahead.

    possible_loads = zeros(n_bus,1);
    possible_loads(target_idx) = 1; % only buses with active load are candidates

    % Condition 1 — cooldown: bus i is blocked if cap_availability(i) > 0
    possible_loads = max(possible_loads - cap_availability, zeros(n_bus,1));

    % Condition 2 — past stability: voltage at each bus must have stayed within
    % var_v_cap_last of the current estimate over the last t_cap_last steps
    last_moments = zeros(n_bus, t_cap_last);
    last_moments(load_idx,:) = X_log(1:n_load, k-t_cap_last+2:k+1);
    V_load_estimated_k1 = zeros(n_bus,1);
    V_load_estimated_k1(load_idx) = x_estimated_k1(1:n_load);
    % voltage_last_moments(i,j) = 1 if bus i was within tolerance at step k-t_cap_last+1+j
    voltage_last_moments = abs(last_moments - V_load_estimated_k1) < var_v_cap_last;
    load_last_moments = sum(voltage_last_moments, 2) == t_cap_last; % 1 if stable over all t_cap_last steps
    possible_loads = min(possible_loads, load_last_moments);

    % Condition 3 — future reference stability: the voltage reference must not change
    % more than var_v_cap_next over the next t_cap_next steps
    next_moments = zeros(n_bus, t_cap_next);
    next_moments(target_idx,:) = Yc_V(:, k+2:k+t_cap_next+1); % reference at k+2..k+t_cap_next+1
    current_objective = zeros(n_bus,1);
    current_objective(target_idx) = Yc_V(:, k+1);             % reference at k+1 (current target)
    voltage_next_moments = abs(next_moments - current_objective) < var_v_cap_next;
    load_next_moments = sum(voltage_next_moments, 2) == t_cap_next;
    possible_loads = min(possible_loads, load_next_moments);

    possible_loads_idx  = find(possible_loads);        % bus indices of eligible buses
    n_possible_loads    = size(possible_loads_idx, 1);

    if n_possible_loads==0
        % No eligible bus: only decrement cooldown counters, log unchanged state
        cap_availability = max(cap_availability - ones(n_bus,1), zeros(n_bus,1));
        if isempty(cap_log)
            cap_log = zeros(n_bus,1);
        else
            cap_log = [cap_log, cap_log(:,end)];
        end
    else
%% Enumerate candidate actions and select the best
        Q_gen_estimated_k1 = x_estimated_k1(n_load+1:n_load+n_gen);
        V_load_estimated_k1= x_estimated_k1(1:n_load);
        Q_gen_estimated_k2 = x_estimated_k2(n_load+1:n_load+n_gen);
        V_load_estimated_k2= x_estimated_k2(1:n_load);

        % Baseline criterion: sum of Q^2 for generators currently outside [-Q_lim, Q_lim] at k+2
        gen_out     = Q_gen_estimated_k2 > Q_lim | Q_gen_estimated_k2 < -Q_lim;
        current_crit= sum(Q_gen_estimated_k2(gen_out).^2);

        % criterion_mat: one row per candidate action
        % Columns: [bus_number, action (+1=capacitor / -1=coil), criterion]
        criterion_mat = zeros(2*n_possible_loads, 3);

        for i = 1:n_possible_loads
            bus_number = possible_loads_idx(i);
            criterion_mat(2*i-1, 1:2) = [bus_number,  1]; % capacitor
            criterion_mat(2*i,   1:2) = [bus_number, -1]; % coil

            for action = [1, -1]
                row = 2*i - (action==1); % row 2*i-1 for cap, 2*i for coil

                % Predict V_load and Q_gen at both horizons if this action is applied
                % (sign: +1 for capacitor injects Q, -1 for coil absorbs Q)
                new_Q_gen_k1 = Q_gen_estimated_k1 + action * Cq_cap_k1(pv_idx,  bus_number);
                new_V_load_k1= V_load_estimated_k1 + action * Cv_cap_k1(load_idx, bus_number);
                new_Q_gen_k2 = Q_gen_estimated_k2 + action * Cq_cap_k2(pv_idx,  bus_number);
                new_V_load_k2= V_load_estimated_k2 + action * Cv_cap_k2(load_idx, bus_number);

                % Reject if hard constraints violated at either k+1 or k+2
                feasible = all(new_V_load_k1 < V_load_max) && all(new_V_load_k1 > V_load_min) && ...
                           all(new_Q_gen_k1  < Q_gen_max)  && all(new_Q_gen_k1  > Q_gen_min)  && ...
                           all(new_V_load_k2 < V_load_max) && all(new_V_load_k2 > V_load_min) && ...
                           all(new_Q_gen_k2  < Q_gen_max)  && all(new_Q_gen_k2  > Q_gen_min);

                if feasible
                    new_gen_out = new_Q_gen_k2 > Q_lim | new_Q_gen_k2 < -Q_lim;
                    criterion_mat(row, 3) = sum(new_Q_gen_k2(new_gen_out).^2);
                else
                    criterion_mat(row, 3) = inf;
                end
            end
        end

        % Select best action — only if it strictly improves the baseline criterion
        [best_new_crit, idx] = min(criterion_mat(:,3));
        if best_new_crit < current_crit
            load_and_action = criterion_mat(idx, 1:2); % [bus_number, +1 or -1]
        else
            load_and_action = [];
        end

%% Apply the selected action
        cap_availability = max(cap_availability - ones(n_bus,1), zeros(n_bus,1)); % decrement all cooldowns

        if isempty(load_and_action)
            % No improvement found: log unchanged state
            if isempty(cap_log)
                cap_log = zeros(n_bus,1);
            else
                cap_log = [cap_log, cap_log(:,end)];
            end
        else
            % Apply shunt susceptance change (takes effect at the next power flow, i.e. k+1)
            mpc.bus(load_and_action(1), BS) = mpc.bus(load_and_action(1), BS) + load_and_action(2)*cap;
            cap_availability(load_and_action(1)) = t_cap; % reset cooldown for this bus
            if isempty(cap_log)
                cap_log = zeros(n_bus,1);
                cap_log(load_and_action(1)) = load_and_action(2);
            else
                cap_log = [cap_log, cap_log(:,end)];
                cap_log(load_and_action(1), end) = cap_log(load_and_action(1), end) + load_and_action(2);
            end
        end
    end
end