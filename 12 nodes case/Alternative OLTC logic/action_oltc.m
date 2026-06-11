function [oltc_idx,direction,target_v,x_estimated_k2,oltc_availability] = action_oltc(mpc,JAC,x_estimated_k2, transformer_indices, step_size, tol, alpha_tap, Q_lim, Yc_V, k, tap_min, tap_max, oltc_availability, t_oltc, var_v_oltc_last, t_oltc_last, var_v_oltc_next, t_oltc_next, X_log);
% action_oltc  Selects the best single OLTC tap action at step k (intermediate priority actuator).
%
% Called after the MPC (priority 1) and before the capacitor logic (priority 3).
% Evaluates up to 2*n_oltc candidate actions (+1 or -1 step on each transformer).
% The action minimising the criterion
%     J = sum((V_load_target - V_ref)^2) + alpha_tap * sum(Q_gen_out^2)
% is selected if it strictly improves over doing nothing AND passes all hard constraints.
% At most ONE transformer moves per call.
%
% Tap delay: run_pf_oltc_step advances the tap AFTER the power flow, so a tap decided
% at step k first affects the power flow at k+1, and x(k+1) was computed without it —
% the effective delay is one step. All evaluations are therefore done on x_estimated_k2
% (state predicted at k+2 by main.m). If an action is selected, x_estimated_k2 is
% updated with the predicted tap effect so the capacitor logic receives a consistent state.
%
% Inputs:
%   mpc                 : current MATPOWER case (tap state from previous step)
%   JAC                 : full network Jacobian (shared with all sensitivity functions)
%   x_estimated_k2      : predicted state at k+2 [V_load; Q_gen; V_gen], already
%                         corrected for MPC action in main.m
%   transformer_indices : branch indices of OLTC transformers (n_oltc x 1)
%   step_size           : tap increment per action (p.u.); must match run_pf_oltc_step
%   tol                 : voltage deadband for direction encoding; must match run_pf_oltc_step
%   alpha_tap           : weight on reactive power in the OLTC criterion:
%                         J = ||V_load_target - V_ref||^2 + alpha_tap*||Q_gen_out||^2
%   Q_lim               : soft reactive power limit (p.u.); generators beyond [-Q_lim, Q_lim]
%                         contribute to the criterion
%   Yc_V                : (n_target x T) voltage reference trajectory; column k+2 is used
%   k                   : current step (0-based)
%   tap_min, tap_max    : physical tap bounds; must match run_pf_oltc_step ([0.9, 1.1])
%   oltc_availability   : (n_oltc x 1) cooldown counters; transformer j is blocked if > 0
%   t_oltc              : cooldown duration (steps) after each tap action
%   var_v_oltc_last     : max allowed voltage deviation (p.u.) at secondary bus over last t_oltc_last steps
%   t_oltc_last         : number of past steps checked for voltage stability
%   var_v_oltc_next     : max allowed reference variation (p.u.) over next t_oltc_next steps
%   t_oltc_next         : number of future steps checked for reference stability
%   X_log               : (n_load x (k+1)) logged state history; columns are x(0)..x(k)
%
% Outputs:
%   oltc_idx            : local index (1..n_oltc) of the transformer that moved, 0 if none
%   direction           : tap direction applied (-1 = decrease, +1 = increase, 0 = none)
%   target_v            : (n_oltc x 1) target voltages for run_pf_oltc_step
%                         Encodes the tap direction via the deadband logic:
%                           target = V_reg         -> no movement (stays inside deadband)
%                           target = V_reg + 2*tol -> tap decrease (V_LV too high)
%                           target = V_reg - 2*tol -> tap increase (V_LV too low)
%   x_estimated_k2      : updated predicted state at k+2, corrected for the tap effect if one was applied
%   oltc_availability   : updated cooldown counters

    define_constants;

    baseMVA=mpc.baseMVA;

    load_idx   = find(mpc.bus(:, BUS_TYPE) == PQ);
    pv_idx     = find(mpc.bus(:, BUS_TYPE) == PV);
    target_idx = find(mpc.bus(:, PD) ~= 0); % buses with active load (voltage tracking targets)

    n_load   = length(load_idx);
    n_gen    = length(pv_idx);
    n_oltc   = length(transformer_indices);

    % Hard constraints (from network data)
    V_load_max = mpc.bus(load_idx, VMAX);
    V_load_min = mpc.bus(load_idx, VMIN);
    Q_gen_max  = mpc.gen(pv_idx, QMAX) / baseMVA;
    Q_gen_min  = mpc.gen(pv_idx, QMIN) / baseMVA;

    % Predicted state at k+2 (passed in from main.m, already includes MPC effect)
    V_load_cur = x_estimated_k2(1:n_load);
    Q_gen_cur  = x_estimated_k2(n_load+1 : n_load+n_gen);

    % Position of target buses within load_idx (used to extract the voltage tracking subset)
    target_in_load = find(ismember(load_idx, target_idx));

    % Secondary (LV) buses of each OLTC transformer, and their position within load_idx
    % reg_pos(j): row index in V_load_cur of transformer j's secondary bus
    reg_buses    = mpc.branch(transformer_indices, T_BUS);
    [~, reg_pos] = ismember(reg_buses, load_idx);
    V_reg        = V_load_cur(reg_pos); % current predicted voltage at each secondary bus

    % Sensitivity matrices: dV_load/d_tap and dQ_gen/d_tap, evaluated at x_estimated_k2
    [Cv_tap,Cq_tap] = Cal_Cv_tapCq_tap(mpc, x_estimated_k2, transformer_indices,JAC);

    %% Availability and stability checks
    % A transformer is eligible to act only if:
    %   (1) its cooldown has expired (oltc_availability == 0)
    %   (2) its secondary bus voltage has been stable over the last t_oltc_last steps
    %   (3) the voltage reference for its secondary bus will not change significantly
    %       over the next t_oltc_next steps
    % Conditions (2) and (3) prevent tap hunting during transients.

    oltc_availability = max(oltc_availability - 1, zeros(n_oltc, 1)); % decrement cooldown counters

    % For each transformer, find the position of its secondary bus in target_idx
    % (oltc_in_target(j) = 0 if the secondary bus has no voltage reference)
    oltc_in_target = zeros(n_oltc, 1);
    for j = 1:n_oltc
        idx = find(target_idx == reg_buses(j));
        if ~isempty(idx)
            oltc_in_target(j) = idx;
        end
    end

    eligible = false(n_oltc, 1);
    for j = 1:n_oltc
        if oltc_availability(j) > 0, continue; end % cooldown not expired

        % Past stability: secondary bus voltage must not have deviated more than
        % var_v_oltc_last from the current estimate over the last t_oltc_last steps
        if k >= t_oltc_last
            recent_V_j  = X_log(reg_pos(j), end-t_oltc_last+1:end);
            past_stable = all(abs(recent_V_j - V_load_cur(reg_pos(j))) < var_v_oltc_last);
        else
            past_stable = false; % not enough history yet
        end
        if ~past_stable, continue; end

        % Future stability: voltage reference must not change more than var_v_oltc_next
        % over the next t_oltc_next steps (only checked if the bus has a reference)
        if oltc_in_target(j) > 0
            future_refs   = Yc_V(oltc_in_target(j), k+2 : k+t_oltc_next+1);
            future_stable = all(abs(future_refs - Yc_V(oltc_in_target(j), k+2)) < var_v_oltc_next);
        else
            future_stable = true;
        end
        if ~future_stable, continue; end

        eligible(j) = true;
    end

    %% Baseline criterion — cost of doing nothing, evaluated at k+2
    gen_out      = abs(Q_gen_cur) > Q_lim; % generators currently outside the soft limit
    current_crit = sum((V_load_cur(target_in_load) - Yc_V(:, k+2)).^2) ...
                 + alpha_tap * sum(Q_gen_cur(gen_out).^2);

    %% Enumerate candidate actions: +step or -step on each eligible transformer
    % candidate(row, :) = [local transformer index, direction (-1 or +1), cost]
    candidate = zeros(2*n_oltc, 3);

    for j = 1:n_oltc
        if ~eligible(j) % transformer blocked: mark both actions as infeasible
            candidate(2*j-1, 3) = inf;
            candidate(2*j,   3) = inf;
            continue;
        end
        for d = 1:2
            direction = 2*d - 3;    % d=1 -> -1 (tap down), d=2 -> +1 (tap up)
            delta     = direction * step_size;

            % Predict V_load and Q_gen at k+2 if this tap action is applied
            V_pred = V_load_cur + Cv_tap(load_idx, transformer_indices(j)) * delta;
            Q_pred = Q_gen_cur  + Cq_tap(pv_idx,   transformer_indices(j)) * delta;

            row = 2*(j-1) + d;
            candidate(row, 1:2) = [j, direction];

            % Reject: tap would go out of physical bounds
            current_tap = mpc.branch(transformer_indices(j), TAP);
            if current_tap + delta < tap_min || current_tap + delta > tap_max
                candidate(row, 3) = inf;
                continue;
            end

            % Reject: hard voltage or reactive power constraints violated
            if any(V_pred > V_load_max) || any(V_pred < V_load_min) || ...
               any(Q_pred > Q_gen_max)  || any(Q_pred < Q_gen_min)
                candidate(row, 3) = inf;
                continue;
            end

            % Compute criterion for this candidate action
            new_gen_out       = abs(Q_pred) > Q_lim;
            candidate(row, 3) = sum((V_pred(target_in_load) - Yc_V(:, k+2)).^2) ...
                               + alpha_tap * sum(Q_pred(new_gen_out).^2);
        end
    end

    %% Select best action — only if it strictly improves the baseline criterion
    [best_cost, best_row] = min(candidate(:, 3));

    % Default outputs: no tap movement
    % Returning V_reg as target keeps every transformer inside its deadband in run_pf_oltc_step
    target_v  = V_reg;
    oltc_idx  = 0;
    direction = 0;

    if best_cost < current_crit
        best_j   = candidate(best_row, 1); % local index of the transformer that acts (1..n_oltc)
        best_dir = candidate(best_row, 2); % direction: -1 = tap decrease, +1 = tap increase
        oltc_idx  = best_j;
        direction = best_dir;

        % Encode direction as a target voltage offset for run_pf_oltc_step.
        % The function moves the tap by one step when |V_reg - target| > tol:
        %   tap decrease (a smaller -> V_LV rises): set target above V_reg
        %   tap increase (a larger  -> V_LV falls): set target below V_reg
        if best_dir < 0
            target_v(best_j) = V_reg(best_j) + 2*tol;
        else
            target_v(best_j) = V_reg(best_j) - 2*tol;
        end

        % Update x_estimated_k2 to reflect the tap effect — passed to add_cap (priority 3)
        % so it evaluates capacitor actions on a state already corrected for this tap move
        x_estimated_k2(1:n_load) = x_estimated_k2(1:n_load) ...
            + Cv_tap(load_idx, transformer_indices(best_j)) * best_dir * step_size;
        x_estimated_k2(n_load+1:n_load+n_gen) = x_estimated_k2(n_load+1:n_load+n_gen) ...
            + Cq_tap(pv_idx, transformer_indices(best_j)) * best_dir * step_size;

        oltc_availability(best_j) = t_oltc; % reset cooldown for this transformer
    end
end