function [mpc_out, results] = custom_runpf(mpc_in, u_pv, u_oltc_V, u_cap, parameters, oltc_idx, cap_idx)
% mpc_in   : Network at step k
% u_pv     : PV control signal at step k (V_pv(k+1) = V_pv(k) + u_pv)
% u_oltc_V : OLTC voltage setpoint (0.95, 1.00 or 1.05 pu), size n_oltc
% u_cap    : 63 kV capacitor control signal (+1, 0, -1)
%
% mpc_out  : Network at step k+1, with tap adjusted before next power flow
% results  : Power flow results at step k (before tap adjustment,
%            because of the one-step delay)

% oltc_counter and prev_u_oltc_V are kept as persistent variables so that
% the OLTC state is fully encapsulated inside this function.
% They are initialised automatically on the first call.
persistent oltc_counter prev_u_oltc_V

n_oltc = length(oltc_idx);

if isempty(oltc_counter)
    oltc_counter  = zeros(n_oltc, 1);
    prev_u_oltc_V = u_oltc_V;
end

%% Parameters definition
define_constants;

tap_step_size = parameters(1);
tol_tap       = parameters(2);
tap_min       = parameters(3);
tap_max       = parameters(4);
cap           = parameters(5);
first_delay_oltc = parameters(6);

pv_idx             = find(mpc_in.bus(:, BUS_TYPE) == PV);
secondary_oltc_idx = mpc_in.branch(oltc_idx, T_BUS);

n_bus=size(mpc_in.bus,1);
n_gen=size(mpc_in.gen,1);

% gen_number(i) = row index in mpc.gen of the generator at bus i (0 if none)
gen_number = zeros(n_bus, 1);
for i = 1:n_gen
    gen_number(mpc_in.gen(i,1)) = i;
end

V_pv_min=mpc_in.bus(pv_idx,VMIN);
V_pv_max=mpc_in.bus(pv_idx,VMAX);
%% Voltage setpoints update
mpc_in.gen(gen_number(pv_idx), VG) = min(max(mpc_in.gen(gen_number(pv_idx), VG) + u_pv, V_pv_min), V_pv_max);

%% Capacitor update
mpc_in.bus(cap_idx,BS)=mpc_in.bus(cap_idx,BS)+cap*u_cap;

%%    % --- Run Power Flow  ---
mpopt = mpoption('verbose', 0,'out.all', 0,'pf.enforce_q_lims', 1);
results = runpf(mpc_in, mpopt);

%% --- OLTC control ---
% The tap computed here takes effect at the next power flow call, due to the
% inherent one-step delay of the state-space model.
%
% Each OLTC runs an independent persistence counter (oltc_counter) that tracks
% how many consecutive steps the secondary voltage has remained out of tolerance
% on the same side of the setpoint. The tap action rules are:
%
%   First tap action
%     Fires as soon as |counter| reaches first_delay_oltc. This means the
%     voltage must have stayed out of tolerance on the same side for
%     first_delay_oltc consecutive steps with an unchanged setpoint.
%
%   Subsequent tap actions
%     If the secondary voltage is still out of tolerance on the same side after
%     a tap action, |counter| keeps incrementing and a new tap step is applied
%     every step (no additional delay is needed).
%
%   Counter reset to 0 (full delay required again before next action)
%     - Voltage returns within tolerance.
%     - Voltage crosses to the opposite side of the setpoint.
%     - The voltage setpoint changes.
%
%   Counter reset to +/-1 (side change while already out of tolerance)
%     - Voltage was out of tolerance and crosses to the opposite side.
%       The counter restarts from 1 (or -1) so the full delay is required
%       before the first tap action in the new direction.

mpc_out = mpc_in;
v_actuals = results.bus(secondary_oltc_idx, VM);
diffs     = v_actuals - u_oltc_V;

for k = 1:length(u_oltc_V)

    % --- Setpoint change: reset counter so the full delay is required again ---
    if u_oltc_V(k) ~= prev_u_oltc_V(k)
        oltc_counter(k) = 0;
    end

    % --- Update persistence counter ---
    if abs(diffs(k)) > tol_tap
        new_sign = sign(diffs(k));   % +1 if voltage above setpoint, -1 if below
        if oltc_counter(k) == 0 || sign(oltc_counter(k)) == new_sign
            % Voltage out of tolerance on the same side: increment counter
            oltc_counter(k) = oltc_counter(k) + new_sign;
        else
            % Voltage switched side: restart count from +/-1
            oltc_counter(k) = new_sign;
        end
    else
        % Voltage back within tolerance: reset
        oltc_counter(k) = 0;
    end

    % --- Tap action ---
    % |counter| == first_delay_oltc : first action, after the required delay has elapsed
    % |counter|  > first_delay_oltc : voltage still out of tolerance on the same side,
    %                                 one tap step applied every step (no additional delay)
    if abs(oltc_counter(k)) >= first_delay_oltc
        if diffs(k) < 0
            % TAP in matpower = HV/LV
            % Voltage low -> decrease tap ratio to boost LV voltage (LV = HV/ratio)
            mpc_out.branch(oltc_idx(k), TAP) = mpc_out.branch(oltc_idx(k), TAP) - tap_step_size;
        else
            % Voltage high -> increase tap ratio to lower LV voltage
            mpc_out.branch(oltc_idx(k), TAP) = mpc_out.branch(oltc_idx(k), TAP) + tap_step_size;
        end
    end

end

% Enforce physical tap limits
mpc_out.branch(oltc_idx, TAP) = max(tap_min, min(tap_max, mpc_out.branch(oltc_idx, TAP)));

% Store current setpoint for next step
prev_u_oltc_V = u_oltc_V;

end