function [mpc_out,results] = custom_runpf(mpc_in,u_pv,u_oltc_V,u_cap,parameters,oltc_idx,cap_idx)
%mpc_in: The network at step k
%U_pv: PV control signal at step k (V_pv(k+1)=V_pv(k)+U_pv
%U_oltc_V: OLTC control signal (0.95pu 1pu or 1.05pu)
%U_cap: 63kV capacitors control signal (+1 0 -1)

%results: Network at step k (before tap adjustment
%because of the delay) with powerflow ran
%mpc_out: Network at step k+1 with tap adjusted before powerflow and other
%control signal update

%% Parameters definition
define_constants;

tap_step_size=parameters(1);
tol_tap=parameters(2);
tap_min=parameters(3);
tap_max=parameters(4);
cap=parameters(5);

pv_idx = find(mpc_in.bus(:, BUS_TYPE) == PV);
secondary_oltc_idx = mpc_in.branch(oltc_idx,T_BUS);

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
% Tap is adjusted for the next step because of the delay
mpc_out=mpc_in;
v_actuals = results.bus(secondary_oltc_idx, VM); %Check voltage for tap adjustment
diffs = v_actuals - u_oltc_V;

for k = 1:length(u_oltc_V)
    if abs(diffs(k)) > tol_tap
        if diffs(k) < 0
            % Voltage low -> Decrease Ratio on LV side to boost voltage
            % V2(LV) = V1(HV)/ratio
            mpc_out.branch(oltc_idx(k), TAP) = mpc_out.branch(oltc_idx(k), TAP) - tap_step_size;
        else
            % Voltage high -> Increase Ratio on LV side to lower voltage
            mpc_out.branch(oltc_idx(k), TAP) = mpc_out.branch(oltc_idx(k), TAP) + tap_step_size;
        end
    end
end

% Enforce physical limits [0.9, 1.1]
mpc_out.branch(oltc_idx, TAP) = max(tap_min, min(tap_max, mpc_out.branch(oltc_idx, TAP)));

end



