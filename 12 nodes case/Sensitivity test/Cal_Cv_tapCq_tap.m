function [Cv_tap,Cq_tap] = Cal_Cv_tapCq_tap(mpc, x, transformer_indices,JAC)
% Cv_tap(i,j): sensitivity of voltage at load bus i to a unit increase of tap on transformer branch j
% Cq_tap(i,j): sensitivity of reactive power at generator i to a unit increase of tap on transformer branch j
%
% Derivation: for a lossless transformer with reactance X and tap ratio a (V_primary/V_secondary),
% the reactive injection at the secondary bus is approximately Q ≈ V_secondary^2 / (a^2 * X).
% Differentiating with respect to a at nominal tap (a ≈ 1):
%   dQ_secondary / d_tap ≈ -V_secondary^2 / X
% This direct injection change is then propagated through the network via the Jacobian,
% following the same approach as Cal_CvCq and Cal_Cv_capCq_cap.
% The primary (HV) bus sensitivity is zeroed: the 400kV network is stiff and its
% voltage is effectively set by the slack, so the tap action is only relevant on the LV side.

define_constants;

load_idx         = find(mpc.bus(:,BUS_TYPE)==PQ);
pv_and_slack_idx = find(ismember(mpc.bus(:,BUS_TYPE),[PV,REF]));
secondary_bus_idx = mpc.branch(transformer_indices, T_BUS);
primary_bus_idx   = mpc.branch(transformer_indices, F_BUS);

n_branch  = size(mpc.branch, 1);
n_bus     = size(mpc.bus, 1);
n_load    = size(load_idx, 1);
n_transfo = size(transformer_indices, 1);

secondary_bus = find(ismember(load_idx, secondary_bus_idx)); % positions within load_idx
primary_bus   = find(ismember(load_idx, primary_bus_idx));

X = mpc.branch(transformer_indices, BR_X); % Transformer reactances (p.u.)

% Direct reactive injection change at each secondary bus: dQ/d_tap = -V^2 / X
dQldTap = zeros(n_load, n_transfo);
dQldTap(secondary_bus, :) = -diag((x(secondary_bus).^2) ./ X);

% Voltage propagation through the network (same derivation as Cal_CvCq)
dQldVl = JAC(n_bus+load_idx,         n_bus+load_idx);
dQgdVl = JAC(n_bus+pv_and_slack_idx, n_bus+load_idx);

Cv_tap_temp = sparse(dQldVl \ dQldTap);
Cv_tap_temp(primary_bus, :) = 0; % HV side: stiff 400kV network, tap has no effective voltage impact
Cq_tap_temp = sparse(dQgdVl * Cv_tap_temp);

%%
Cv_tap = sparse(n_bus, n_branch);
Cv_tap(load_idx, transformer_indices) = Cv_tap_temp;

Cq_tap = sparse(n_bus, n_branch);
Cq_tap(pv_and_slack_idx, transformer_indices) = Cq_tap_temp;

end
