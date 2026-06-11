function [Cv,Cq] = Cal_CvCq(mpc,JAC)
% Cv(i,j): sensitivity of voltage at load bus i to a unit increase of voltage setpoint at generator j
% Cq(i,j): sensitivity of reactive power at generator i to a unit increase of voltage setpoint at generator j
%
% Derivation from the power flow Jacobian (makeJac, full form):
%   The full Jacobian J maps [dTheta; dV] -> [dP; dQ]
%   Submatrices used (Q rows only, voltage columns only, indices shifted by n_bus):
%     dQl/dVl : effect of load voltages on load reactive injections
%     dQl/dVg : effect of generator voltages on load reactive injections
%     dQg/dVg : effect of generator voltages on generator reactive injections
%     dQg/dVl : effect of load voltages on generator reactive injections
%
%   At PQ buses, imposing dQl=0 (fixed load) gives:
%     Cv = -(dQl/dVl)^{-1} * (dQl/dVg)
%   Generator reactive power change follows from voltage propagation:
%     Cq = dQg/dVg + (dQg/dVl)*Cv

define_constants;

pv_and_slack_idx=find(ismember(mpc.bus(:,BUS_TYPE),[PV,REF]));
load_idx=find(mpc.bus(:,BUS_TYPE)==PQ);
n_bus=size(mpc.bus,1);

% Extract voltage-magnitude submatrices of the Q block
dQldVg = JAC(n_bus+load_idx,            n_bus+pv_and_slack_idx); % dQl/dVg
dQldVl = JAC(n_bus+load_idx,            n_bus+load_idx);          % dQl/dVl
Cv_temp = -dQldVl \ dQldVg;                                     % dVl/dVg, from dQl=0

dQgdVg = JAC(n_bus+pv_and_slack_idx,    n_bus+pv_and_slack_idx); % dQg/dVg
dQgdVl = JAC(n_bus+pv_and_slack_idx,    n_bus+load_idx);          % dQg/dVl
Cq_temp = dQgdVg + dQgdVl*Cv_temp;                              % total dQg/dVg

%%
Cv = sparse(n_bus, n_bus);
Cv(load_idx, pv_and_slack_idx) = Cv_temp;

Cq = sparse(n_bus, n_bus);
Cq(pv_and_slack_idx, pv_and_slack_idx) = Cq_temp;
end