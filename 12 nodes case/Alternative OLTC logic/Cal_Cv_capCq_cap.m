function [Cv_cap,Cq_cap] = Cal_Cv_capCq_cap(mpc,x,cap,JAC)
% Cv_cap(i,j): sensitivity of voltage at load bus i to adding one shunt capacitor on load bus j
% Cq_cap(i,j): sensitivity of reactive power at generator i to adding one shunt capacitor on load bus j

define_constants;
baseMVA=mpc.baseMVA;

load_idx=find(mpc.bus(:,BUS_TYPE)==PQ);
pv_and_slack_idx=find(ismember(mpc.bus(:,BUS_TYPE),[PV,REF]));
target_idx=find(mpc.bus(:,PD)~=0); %Bus indices of targeted loads

n_target=size(target_idx,1);
n_bus=size(mpc.bus,1);
n_load=size(load_idx,1);

target=find(ismember(load_idx,target_idx)); % indices of target buses within load_idx

cap=cap/baseMVA;

% Reactive injection from shunt susceptance B at voltage V: Q = B*V^2
% Adding cap (delta_B = cap) on target bus j changes reactive injection by: delta_Q_j = cap * V_j^2
% dQl/dCap is diagonal on target buses (each cap only directly affects its own bus)
dQlDCap = zeros(n_load, n_target);
dQlDCap(target, :) = diag(x(target))^2 * cap;

% Voltage propagation to all load buses (same derivation as Cal_CvCq)
dQldVl = JAC(n_bus+load_idx,            n_bus+load_idx);
dQgdVl = JAC(n_bus+pv_and_slack_idx,    n_bus+load_idx);

Cv_cap_temp = dQldVl \ dQlDCap;      % dVl/dCap
Cq_cap_temp = dQgdVl * Cv_cap_temp;  % dQg/dCap

%%
Cv_cap = sparse(n_bus, n_bus);
Cv_cap(load_idx, target_idx) = Cv_cap_temp;

Cq_cap = sparse(n_bus, n_bus);
Cq_cap(pv_and_slack_idx, target_idx) = Cq_cap_temp;
end