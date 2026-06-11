function [Cv_cap,Cq_cap] = Cal_Cv_capCq_cap(mpc,JAC,voltage)
%Cv(i,j) reflects the increase in voltage in pq bus i when the shunt susceptance of the
%bus j is increased by 1pu
%Cq(i,j) reflects the increase in reactive power in pv bus i when the shunt susceptance of the
% bus j is increased by 1 pu
%voltage=voltage in every bus

define_constants;

pq_idx=find(mpc.bus(:,BUS_TYPE)==PQ);
pv_and_slack_idx=find(ismember(mpc.bus(:,BUS_TYPE),[PV,REF]));

n_bus=size(mpc.bus,1);

dQpqDCAPpq=diag(voltage(pq_idx))^2;
dQpqdVpq=JAC(n_bus+pq_idx,n_bus+pq_idx);
dQpvdVpq=JAC(n_bus+pv_and_slack_idx,n_bus+pq_idx);

Cv_cap_pq_pq=dQpqdVpq\dQpqDCAPpq;
Cq_cap_pq=dQpvdVpq*Cv_cap_pq_pq;
Cq_cap_pv=-diag(voltage(pv_and_slack_idx))^2;

%%
Cv_cap = sparse(n_bus, n_bus);
Cv_cap(pq_idx, pq_idx) = Cv_cap_pq_pq;

Cq_cap = sparse(n_bus, n_bus);
Cq_cap(pv_and_slack_idx, pq_idx) = Cq_cap_pq;
Cq_cap(pv_and_slack_idx, pv_and_slack_idx) = Cq_cap_pv;

%Cv_cap(abs(Cv_cap)<0.01)=0;
%Cq_cap(abs(Cq_cap)<0.01)=0;
end