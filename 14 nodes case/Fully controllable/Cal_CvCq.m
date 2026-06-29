function [Cv,Cq] = Cal_CvCq(mpc,JAC)
%Cv(i,j) reflects the increase in voltage in pq bus i when the voltage in
%pv bus j is increased by 1pu
%Cq(i,j) reflects the increase in reactive power in pv bus i when the voltage in
%pv bus j is increased by 1 pu

%% Compute analytical voltage and reactive power sensitivities.
define_constants;

pv_and_slack_idx=find(ismember(mpc.bus(:,BUS_TYPE),[PV,REF]));
pq_idx=find(mpc.bus(:,BUS_TYPE)==PQ);
n_bus=size(mpc.bus,1);

dQpqdVpv=JAC(n_bus+pq_idx,n_bus+pv_and_slack_idx);
dQpqdVpq=JAC(n_bus+pq_idx,n_bus+pq_idx);
Cv_temp = -dQpqdVpq\dQpqdVpv;

dQpvdVpv=JAC(n_bus+pv_and_slack_idx,n_bus+pv_and_slack_idx);
dQpvdVpq=JAC(n_bus+pv_and_slack_idx,n_bus+pq_idx);
Cq_temp = (dQpvdVpv + dQpvdVpq*Cv_temp);

%%
Cv = sparse(n_bus, n_bus);
Cv(pq_idx, pv_and_slack_idx) = Cv_temp;
%Cv(abs(Cv)<0.01)=0;

Cq = sparse(n_bus, n_bus);
Cq(pv_and_slack_idx, pv_and_slack_idx) = Cq_temp;
%Cq(abs(Cq)<0.01)=0;

end