function [Cv_tap,Cq_tap] = Cal_Cv_tapCq_tap(mpc,JAC,voltage,TAP_estimated)
%Cv_tap(i,j) reflects the increase in voltage in pq bus i when the tap  of oltc branch j is
%increased by 1pu 
%Cq_tap(i,j) reflects the increased in reactive power in pv bus i when the tap of oltc branch j
%is increased by 1pu
%voltage: Voltage of every bus

define_constants;

pq_idx=find(mpc.bus(:,BUS_TYPE)==PQ);
pv_and_slack_idx=find(ismember(mpc.bus(:,BUS_TYPE),[PV,REF]));

n_branch=size(mpc.branch,1);
n_bus=size(mpc.bus,1);

dQpvdVpq=JAC(n_bus+pv_and_slack_idx,n_bus+pq_idx);
dQpqdVpq=JAC(n_bus+pq_idx,n_bus+pq_idx);

dQdTAP=zeros(n_bus,n_branch);
for k=1:n_branch
    Fbus=mpc.branch(k,F_BUS);
    Tbus=mpc.branch(k,T_BUS);
    X=mpc.branch(k,BR_X);
    TAP_k=TAP_estimated(k);
    dQdTAP(Fbus,k)=2*voltage(Fbus)^2/(X*TAP_k^3)-voltage(Fbus)*voltage(Tbus)/(X*TAP_k^2);
    dQdTAP(Tbus,k)=-voltage(Tbus)*voltage(Fbus)/(X*TAP_k^2);
end
dQpqdTAP = dQdTAP(pq_idx, :);
Cv_tap_temp=sparse(dQpqdVpq\dQpqdTAP);

Cq_tap_temp = dQpvdVpq * Cv_tap_temp - dQdTAP(pv_and_slack_idx, :);

%%
Cv_tap=zeros(n_bus,n_branch);
Cv_tap(pq_idx,:)=Cv_tap_temp;

Cq_tap=zeros(n_bus,n_branch);
Cq_tap(pv_and_slack_idx,:)=Cq_tap_temp;

Cv_tap=sparse(Cv_tap);
Cq_tap=sparse(Cq_tap);

%Cv_tap(abs(Cv_tap)<0.01)=0;
%Cq_tap(abs(Cq_tap)<0.01)=0;


end