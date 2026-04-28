function [Cv_cap,Cq_cap] = Cal_Cv_capCq_cap(mpc,x,cap)
%Cv_cap(i,j) reflects the increase in voltage in load i when a
%capacitor is added on load j
%Cq_cap(i,j) reflects the icrease in reactive power in generator i
%when a capacitor is added on the load j

define_constants;
baseMVA=mpc.baseMVA;

load_idx=find(mpc.bus(:,BUS_TYPE)==PQ);
pv_and_slack_idx=find(ismember(mpc.bus(:,BUS_TYPE),[PV,REF]));
target_idx=find(mpc.bus(:,PD)~=0); %Bus indices of targeted loads

n_target=size(target_idx,1);
n_bus=size(mpc.bus,1);
n_load=size(load_idx,1);

target=find(ismember(load_idx,target_idx));

cap=cap/baseMVA;

fullJac=true;
J = makeJac(mpc,fullJac);

dQlDCap=zeros(n_load,n_target);
dQlDCap(target,:)=diag(x(target))^2*cap;
dQldVl=J(n_bus+load_idx,n_bus+load_idx);
dQgdVl=J(n_bus+pv_and_slack_idx,n_bus+load_idx);

Cv_cap_temp=dQldVl\dQlDCap;
Cq_cap_temp=dQgdVl*Cv_cap_temp;

%%
Cv_cap = sparse(n_bus, n_bus);
Cv_cap(load_idx, target_idx) = Cv_cap_temp;

Cq_cap = sparse(n_bus, n_bus);
Cq_cap(pv_and_slack_idx, target_idx) = Cq_cap_temp;
end