function [Cv,Cq] = Cal_CvCq(mpc)
%Cv(i,j) reflects the increase in voltage in load i when the voltage in
%generator j is increased
%Cq(i,j) reflects the increase in voltage in generator i when the voltage in
%generator j is increased

%% Compute analytical voltage and reactive power sensitivities.
define_constants;

pv_and_slack_idx=find(ismember(mpc.bus(:,BUS_TYPE),[PV,REF]));
load_idx=find(mpc.bus(:,BUS_TYPE)==PQ);
n_bus=size(mpc.bus,1);

fullJac=true;
J = makeJac(mpc,fullJac);

dQldVg=J(n_bus+load_idx,n_bus+pv_and_slack_idx);
dQldVl=J(n_bus+load_idx,n_bus+load_idx);
Cv_temp = -dQldVl\dQldVg;

dQgdVg=J(n_bus+pv_and_slack_idx,n_bus+pv_and_slack_idx);
dQgdVl=J(n_bus+pv_and_slack_idx,n_bus+load_idx);
Cq_temp = (dQgdVg + dQgdVl*Cv_temp);

%%
Cv = sparse(n_bus, n_bus);
Cv(load_idx, pv_and_slack_idx) = Cv_temp;

Cq = sparse(n_bus, n_bus);
Cq(pv_and_slack_idx, pv_and_slack_idx) = Cq_temp;
end