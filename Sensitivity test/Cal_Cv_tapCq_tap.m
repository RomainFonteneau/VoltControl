function [Cv_tap,Cq_tap] = Cal_Cv_tapCq_tap(mpc,x)
define_constants;

load_idx=find(mpc.bus(:,BUS_TYPE)==PQ);
pv_and_slack_idx=find(ismember(mpc.bus(:,BUS_TYPE),[PV,REF]));
transformer_indices = [16; 17]; % Branch indices of transformes (Buses 10-5 and 11-6)
secondary_bus_idx=mpc.branch(transformer_indices,T_BUS);
primary_bus_idx=mpc.branch(transformer_indices,F_BUS);

n_branch=size(mpc.branch,1);
n_bus=size(mpc.bus,1);
n_load=size(load_idx,1);
n_transfo=size(transformer_indices,1);

secondary_bus=find(ismember(load_idx,secondary_bus_idx));
primary_bus=find(ismember(load_idx,primary_bus_idx));

fullJac=true;
J = makeJac(mpc,fullJac);

X=mpc.branch(transformer_indices,BR_X);%Reactance
dQldTap=zeros(n_load,n_transfo);
dQldTap(secondary_bus,:)=-diag((x(secondary_bus).^2)./X);
dQldVl=J(n_bus+load_idx,n_bus+load_idx);
dQgdVl=J(n_bus+pv_and_slack_idx,n_bus+load_idx);

Cv_tap_temp=sparse(dQldVl\dQldTap);
Cv_tap_temp(primary_bus,:)=zeros(n_transfo,n_transfo);
Cq_tap_temp=sparse(dQgdVl*Cv_tap_temp);

%%
Cv_tap=zeros(n_bus,n_branch);
for i=1:n_load
    for j=1:n_transfo
        Cv_tap(load_idx(i),transformer_indices(j))=Cv_tap_temp(i,j);
    end
end

Cq_tap=zeros(n_bus,n_branch);
for i=1:size(pv_and_slack_idx,1)
    for j=1:n_transfo
        Cq_tap(pv_and_slack_idx(i),transformer_indices(j))=Cq_tap_temp(i,j);
    end
end

Cv_tap=sparse(Cv_tap);
Cq_tap=sparse(Cq_tap);

end