clear 'all'
mpc = case9xx_Bsh();
define_constants;
mpopt = mpoption('verbose',0,'out.all',0,'pf.enforce_q_lims',1);
baseMVA = mpc.baseMVA;
load_idx=find(mpc.bus(:,BUS_TYPE)==PQ); % Bus indices of loads 
pv_and_slack_idx=find(ismember(mpc.bus(:,BUS_TYPE),[PV,REF])); % Bus indices of generators and slacks
target_idx=find(mpc.bus(:,PD)~=0); %Bus indices of targeted loads
transformer_indices = [16; 17]; % Branch indices of transformes (Buses 10-5 and 11-6)

n_bus=size(mpc.bus,1);
n_load=size(load_idx,1);
n_gen=size(mpc.gen,1);
n_branch=size(mpc.branch,1);
n_transfo=size(transformer_indices,1);
n_target=size(target_idx,1);
cap=5;

gen_number=zeros(n_bus,1); %gen_number(i) = idx in mpc.gen of generator of bus i 
for i =1:size(mpc.gen,1)
    gen_number(mpc.gen(i,1))=i;
end

[Cv,Cq]=Cal_CvCq(mpc);

results0=runpf(mpc,mpopt);
x=results0.bus(load_idx,VM);
[Cv_cap,Cq_cap] = Cal_Cv_capCq_cap(mpc,x,cap);

[Cv_tap,Cq_tap]=Cal_Cv_tapCq_tap(mpc,x);

disp('Percentage of error')

%% Cv
x0=results0.bus(load_idx,VM);
Cv_exp_temp=zeros(n_load,n_gen);
Cv_exp=zeros(n_bus,n_bus);

for i= 1:n_gen
    mpc1=mpc;
    mpc1.gen(i,VG)=mpc1.gen(i,VG)+0.1;
    results1=runpf(mpc1,mpopt);
    x1=results1.bus(load_idx,VM);
    Cv_exp_temp(:,i)=(x1-x0)*100;
end 
for i=1:n_load
    for j=1:n_gen
        Cv_exp(load_idx(i),pv_and_slack_idx(j))=Cv_exp_temp(i,j);
    end
end

Cv_exp=sparse(Cv_exp);
testCv=(Cv_exp-Cv)./(Cv_exp+10e-10)*100;
disp('--------------------------')
disp('Cv')
disp(testCv)

%% Cq
x0=results0.gen(gen_number(pv_and_slack_idx),QG)/baseMVA;
Cq_exp_temp=zeros(n_gen,n_gen);
Cq_exp=zeros(n_bus,n_bus);
for i = 1:n_gen
    mpc1=mpc;
    mpc1.gen(i,VG)=mpc1.gen(i,VG)+0.1;
    results1=runpf(mpc1,mpopt);
    x1=results1.gen(gen_number(pv_and_slack_idx),QG)/baseMVA;
    Cq_exp_temp(:,i)=(x1-x0)*100;
end
for i=1:n_gen
    for j=1:n_gen
        Cq_exp(pv_and_slack_idx(i),pv_and_slack_idx(j))=Cq_exp_temp(i,j);
    end
end

Cq_exp=sparse(Cq_exp);
testCq=(Cq_exp-Cq)./(Cq_exp+10e-10)*100;

disp('--------------------------')
disp('Cq')
disp(testCq)

%% Cv_cap
x0=results0.bus(load_idx,VM);
Cv_cap_exp_temp=zeros(n_load,n_target);
Cv_cap_exp=zeros(n_bus,n_bus);

for i=1:n_target
    mpc1=mpc;
    mpc1.bus(target_idx(i),BS)=mpc1.bus(target_idx(i),BS)+cap;
    results1=runpf(mpc1,mpopt);
    x1=results1.bus(load_idx,VM);
    Cv_cap_exp_temp(:,i)=(x1-x0);
end
for i =1:n_load
    for j=1:n_target
        Cv_cap_exp(load_idx(i),target_idx(j))=Cv_cap_exp_temp(i,j);
    end
end

Cv_cap_exp=sparse(Cv_cap_exp);
testCv_cap=(Cv_cap_exp-Cv_cap)./(Cv_cap_exp+10e-10)*100;

disp('--------------------------')
disp('Cv_cap')
disp(testCv_cap)
%% Cq_cap
x0=results0.gen(gen_number(pv_and_slack_idx),QG)/baseMVA;
Cq_cap_exp_temp=zeros(n_gen,n_target);
Cq_cap_exp=zeros(n_bus,n_bus);

for i=1:n_target
    mpc1=mpc;
    mpc1.bus(target_idx(i),BS)=mpc1.bus(target_idx(i),BS)+cap;
    results1=runpf(mpc1,mpopt);
    x1=results1.gen(gen_number(pv_and_slack_idx),QG)/baseMVA;
    Cq_cap_exp_temp(:,i)=(x1-x0);
end
for i=1:n_gen
    for j=1:n_target
        Cq_cap_exp(pv_and_slack_idx(i),target_idx(j))=Cq_cap_exp_temp(i,j);
    end
end

Cq_cap_exp=sparse(Cq_cap_exp);
testCq_cap=(Cq_cap_exp-Cq_cap)./(Cq_cap_exp+10e-10)*100;

disp('--------------------------')
disp('Cq_cap')
disp(testCq_cap)

%% Cv_tap
x0=results0.bus(load_idx,VM);
step_size = 0.00625; % 0.625% fixed step size
Cv_tap_exp_temp=zeros(n_load,n_transfo);
Cv_tap_exp=zeros(n_bus,n_branch);

for i=1:n_transfo
    mpc1=mpc;
    mpc1.branch(transformer_indices(i),TAP)=mpc1.branch(transformer_indices(i),TAP)+step_size;
    results1=runpf(mpc1,mpopt);
    x1=results1.bus(load_idx,VM);
    Cv_tap_exp_temp(:,i)=(x1-x0)/step_size;
end

for i=1:n_load
    for j=1:n_transfo
        Cv_tap_exp(load_idx(i),transformer_indices(j))=Cv_tap_exp_temp(i,j);
    end
end

Cv_tap_exp=sparse(Cv_tap_exp);
testCv_tap=(Cv_tap_exp-Cv_tap)./(Cv_tap_exp+10e-10)*100;

disp('--------------------------')
disp('Cv_tap')
disp(testCv_tap)

%% Cq_tap
x0=results0.gen(gen_number(pv_and_slack_idx),QG)/baseMVA;
step_size = 0.00625; % 0.625% fixed step size
Cq_tap_exp_temp=zeros(n_gen,n_transfo);
Cq_tap_exp=zeros(n_bus,n_branch);

for i=1:n_transfo
    mpc1=mpc;
    mpc1.branch(transformer_indices(i),TAP)=mpc1.branch(transformer_indices(i),TAP)+step_size;
    results1=runpf(mpc1,mpopt);
    x1=results1.gen(gen_number(pv_and_slack_idx),QG)/baseMVA;
    Cq_tap_exp_temp(:,i)=(x1-x0)/step_size;
end

for i=1:n_gen
    for j=1:n_transfo
        Cq_tap_exp(pv_and_slack_idx(i),transformer_indices(j))=Cq_tap_exp_temp(i,j);
    end
end

Cq_tap_exp=sparse(Cq_tap_exp);
testCq_tap=(Cq_tap_exp-Cq_tap)./(Cq_tap_exp+10e-10)*100;

disp('--------------------------')
disp('Cq_tap')
disp(testCq_tap)