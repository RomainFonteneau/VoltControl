function u_cap_k = add_cap(mpc,x_estimated_k1,Cv_cap,Cq_cap,cap_availability,cap_idx,cap,Q_lim_rate,cap_k)
%This function add or remove the capacitor which minimized the most the
%criterion. The criterion is here defined by the reactive power in
%generator out of the soft limits. This function acts to preserve reactive
%margins.
%Capacitors are present on the 20kV network and the 63kV one.
%Reactances are only present on the 63kV network.

%cap_k: number of capacitors activated at step k

%% Preparation
define_constants;

baseMVA=mpc.baseMVA;

pv_idx    = find(mpc.bus(:,BUS_TYPE)==PV);  
pq_idx  = find(mpc.bus(:,BUS_TYPE)==PQ);  

n_gen = size(mpc.gen,1);
n_bus = size(mpc.bus,1);
n_cap = size(cap_idx,1);
n_pq  = size(pq_idx,1);
n_pv  = size(pv_idx,1);

gen_number = zeros(n_bus, 1);
for i = 1:n_gen
    gen_number(mpc.gen(i,1)) = i;
end

% Hard constraints 
V_pq_max=mpc.bus(pq_idx,VMAX);
V_pq_min=mpc.bus(pq_idx,VMIN);
Q_pv_max =mpc.gen(gen_number(pv_idx),QMAX)/baseMVA;
Q_pv_min =mpc.gen(gen_number(pv_idx),QMIN)/baseMVA;

% Soft constraints
Q_lim_up=Q_lim_rate*Q_pv_max;
Q_lim_down=Q_lim_rate*Q_pv_min;

%% Logic
V_pq_estimated_k1=x_estimated_k1(1:n_pq,:);
Q_pv_estimated_k1=x_estimated_k1(n_pq+1:n_pq+n_pv,:);

% Baseline criterion: sum of Q^2 for generators currently outside [-Q_lim_up,
% Q_lim_down] at k+1
pv_out     = Q_pv_estimated_k1 > Q_lim_up | Q_pv_estimated_k1 < Q_lim_down;
best_crit= sum(Q_pv_estimated_k1(pv_out).^2);
action=0;

for i=1:n_cap
    if cap_availability(cap_idx(i))==0
        %Try to add a capacitor
        new_Q_pv_estimated_k1 = Q_pv_estimated_k1 + cap/baseMVA*Cq_cap(pv_idx, cap_idx(i));
        new_V_pq_estimated_k1 = V_pq_estimated_k1 + cap/baseMVA*Cv_cap(pq_idx, cap_idx(i));
        new_pv_out  = new_Q_pv_estimated_k1 > Q_lim_up | new_Q_pv_estimated_k1 < Q_lim_down;
        new_crit= sum(new_Q_pv_estimated_k1(new_pv_out).^2);
        if all([new_Q_pv_estimated_k1<Q_pv_max;new_Q_pv_estimated_k1>Q_pv_min;new_V_pq_estimated_k1<V_pq_max;new_V_pq_estimated_k1>V_pq_min]) && best_crit>new_crit
           best_crit=new_crit;
           action=1;
           bus_action=i;
        end
        %Try to remove a capacitor
        if mpc.bus(cap_idx(i),BASE_KV)==63 || cap_k(cap_idx(i))>0   
            new_Q_pv_estimated_k1 = Q_pv_estimated_k1 - cap/baseMVA*Cq_cap(pv_idx, cap_idx(i));
            new_V_pq_estimated_k1 = V_pq_estimated_k1 - cap/baseMVA*Cv_cap(pq_idx, cap_idx(i));
            new_pv_out  = new_Q_pv_estimated_k1 > Q_lim_up | new_Q_pv_estimated_k1 < Q_lim_down;
            new_crit= sum(new_Q_pv_estimated_k1(new_pv_out).^2);
            if all([new_Q_pv_estimated_k1<Q_pv_max;new_Q_pv_estimated_k1>Q_pv_min;new_V_pq_estimated_k1<V_pq_max;new_V_pq_estimated_k1>V_pq_min]) && new_crit<best_crit
               best_crit=new_crit;
               action=-1;
               bus_action=i;            
            end
        end
    end
end

u_cap_k=zeros(n_cap,1);
if action ==1
    u_cap_k(bus_action)=1;
elseif action==-1
    u_cap_k(bus_action)=-1;
end

end