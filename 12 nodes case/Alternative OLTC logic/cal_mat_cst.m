function [GAMMA_cst,f1_cst,f2_cst,b1,b2,Hu,Hx,PSI,C_tilde,A,C,J1] = cal_mat_cst(mpc,N,u_max,alpha)
% Computes MPC matrices that do not depend on the current network state.
% These only need to be computed once at initialization.
%
% State:  x(k) = [V_load (n_load x 1) ; Q_gen (n_gen x 1) ; V_gen (n_gen x 1)]
% Output: y(k) = [V_load_target (n_target x 1) ; Q_gen (n_gen x 1)]
% Input:  u(k) = delta_V_gen (n_gen x 1), voltage setpoint increments
%
% Outputs:
%   GAMMA_cst : base convolution matrix (n_state*N x n_state*N), A^(i-j) blocks;
%               multiply by kron(eye(N),B) in Cal_mat_var to get the true GAMMA
%   f1_cst, f2_cst : partial gradient terms (input-independent parts)
%   b1, b2    : QP constraint vectors: Aineq*U <= b1 + b2*x(k)
%   Hu, Hx    : constraint selection matrices for inputs and states
%   PSI       : block-diagonal weighting matrix (1/alpha on V, 1 on Q, scaled for conditioning)
%   C_tilde   : block-diagonal output matrix over the horizon
%   A         : state transition matrix (identity: linearized model assumption)
%   C         : output matrix extracting [V_load_target ; Q_gen] from x
%   J1        : C_tilde * PI, used to compute the open-loop predicted output

%% Preparation
define_constants;
baseMVA = mpc.baseMVA;

pv_idx=find(mpc.bus(:,BUS_TYPE)==PV); % Bus indices of generators
load_idx=find(mpc.bus(:,BUS_TYPE)==PQ); % Bus indices of loads 
target_idx=find(mpc.bus(:,PD)~=0); %Bus indices of targeted loads

n_load=size(load_idx,1);
n_gen=size(pv_idx,1);
n_target=size(target_idx,1);

Q_gen_max=mpc.gen(pv_idx,QMAX)/baseMVA;
Q_gen_min=mpc.gen(pv_idx,QMIN)/baseMVA;

V_gen_max=mpc.bus(pv_idx,VMAX);
V_gen_min=mpc.bus(pv_idx,VMIN);

V_load_max=mpc.bus(load_idx,VMAX);
V_load_min=mpc.bus(load_idx,VMIN);


%% Computation

% A=I: linearized model assumption (sensitivity matrices capture all dynamics)
% x evolves as x(k+1) = A*x(k) + B*u(k), with B computed in Cal_mat_var
A=eye(n_load+n_gen+n_gen);

% C extracts the controlled outputs from the state:
%   rows 1..n_target : select V_load at targeted buses
%   rows n_target+1..end : select Q_gen
C=zeros(n_target+n_gen,n_load+n_gen+n_gen);
[~,target_in_load]=ismember(target_idx,load_idx); %positions within load_idx
for i =1:n_target
    C(i,target_in_load(i))=1;
end
C(n_target+1:end,n_load+1:n_load+n_gen)=eye(n_gen);

% Build stacked prediction matrices over horizon N:
%   X = PI*x(k) + GAMMA*U
% where X = [x(k+1);...;x(k+N)], U = [u(k);...;u(k+N-1)]
% PI(i,:)    = A^(i+1)               (free response from initial state)
% GAMMA(i,j) = A^(i-j)*B for j<=i   (forced response; B factored out via kron in Cal_mat_var)
PI = zeros((n_load+n_gen+n_gen)*N,n_load+n_gen+n_gen);
GAMMA_cst = zeros((n_load+n_gen+n_gen)*N,(n_load+n_gen+n_gen)*N);
for i = 0:N-1
    row = i*(n_load+n_gen+n_gen)+1:(i+1)*(n_load+n_gen+n_gen);
    PI(row,:)=A^(i+1);
    for j=0:i
        col = j*(n_load+n_gen+n_gen)+1:(j+1)*(n_load+n_gen+n_gen);
        GAMMA_cst(row,col)=A^(i-j);
    end
end

% Inequality constraints in the form: Aineq*U <= b1 + b2*x(k)
% - Input constraints:  -u_max <= u(k) <= u_max  at each horizon step  (Hu, bu)
% - State constraints:  V_load, Q_gen, V_gen within their bounds at each predicted step (Hx, bx)
% Hbx_* are selection matrices that extract one variable type from the state vector x
Hu=[eye(N*n_gen);-eye(N*n_gen)];
Hbx_vload=[eye(n_load),zeros(n_load,2*n_gen)];
Hbx_q=[zeros(n_gen,n_load),eye(n_gen),zeros(n_gen,n_gen)];
Hbx_vgen=[zeros(n_gen,n_load+n_gen),eye(n_gen)];
Hx=[kron(eye(N),Hbx_vload);-kron(eye(N),Hbx_vload);kron(eye(N),Hbx_q);-kron(eye(N),Hbx_q);kron(eye(N),Hbx_vgen);-kron(eye(N),Hbx_vgen)];
bu=u_max*ones(2*n_gen*N,1);
bx=[repmat(V_load_max,N,1);-repmat(V_load_min,N,1);repmat(Q_gen_max,N,1);-repmat(Q_gen_min,N,1);repmat(V_gen_max,N,1);-repmat(V_gen_min,N,1)];

b1=[bu;bx];
b2=[zeros(2*N*n_gen,n_load+n_gen+n_gen);Hx*PI];

% Weighting matrices for the MPC cost: J = sum_k [ ||V_target - V_load||^2 + alpha*||Q_gen||^2 ]
% All weights are scaled by 1/alpha to improve Hessian conditioning (eigenvalues ~1 instead of ~alpha).
% This yields an equivalent optimisation problem with the same optimal U*.
Q=diag([ones(1,n_target)/alpha, ones(1,n_gen)]);
C_tilde=kron(eye(N),C);
PSI=kron(eye(N),Q);

% QP cost pre-computation (input-independent parts):
%   min_U  0.5*U'*H*U + (f1*x - f2*Yc)'*U
% f1_cst and f2_cst are completed in Cal_mat_var once GAMMA is known
f1_cst=C_tilde'*PSI*C_tilde*PI;
f2_cst=C_tilde'*PSI';

J1=C_tilde*PI;

end