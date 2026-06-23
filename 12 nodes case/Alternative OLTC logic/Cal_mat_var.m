function [Aineq,H,f1,f2,B] = Cal_mat_var(mpc,N,GAMMA_cst,Hu,Hx,PSI,C_tilde,f1_cst,f2_cst,JAC)
% Computes MPC matrices that depend on the current network state (sensitivity matrices).
% Must be called at least once before each QP solve; can be called less frequently
% if the network operating point changes slowly (controlled by t_sen in main.m).
%
% B: linearized input-to-state matrix, maps delta_V_gen -> delta_x
%   rows 1..n_load             : voltage sensitivity at load buses      (from Cv)
%   rows n_load+1..n_load+n_gen: reactive power sensitivity at generators (from Cq)
%   rows n_load+n_gen+1..end   : identity (V_gen tracks its setpoint directly)
define_constants;

load_idx=find(mpc.bus(:,BUS_TYPE)==PQ); % Bus indices of loads
pv_idx=find(mpc.bus(:,BUS_TYPE)==PV);
n_gen=size(pv_idx,1);

[Cv, Cq] = Cal_CvCq(mpc,JAC);

B=[Cv(load_idx,pv_idx);Cq(pv_idx,pv_idx);eye(n_gen)];

% GAMMA = GAMMA_cst * kron(eye(N),B): full forced-response matrix over the horizon
GAMMA=GAMMA_cst*kron(eye(N),B);

% QP inequality constraints: Aineq*U <= b1 + b2*x(k)
Aineq=[Hu;Hx*GAMMA];

% QP cost matrices: min_U  0.5*U'*H*U + (f1*x - f2*Yc)'*U
H=GAMMA'*C_tilde'*PSI*C_tilde*GAMMA;
H=(H+H')/2; % Enforce exact symmetry to avoid numerical issues in quadprog
f1=GAMMA'*f1_cst;
f2=GAMMA'*f2_cst;

end