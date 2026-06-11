function [J_Q,J_V,J_VQ] = Cal_criterion(mpc,y,Yc_V,Yc_Q,k,alpha)
% Computes performance criteria at time step k (for logging and display only, not used by MPC).
%
% Inputs:
%   y     : current output vector [V_load_target ; Q_gen]
%   Yc_V  : voltage targets over time  (n_target x T)
%   Yc_Q  : reactive power targets (zeros) over time (n_gen x T)
%   k     : current time step (1-based)
%   alpha : weight on Q_gen (matches MPC cost: J = ||V_error||^2 + alpha*||Q_gen||^2)
%
% Outputs:
%   J_Q  : sum of squared reactive power deviations
%   J_V  : sum of squared voltage tracking errors
%   J_VQ : weighted sum: J_V + alpha*J_Q
define_constants;

target_idx=find(mpc.bus(:,PD)~=0);
n_target=size(target_idx,1);

J_Q=sum((y(n_target+1:end,1)-Yc_Q(:,k)).^2);
J_V=sum((y(1:n_target,1)-Yc_V(:,k)).^2);
J_VQ=J_V + alpha*J_Q;
end