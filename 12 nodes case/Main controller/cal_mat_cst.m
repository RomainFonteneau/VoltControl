function [C, PSI, C_tilde, Hu, Hx, b1, H_gamma] = cal_mat_cst(mpc, N, u_max, alpha, beta, n_oltc, step_size)
% Computes MPC matrices that do not depend on the current network state.
% These only need to be computed once at initialisation.
%
% --- Augmented state formulation ---
% The OLTC tap commands are now MPC inputs. To handle the two-step delay
% (tap decided at k takes effect at k+2), the previous OLTC command is
% stored as part of the state:
%
%   Augmented state : x_tilde(k) = [V_load (n_load x 1)     ]
%                                   [Q_gen  (n_gen  x 1)     ]
%                                   [V_gen  (n_gen  x 1)     ]
%                                   [u_OLTC(k-1) (n_oltc x 1)]
%
%   Extended input  : u(k) = [Delta_V_gen (n_gen  x 1)]
%                             [u_OLTC     (n_oltc x 1)]
%
%   Dynamics: x_tilde(k+1) = A_tilde * x_tilde(k) + B_tilde * u(k)
%
%     A_tilde = [ I       B_OLTC ]    B_tilde = [ B_gen   0      ]
%               [ 0       0      ]              [ 0       I_oltc ]
%
%   B_OLTC and B_gen are operating-point-dependent and are filled in
%   Cal_mat_var. Here only the structural (constant) matrices are built.
%
%   Key property used in Cal_mat_var: A_tilde^n = A_tilde for all n >= 1,
%   which makes PI and GAMMA trivial to compute from A_tilde and B_tilde.
%
% --- Cost function ---
%   J = sum_{i=1}^{N} [ ||V_load_target - V_ref||^2
%                     + alpha * ||Q_gen||^2
%                     + beta  * ||u_OLTC||^2 ]
%
%   V term has coefficient 1 (dominant objective).
%   alpha penalises reactive power deviation from zero.
%   beta  penalises tap magnitude -> tap used only when generators insufficient.
%
% --- What moved to Cal_mat_var ---
%   GAMMA_cst, f1_cst, f2_cst, b2, A, J1 are no longer returned here because
%   A_tilde depends on B_OLTC (operating point). Cal_mat_var now builds PI,
%   GAMMA, H, f1, f2, b2, J1 entirely.
%
% Inputs:
%   mpc          : MATPOWER case (used for dimensions and bounds only)
%   N            : MPC prediction horizon (steps)
%   u_max        : max generator voltage increment per step (p.u.)
%   alpha        : weight on generator reactive power (V term has coefficient 1)
%   beta         : weight on OLTC command magnitude (>0 makes tap costly)
%   n_oltc       : number of OLTC transformers
%   step_size    : tap step size (p.u.); defines the u_OLTC input bound
%
% Outputs:
%   C            : output matrix  (n_output x n_state_aug)
%                  extracts [V_load_target; Q_gen] from x_tilde
%   PSI          : block-diagonal cost weighting matrix over the horizon
%   C_tilde      : kron(eye(N), C)
%   Hu           : input bound constraint matrix  [I; -I]
%                  (2*N*n_input x N*n_input)
%   Hx           : state constraint selection matrix
%                  rows select V_load, Q_gen, V_gen from predicted states
%   b1           : constant RHS of QP constraints [bu; bx]
%                  (variable part b2 = [0; Hx*PI] is built in Cal_mat_var)
%   H_gamma      : constant H contribution from beta term
%                  H_gamma = beta * S_oltc' * S_oltc

%% Preparation
define_constants;
baseMVA = mpc.baseMVA;

pv_idx     = find(mpc.bus(:, BUS_TYPE) == PV);  % generator buses
load_idx   = find(mpc.bus(:, BUS_TYPE) == PQ);  % load buses
target_idx = find(mpc.bus(:, PD) ~= 0);          % buses with active load (tracked)

n_load   = length(load_idx);
n_gen    = length(pv_idx);
n_target = length(target_idx);

% Dimensions
n_state     = n_load + n_gen + n_gen;  % original state size
n_state_aug = n_state + n_oltc;         % augmented state (adds u_OLTC(k-1))
n_input     = n_gen + n_oltc;           % extended input  (adds u_OLTC)
n_output    = n_target + n_gen;         % outputs: [V_load_target; Q_gen]

% Bounds (from network data)
Q_gen_max  = mpc.gen(pv_idx, QMAX) / baseMVA;
Q_gen_min  = mpc.gen(pv_idx, QMIN) / baseMVA;
V_gen_max  = mpc.bus(pv_idx, VMAX);
V_gen_min  = mpc.bus(pv_idx, VMIN);
V_load_max = mpc.bus(load_idx, VMAX);
V_load_min = mpc.bus(load_idx, VMIN);

%% Output matrix C  (n_output x n_state_aug)
% Extracts [V_load_target (n_target x 1); Q_gen (n_gen x 1)] from x_tilde.
% The last n_oltc columns (u_OLTC(k-1) part of x_tilde) are not observed -> zero.
C = zeros(n_output, n_state_aug);
[~, target_in_load] = ismember(target_idx, load_idx);    % positions within load_idx
for i = 1:n_target
    C(i, target_in_load(i)) = 1;                         % V_load at target bus i
end
C(n_target+1:end, n_load+1:n_load+n_gen) = eye(n_gen);  % Q_gen

%% Weighting matrices
% Q_w: per-step cost matrix on [V_load_target error; Q_gen]
% PSI = kron(eye(N), Q_w): cost over the full horizon
%
% The cost J = ||V_error||^2 + alpha*||Q||^2 + beta*||u_OLTC||^2 is normalised
% by 1/alpha so the Hessian magnitude stays O(1/alpha) regardless of the
% absolute value of alpha. The optimal U* is unchanged by this scaling.
%   V weight   : 1/alpha
%   Q weight   : 1          (= alpha * 1/alpha)
%   OLTC weight: beta/alpha  (applied in H_gamma below)
Q_w     = diag([1/alpha * ones(1, n_target), ones(1, n_gen)]);
PSI     = kron(eye(N), Q_w);
C_tilde = kron(eye(N), C);  % stacks C for all N prediction steps

%% Input bound constraints Hu  (2*N*n_input x N*n_input)
% |Delta_V_gen(k+i)| <= u_max     for all i = 0..N-1
% |u_OLTC(k+i)|      <= step_size  for all i = 0..N-1
% Hard bounds on tap physical limits [tap_min, tap_max] are handled
% dynamically in main.m (unilateral constraints added when at a limit).
Hu = [eye(N*n_input); -eye(N*n_input)];
bu = repmat([u_max * ones(n_gen, 1); step_size * ones(n_oltc, 1)], 2*N, 1);

%% State constraint matrix Hx
% Selects V_load, Q_gen, V_gen from the predicted augmented state.
% The u_OLTC_prev block (last n_oltc rows of x_tilde) is not state-constrained here.
Hbx_vload = [eye(n_load),                    zeros(n_load, 2*n_gen),    zeros(n_load, n_oltc)];
Hbx_q     = [zeros(n_gen, n_load),           eye(n_gen),                zeros(n_gen, n_gen + n_oltc)];
Hbx_vgen  = [zeros(n_gen, n_load + n_gen),   eye(n_gen),                zeros(n_gen, n_oltc)];

Hx = [kron(eye(N),  Hbx_vload); ...   % V_load <= V_load_max
      kron(eye(N), -Hbx_vload); ...   % V_load >= V_load_min
      kron(eye(N),  Hbx_q);    ...   % Q_gen  <= Q_gen_max
      kron(eye(N), -Hbx_q);    ...   % Q_gen  >= Q_gen_min
      kron(eye(N),  Hbx_vgen); ...   % V_gen  <= V_gen_max
      kron(eye(N), -Hbx_vgen)];      % V_gen  >= V_gen_min

bx = [repmat( V_load_max, N, 1); repmat(-V_load_min, N, 1); ...
      repmat( Q_gen_max,  N, 1); repmat(-Q_gen_min,  N, 1); ...
      repmat( V_gen_max,  N, 1); repmat(-V_gen_min,  N, 1)];

% b1: constant part of the RHS.  Full RHS at runtime: b1 - b2*x_tilde
% b2 = [zeros(2*N*n_input, n_state_aug); Hx*PI] is built in Cal_mat_var
% because PI = repmat(A_tilde, N, 1) depends on B_OLTC (operating point).
b1 = [bu; bx];

%% OLTC penalty matrix
% S_oltc  (N*n_oltc x N*n_input)
%   Selects u_OLTC from each block of the stacked input vector U.
%   U = [u(k); u(k+1); ...; u(k+N-1)],  each u(k+i) = [Delta_V_gen; u_OLTC]
S_oltc_block = [zeros(n_oltc, n_gen), eye(n_oltc)];   % (n_oltc x n_input)
S_oltc       = kron(eye(N), S_oltc_block);              % (N*n_oltc x N*n_input)

% H_gamma: constant H contribution from the beta term (quadprog minimises 0.5*U'*H*U + f'*U)
%   0.5 * (beta/alpha) * ||S_oltc * U||^2  ->  H += (beta/alpha) * S_oltc' * S_oltc
H_gamma = (beta/alpha) * (S_oltc' * S_oltc);

end
