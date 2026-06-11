function [Aineq, b2, H, f1, f2, B_tilde, A_tilde, J1] = Cal_mat_var(mpc, x, N, C_tilde, Hu, Hx, PSI, H_gamma, JAC, transformer_indices)
% Computes MPC matrices that depend on the current network state (sensitivity matrices).
% Must be called at least once before each QP solve.
%
% Extended from the original to include OLTC taps as MPC inputs.
% Builds the augmented state matrices A_tilde and B_tilde, then computes
% PI and GAMMA analytically using the property A_tilde^n = A_tilde for n >= 1.
%
% --- Augmented model ---
%   Augmented state : x_tilde = [V_load; Q_gen; V_gen; u_OLTC(k-1)]
%   Extended input  : u       = [Delta_V_gen; u_OLTC]
%
%   x_tilde(k+1) = A_tilde * x_tilde(k) + B_tilde * u(k)
%
%   A_tilde = [ I        B_OLTC_x ]    B_tilde = [ B_gen   0      ]
%             [ 0        0        ]              [ 0       I_oltc ]
%
%   B_gen    = [Cv(load,gen); Cq(gen,gen); I_gen]  (from Cal_CvCq)
%   B_OLTC_x = [Cv_tap(load,transfo); Cq_tap(gen,transfo); 0]  (from Cal_Cv_tapCq_tap)
%
%   The tap delay (two steps) is captured by the augmented state:
%   u_OLTC(k) is stored in x_tilde(k+1) via B_tilde, then propagated into x
%   at step k+2 via the B_OLTC_x block of A_tilde.
%
% --- Key property ---
%   A_tilde^n = A_tilde  for all n >= 1
%   Proof: A_tilde^2 = [I, B_OLTC_x; 0, 0]^2 = [I, B_OLTC_x; 0, 0] = A_tilde (by induction).
%   This gives closed-form expressions:
%     PI    = repmat(A_tilde, N, 1)                              (N*n_state_aug x n_state_aug)
%     GAMMA = kron(eye(N), B_tilde) + kron(tril(ones(N),-1), A_tilde*B_tilde)
%             block(i,j): B_tilde if i==j,  A_tilde*B_tilde if i>j   (N*n_state_aug x N*n_input)
%
% --- QP cost (quadprog minimises 0.5*U'*H*U + f'*U) ---
%   H  = GAMMA'*C_tilde'*PSI*C_tilde*GAMMA + H_gamma
%   f1 = GAMMA'*C_tilde'*PSI*C_tilde*PI                      (so f = f1*x_tilde - f2*Yc)
%   f2 = GAMMA'*C_tilde'*PSI
%
% --- QP constraints (Aineq*U <= b1 - b2*x_tilde) ---
%   b2 = [zeros(2*N*n_input, n_state_aug); Hx*PI]
%   Built here because PI depends on A_tilde (operating-point-dependent).
%
% Inputs:
%   mpc                 : current MATPOWER case 
%   x                   : current state  
%   N                   : MPC prediction horizon
%   C_tilde             : kron(eye(N), C)                    (from cal_mat_cst)
%   Hu                  : input bound constraint matrix       (from cal_mat_cst)
%   Hx                  : state constraint selection matrix   (from cal_mat_cst)
%   PSI                 : block-diagonal weighting matrix     (from cal_mat_cst)
%   b1                  : constant RHS of QP constraints      (from cal_mat_cst)
%   H_gamma             : constant H term from beta (OLTC magnitude) penalty (from cal_mat_cst)
%   JAC                 : full network Jacobian (makeJac), shared with all sensitivity functions
%   transformer_indices : branch indices of OLTC transformers (n_oltc x 1)
%
% Outputs:
%   Aineq   : QP inequality constraint matrix  [Hu; Hx*GAMMA]
%   b2      : state-dependent RHS term [0; Hx*PI]  (used as b1 - b2*x_tilde in quadprog)
%   H       : QP cost Hessian (symmetric, enforced numerically)
%   f1      : gradient matrix for state     (N*n_input x n_state_aug)
%   f2      : gradient matrix for reference  (N*n_input x N*n_output)
%   B_tilde : augmented input-to-state matrix (for state prediction in main.m)
%   A_tilde : augmented state-transition matrix (for state prediction in main.m)
%   J1      : C_tilde * PI  (for MPC cost logging in main.m)

    define_constants;

    load_idx  = find(mpc.bus(:, BUS_TYPE) == PQ);
    pv_idx    = find(mpc.bus(:, BUS_TYPE) == PV);

    n_load  = length(load_idx);
    n_gen   = length(pv_idx);
    n_oltc  = length(transformer_indices);

    n_state     = n_load + n_gen + n_gen;
    n_state_aug = n_state + n_oltc;
    n_input     = n_gen + n_oltc;

    %% Sensitivity matrices at the current operating point

    % Generator voltage sensitivities (same as original)
    [Cv, Cq] = Cal_CvCq(mpc, JAC);
    B_gen = [Cv(load_idx, pv_idx); ...   % dV_load / dDelta_V_gen
             Cq(pv_idx,   pv_idx); ...   % dQ_gen  / dDelta_V_gen
             eye(n_gen)];                % V_gen tracks its setpoint directly

    % Tap sensitivities, linearised at the current operating point.
    [Cv_tap, Cq_tap] = Cal_Cv_tapCq_tap(mpc, x, transformer_indices, JAC);

    % B_OLTC_x: columns of the state-transition matrix corresponding to u_OLTC stored
    % in x_tilde.  Describes the effect on [V_load; Q_gen; V_gen] one step later.
    % V_gen is not affected by the tap (HV side is stiff, voltage set by slack).
    B_OLTC_x = [Cv_tap(load_idx, transformer_indices); ...   % dV_load / d_tap
                Cq_tap(pv_idx,   transformer_indices); ...   % dQ_gen  / d_tap
                zeros(n_gen, n_oltc)];                        % V_gen unaffected

    %% Augmented matrices A_tilde and B_tilde

    % A_tilde (n_state_aug x n_state_aug):
    %   Top-left  I       : original state propagation (A = I, linearised model)
    %   Top-right B_OLTC_x: u_OLTC stored at step k-1 acts on x at step k
    %                        (implements the one-step propagation of the tap delay)
    %   Bottom    0       : u_OLTC storage resets; updated by B_tilde with new u_OLTC
    A_tilde = [eye(n_state),           B_OLTC_x; ...
               zeros(n_oltc, n_state), zeros(n_oltc, n_oltc)];

    % B_tilde (n_state_aug x n_input):
    %   Top-left  B_gen : Delta_V_gen immediately affects [V_load; Q_gen; V_gen]
    %   Top-right 0     : u_OLTC does not immediately affect x (acts via A_tilde next step)
    %   Bottom-left  0  : Delta_V_gen does not write to the u_OLTC storage slot
    %   Bottom-right I  : u_OLTC(k) is stored in x_tilde for propagation at k+1
    B_tilde = [B_gen,               zeros(n_state, n_oltc); ...
               zeros(n_oltc, n_gen), eye(n_oltc)];

    %% Prediction matrices PI and GAMMA
    % A_tilde^n = A_tilde for all n >= 1, so:
    %   PI block i          = A_tilde^(i+1) = A_tilde          for i = 0..N-1
    %   GAMMA block (i,j)   = A_tilde^(i-j) * B_tilde
    %                       = B_tilde           if i == j   (A_tilde^0 = I)
    %                       = A_tilde * B_tilde  if i >  j   (A_tilde^n = A_tilde)
    PI    = repmat(A_tilde, N, 1);                              % (N*n_state_aug x n_state_aug)
    AB    = A_tilde * B_tilde;                                  % precompute once
    GAMMA = kron(eye(N),            B_tilde) + ...              % diagonal blocks
            kron(tril(ones(N), -1), AB);                        % subdiagonal blocks
    % GAMMA is (N*n_state_aug x N*n_input)

    %% QP cost matrices
    % quadprog minimises  0.5*U'*H*U + (f1*x_tilde - f2*Yc)'*U
    CtPsiC = C_tilde' * PSI * C_tilde;    % shared factor  (N*n_state_aug x N*n_state_aug)

    H  = GAMMA' * CtPsiC * GAMMA + H_gamma;
    H  = (H + H') / 2;    % enforce exact symmetry to avoid numerical issues in quadprog

    % f1 (N*n_input x n_state_aug): tracking contribution only
    f1 = GAMMA' * CtPsiC * PI;

    % f2 (N*n_input x N*n_output)
    f2 = GAMMA' * C_tilde' * PSI;

    %% QP constraint matrices
    % Aineq * U <= b1 - b2 * x_tilde
    %
    % b2 (state-dependent RHS component):
    %   Input constraint rows: no x_tilde dependence -> zeros
    %   State constraint rows: Hx*(PI*x_tilde + GAMMA*U) <= bx
    %                          -> Hx*GAMMA*U <= bx - Hx*PI*x_tilde
    %                          -> b2 block: Hx*PI
    n_input_cst_rows = size(Hu, 1);   % 2*N*n_input
    b2 = [zeros(n_input_cst_rows, n_state_aug); ...
          Hx * PI];

    Aineq = [Hu; Hx * GAMMA];

    %% J1: free-response output prediction (for MPC cost logging in main.m)
    % J1 * x_tilde gives the predicted output trajectory when U = 0.
    % Used in main.m: JMPC = J_qp + (J1*x_tilde - Yc)' * PSI * (J1*x_tilde - Yc)
    J1 = C_tilde * PI;    % (N*n_output x n_state_aug)

end
