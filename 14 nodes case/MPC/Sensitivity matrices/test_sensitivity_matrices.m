% Validates the six sensitivity matrices Cv,Cq,Cv_tap,Cq_tap,Cv_cap,Cq_cap
% by comparing them against finite-difference approximations.
%
% For each generator, the voltage setpoint is perturbed by +dV (pu),
% a new power flow is run, and the resulting changes in load voltages and
% generator reactive powers are divided by dV to obtain numerical sensitivities.
%
% The same approach is applied for tap and capacitors, the tap is modified,
% a capacitor is added.
%
% For each sensitivity matrix, the script displays the matrix of element-wise
% percentage differences: 100 * (analytical - numerical) / |numerical|.
% A result near zero everywhere confirms that the analytical formula is correct.

clear all
define_constants;

% -------------------------------------------------------------------------
% Configuration
% -------------------------------------------------------------------------
mpc      = Extented_case();        % MATPOWER case
dV       = 0.01;       % Voltage perturbation magnitude (pu)
cap      = 0.05;       % Susceptance of the added capacitor (pu)
dTAP     = 0.00625;    % Tap modification

% -------------------------------------------------------------------------
% Base-case power flow and Jacobian
% -------------------------------------------------------------------------
mpopt = mpoption('verbose', 0, 'out.all', 0, 'pf.enforce_q_lims', 0);

results0 = runpf(mpc, mpopt);

voltage = results0.bus(:,VM);

JAC = makeJac(results0,true);

baseMVA=mpc.baseMVA;

% Bus classification (internal numbering)
pv_and_slack_idx = find(ismember(mpc.bus(:, BUS_TYPE), [PV, REF]));
pv_idx           = find(mpc.bus(:, BUS_TYPE) == PV);
pq_idx           = find(mpc.bus(:, BUS_TYPE) == PQ);
slack_idx        = find(mpc.bus(:, BUS_TYPE) == REF);
oltc_idx         = [5;6;20;21;22;23]; 

n_bus            = size(mpc.bus, 1);
n_pv             = size(pv_idx, 1);
n_slack          = size(slack_idx, 1);
n_gen            = size(mpc.gen,1);
n_oltc           = size(oltc_idx,1);
n_branch         =size(mpc.branch,1);

% gen_number(i) = row index in mpc.gen of the generator at bus i (0 if none)
gen_number = zeros(n_bus, 1);
for i = 1:n_gen
    gen_number(mpc.gen(i,1)) = i;
end

% -------------------------------------------------------------------------
% Analytical sensitivities (Jacobian-based)
% -------------------------------------------------------------------------
[Cv_ana,Cq_ana] = Cal_CvCq(mpc,JAC);
[Cv_tap_ana,Cq_tap_ana] = Cal_Cv_tapCq_tap(mpc,JAC,voltage);
[Cv_cap_ana,Cq_cap_ana] = Cal_Cv_capCq_cap(mpc,JAC,voltage);


Vl0    = results0.bus(pq_idx, VM);             % pq bus voltages at base
Qg0    = results0.gen(gen_number(pv_and_slack_idx),QG);  % gen reactive powers at base
% -------------------------------------------------------------------------
% Numerical Cv and Cq
% One power flow per generator bus: perturb Vg setpoint by +dV
% -------------------------------------------------------------------------

Cv_num = sparse(n_bus, n_bus);
Cq_num = sparse(n_bus, n_bus);

for j = 1:n_pv+n_slack
    gen_j   = gen_number(pv_and_slack_idx(j));
    mpc_j  = results0;

    mpc_j.gen(gen_j, VG) = mpc_j.gen(gen_j, VG) + dV;
    
    res_j = runpf(mpc_j, mpopt);

    dVpq = res_j.bus(pq_idx, VM) - Vl0;
    Cv_num(pq_idx, pv_and_slack_idx(j)) = dVpq / dV;

    dQg = res_j.gen(gen_number(pv_and_slack_idx), QG)/baseMVA - Qg0/baseMVA;
    Cq_num(pv_and_slack_idx, pv_and_slack_idx(j)) = dQg / dV;
end

Cv_num(abs(Cv_num)<0.01)=0;
Cq_num(abs(Cq_num)<0.01)=0;

% -------------------------------------------------------------------------
% Numerical Cv_tap and Cq_tap
% -------------------------------------------------------------------------

Cv_tap_num = sparse(n_bus, n_branch);
Cq_tap_num = sparse(n_bus, n_branch);

for j = 1:n_oltc
    branch_j  = oltc_idx(j);
    mpc_j = results0;

    mpc_j.branch(branch_j, TAP) = mpc_j.branch(branch_j,TAP) +dTAP;

    res_j = runpf(mpc_j, mpopt);

    dVpq = res_j.bus(pq_idx, VM) - Vl0;
    Cv_tap_num(pq_idx, branch_j) = dVpq /dTAP;

    dQg = res_j.gen(gen_number(pv_and_slack_idx), QG)/baseMVA - Qg0/baseMVA;
    Cq_tap_num(pv_and_slack_idx, branch_j)   = dQg / dTAP;
end

Cv_tap_num(abs(Cv_tap_num)<0.01)=0;
Cq_tap_num(abs(Cq_tap_num)<0.01)=0;

% -------------------------------------------------------------------------
% Numerical Cv_cap and Cq_cap
% -------------------------------------------------------------------------

Cv_cap_num = sparse(n_bus, n_bus);
Cq_cap_num = sparse(n_bus, n_bus);

for j = 1:n_bus
    bus_j  = j;
    mpc_j = results0;

    mpc_j.bus(bus_j, BS) = mpc_j.bus(bus_j,BS) + cap * baseMVA;

    res_j = runpf(mpc_j, mpopt);

    dVpq = res_j.bus(pq_idx, VM) - Vl0;
    Cv_cap_num(pq_idx, j) = dVpq / cap;

    dQg = res_j.gen(gen_number(pv_and_slack_idx), QG)/baseMVA - Qg0/baseMVA;
    Cq_cap_num(pv_and_slack_idx, j)   = dQg / cap;
end

Cv_cap_num(abs(Cv_cap_num)<0.01)=0;
Cq_cap_num(abs(Cq_cap_num)<0.01)=0;

% -------------------------------------------------------------------------
% Display percentage differences (analytical vs numerical)
% Rows and columns are restricted to the buses where each matrix is non-trivial
% % -------------------------------------------------------------------------
fprintf('\n=== Cv : %%diff  (rows = load buses, cols = gen+slack buses) ===\n');
disp(pct_diff( ...
    Cv_ana, ...
         Cv_num));

fprintf('\n=== Cq : %%diff  (rows = gen+slack buses, cols = gen+slack buses) ===\n');
disp(pct_diff( ...
    Cq_ana, ...
         Cq_num));

fprintf('\n=== Cv_tap : %%diff  (rows = load buses, cols = load buses) ===\n');
disp(pct_diff( ...
    Cv_tap_ana(:,oltc_idx), ...
         Cv_tap_num(:,oltc_idx)));

fprintf('\n=== Cq_tap : %%diff  (rows = PV buses, cols = load buses) ===\n');
disp(pct_diff( ...
    Cq_tap_ana(:,oltc_idx), ...
         Cq_tap_num(:,oltc_idx)));

fprintf('\n=== Cv_cap : %%diff  (rows = load buses, cols = all buses) ===\n');
disp(pct_diff( ...
    Cv_cap_ana, ...
    Cv_cap_num));

fprintf('\n=== Cq_cap : %%diff  (rows = PV buses, cols = all buses) ===\n');
disp(pct_diff( ...
    Cq_cap_ana, ...
    Cq_cap_num));


% =========================================================================
% Local functions
% =========================================================================

function pct = pct_diff(A, B)
% Element-wise percentage difference of A relative to B:
%   pct(i,j) = 100 * (A(i,j) - B(i,j)) / |B(i,j)|
% Returns 0 where both values are zero (no error), and Inf where B is zero
% but A is not (undefined relative error).
    pct               = 100 * (A - B) ./ abs(B);
    pct(A == 0 & B == 0) = 0;
end
