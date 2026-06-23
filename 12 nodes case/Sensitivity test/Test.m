% Test.m
% Validates sensitivity matrices (Cv, Cq, Cv_cap, Cq_cap, Cv_tap, Cq_tap)
% by comparing analytical values 
% against finite-difference estimates obtained by perturbing the base case.
%
% For each sensitivity matrix X, the relative error is computed as:
%   (X_fd - X_analytical) / (X_fd + epsilon) * 100   [%]
% where epsilon avoids division by zero on near-zero entries.
%
% Sensitivity definitions:
%   Cv(i,j)      : dV_load(i) / dVg_setpoint(j)       [pu/pu]
%   Cq(i,j)      : dQg(i)     / dVg_setpoint(j)       [pu/pu]
%   Cv_cap(i,j)  : dV_load(i) / dBsh(j)               [pu/pu]
%   Cq_cap(i,j)  : dQg(i)     / dBsh(j)               [pu/pu]
%   Cv_tap(i,j)  : dV_load(i) / dtap(j)               [pu/pu]
%   Cq_tap(i,j)  : dQg(i)     / dtap(j)               [pu/pu]

clear all

%% Network setup
mpc = case9xx_Bsh();
define_constants;
mpopt = mpoption('verbose', 0, 'out.all', 0, 'pf.enforce_q_lims', 0);

baseMVA = mpc.baseMVA;

% Bus index classification
load_idx         = find(mpc.bus(:, BUS_TYPE) == PQ);          % PQ load buses
pv_and_slack_idx = find(ismember(mpc.bus(:, BUS_TYPE), [PV, REF])); % generator and slack buses
target_idx       = find(mpc.bus(:, PD) ~= 0);                % buses with non-zero active load

% Branch indices of tap-changing transformers (buses 10-5 and 11-6)
transformer_indices = [16; 17];

% Dimensions
n_bus    = size(mpc.bus, 1);
n_load   = size(load_idx, 1);
n_gen    = size(mpc.gen, 1);
n_branch = size(mpc.branch, 1);
n_transfo = size(transformer_indices, 1);
n_target = size(target_idx, 1);

% Capacitor susceptance value
cap = 5;

% Map each bus to its row index in mpc.gen (0 if no generator on that bus)
gen_number = zeros(n_bus, 1);
for i = 1:n_gen
    gen_number(mpc.gen(i, 1)) = i;
end

%% Base-case power flow and Jacobian
results0 = runpf(mpc, mpopt);
JAC      = makeJac(results0, true); % full Jacobian (angles and magnitudes)

%% Analytical sensitivity matrices
[Cv,     Cq    ] = Cal_CvCq(        mpc, JAC);
[Cv_cap, Cq_cap] = Cal_Cv_capCq_cap(mpc, results0.bus(load_idx, VM), cap, JAC);
[Cv_tap, Cq_tap] = Cal_Cv_tapCq_tap(mpc, results0.bus(load_idx, VM), transformer_indices, JAC);

%% -----------------------------------------------------------------------
%  Finite-difference validation
%  Perturbation: +0.1 pu on generator voltage setpoint
%  Finite-difference estimate: (x_perturbed - x_base) / 0.1
% -----------------------------------------------------------------------

dVg = 0.1; % generator voltage perturbation [pu]

disp('Percentage of error between finite-difference and analytical sensitivities')

%% Cv  (load voltage / generator voltage setpoint)
x0 = results0.bus(load_idx, VM);

Cv_fd_temp = zeros(n_load, n_gen);
for i = 1:n_gen
    mpc1 = mpc;
    mpc1.gen(i, VG) = mpc1.gen(i, VG) + dVg;
    results1 = runpf(mpc1, mpopt);
    Cv_fd_temp(:, i) = (results1.bus(load_idx, VM) - x0) / dVg;
end

Cv_fd = sparse(n_bus, n_bus);
for i = 1:n_load
    for j = 1:n_gen
        Cv_fd(load_idx(i), pv_and_slack_idx(j)) = Cv_fd_temp(i, j);
    end
end

disp('--------------------------')
disp('Cv')
disp((Cv_fd - Cv) ./ (Cv_fd + 1e-9) * 100)

%% Cq  (generator reactive power / generator voltage setpoint)
x0 = results0.gen(gen_number(pv_and_slack_idx), QG) / baseMVA;

Cq_fd_temp = zeros(n_gen, n_gen);
for i = 1:n_gen
    mpc1 = mpc;
    mpc1.gen(i, VG) = mpc1.gen(i, VG) + dVg;
    results1 = runpf(mpc1, mpopt);
    Cq_fd_temp(:, i) = (results1.gen(gen_number(pv_and_slack_idx), QG) / baseMVA - x0) / dVg;
end

Cq_fd = sparse(n_bus, n_bus);
for i = 1:n_gen
    for j = 1:n_gen
        Cq_fd(pv_and_slack_idx(i), pv_and_slack_idx(j)) = Cq_fd_temp(i, j);
    end
end

disp('--------------------------')
disp('Cq')
disp((Cq_fd - Cq) ./ (Cq_fd + 1e-9) * 100)

%% Cv_cap  (load voltage / shunt susceptance step)
% Perturbation: adding one capacitor
x0 = results0.bus(load_idx, VM);

Cv_cap_fd_temp = zeros(n_load, n_target);
for i = 1:n_target
    mpc1 = mpc;
    mpc1.bus(target_idx(i), BS) = mpc1.bus(target_idx(i), BS) + cap;
    results1 = runpf(mpc1, mpopt);
    Cv_cap_fd_temp(:, i) = (results1.bus(load_idx, VM) - x0);
end

Cv_cap_fd = sparse(n_bus, n_bus);
for i = 1:n_load
    for j = 1:n_target
        Cv_cap_fd(load_idx(i), target_idx(j)) = Cv_cap_fd_temp(i, j);
    end
end

disp('--------------------------')
disp('Cv_cap')
disp((Cv_cap_fd - Cv_cap) ./ (Cv_cap_fd + 1e-9) * 100)

%% Cq_cap  (generator reactive power / shunt susceptance step)
x0 = results0.gen(gen_number(pv_and_slack_idx), QG) / baseMVA;

Cq_cap_fd_temp = zeros(n_gen, n_target);
for i = 1:n_target
    mpc1 = mpc;
    mpc1.bus(target_idx(i), BS) = mpc1.bus(target_idx(i), BS) + cap;
    results1 = runpf(mpc1, mpopt);
    Cq_cap_fd_temp(:, i) = (results1.gen(gen_number(pv_and_slack_idx), QG) / baseMVA - x0);
end

Cq_cap_fd = sparse(n_bus, n_bus);
for i = 1:n_gen
    for j = 1:n_target
        Cq_cap_fd(pv_and_slack_idx(i), target_idx(j)) = Cq_cap_fd_temp(i, j);
    end
end

disp('--------------------------')
disp('Cq_cap')
disp((Cq_cap_fd - Cq_cap) ./ (Cq_cap_fd + 1e-9) * 100)

%% Cv_tap  (load voltage / transformer tap ratio)
% Perturbation: one discrete tap step on each transformer
x0        = results0.bus(load_idx, VM);
tap_step  = 0.00625; % one tap step = 0.625 % of nominal ratio

Cv_tap_fd_temp = zeros(n_load, n_transfo);
for i = 1:n_transfo
    mpc1 = mpc;
    mpc1.branch(transformer_indices(i), TAP) = mpc1.branch(transformer_indices(i), TAP) + tap_step;
    results1 = runpf(mpc1, mpopt);
    Cv_tap_fd_temp(:, i) = (results1.bus(load_idx, VM) - x0) / tap_step;
end

Cv_tap_fd = sparse(n_bus, n_branch);
for i = 1:n_load
    for j = 1:n_transfo
        Cv_tap_fd(load_idx(i), transformer_indices(j)) = Cv_tap_fd_temp(i, j);
    end
end

disp('--------------------------')
disp('Cv_tap')
disp((Cv_tap_fd - Cv_tap) ./ (Cv_tap_fd + 1e-9) * 100)

%% Cq_tap  (generator reactive power / transformer tap ratio)
x0 = results0.gen(gen_number(pv_and_slack_idx), QG) / baseMVA;

Cq_tap_fd_temp = zeros(n_gen, n_transfo);
for i = 1:n_transfo
    mpc1 = mpc;
    mpc1.branch(transformer_indices(i), TAP) = mpc1.branch(transformer_indices(i), TAP) + tap_step;
    results1 = runpf(mpc1, mpopt);
    Cq_tap_fd_temp(:, i) = (results1.gen(gen_number(pv_and_slack_idx), QG) / baseMVA - x0) / tap_step;
end

Cq_tap_fd = sparse(n_bus, n_branch);
for i = 1:n_gen
    for j = 1:n_transfo
        Cq_tap_fd(pv_and_slack_idx(i), transformer_indices(j)) = Cq_tap_fd_temp(i, j);
    end
end

disp('--------------------------')
disp('Cq_tap')
disp((Cq_tap_fd - Cq_tap) ./ (Cq_tap_fd + 1e-9) * 100)
