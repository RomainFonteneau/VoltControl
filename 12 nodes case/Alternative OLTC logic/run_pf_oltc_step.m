function [mpc_out, results, final_taps] = run_pf_oltc_step(mpc_in, target_v, idx_list)
    % RUN_PF_OLTC_STEP Moves taps by exactly ONE step per call.
    % each 10 seconds
    % mpc_in: The current MATPOWER case struct.
    % target_v: Desired voltage magnitude (pu) for regulated buses.
    % idx_list: Vector of branch indices to be controlled.

    define_constants;
    mpc = mpc_in; 
    
    % --- Configuration ---
    step_size = 0.00625; % 0.625% fixed step size
    tol = 0.005;         % Voltage deadband
        
    % Identify regulated buses (LV bus)
    reg_buses = mpc.branch(idx_list, T_BUS);
    % 1 HV 2 LV, mandatory to Z=constant in pu!
    % but ratio=V1/V2 in matpower

    % --- Run Power Flow to check current state ---
    mpopt = mpoption('verbose', 0,'out.all', 0,'pf.enforce_q_lims', 1);
    results = runpf(mpc, mpopt);
    v_actuals = results.bus(reg_buses, VM);
    diffs = v_actuals - target_v;
    
    % --- Single Step Logic ---
    for k = 1:length(idx_list)
        if abs(diffs(k)) > tol
            if diffs(k) < 0
                % Voltage low -> Decrease Ratio on LV side to boost voltage
                % V2(LV) = V1(HV)/ratio
                mpc.branch(idx_list(k), TAP) = mpc.branch(idx_list(k), TAP) - step_size;
            else
                % Voltage high -> Increase Ratio on LV side to lower voltage
                mpc.branch(idx_list(k), TAP) = mpc.branch(idx_list(k), TAP) + step_size;
            end
            fprintf('trf %d , initial tap %6.3f\n', idx_list(k), mpc.branch(idx_list(k), TAP) );
        end
    end
    
    % Enforce physical limits [0.9, 1.1]
    mpc.branch(idx_list, TAP) = max(0.9, min(1.1, mpc.branch(idx_list, TAP)));
    
    % Taps will change at the next runpf!
    mpc_out = mpc;
    final_taps = mpc.branch(idx_list, TAP);
end