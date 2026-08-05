function params = getParamsAndInit(sizeSeq, t_AMG)
    % GETPARAMSANDINIT Returns solver parameter settings and initial timing matrices
    params.gmres_restart  = 100;
    params.gmres_tol      = 1e-6;
    params.gmres_maxit    = 20;
    params.nstep          = 2;
    params.step_size      = 4;
    params.epss           = 1e-4;
    params.num_strategies = 3;

    % Recomputation step schedules
    params.q  = [4 19 20 26 30 31 36 43 44 49 54 56 57 63 67 69 70 75 79 82 83 88 92 98 99 104 108]; % [cite: 23]
    params.tP = [];
    params.tS = [];

    params.iter    = zeros(sizeSeq, params.num_strategies);
    params.t_setup = zeros(sizeSeq, params.num_strategies);
    params.t_solve = zeros(sizeSeq, params.num_strategies);

    % Pre-assign setup times measured during AMG precomputation for Strategy 1
    params.t_setup(:, 1) = t_AMG(1:sizeSeq);
end