function [M_all, t_AMG] = computePreconditioners(A, TV0, lagrange)
    % COMPUTEPRECONDITIONERS Precomputes preconditioners for all sequence matrices
    sizeSeq = length(A);
    M_all = cell(sizeSeq, 1); 
    t_AMG = zeros(sizeSeq, 1);
    
    for i = 1:sizeSeq
        [M_all{i}, t_AMG(i)] = buildPrecond(A{i}, TV0{i}, lagrange);
        fprintf('Preconditioner computed for sequence step %d/%d\n', i, sizeSeq);
    end
end

function [M, t_prec] = buildPrecond(A_mat, TV0_vec, lagrange)
    % BUILDPRECOND Helper to compute AMG or RACP preconditioner
    t0 = tic;
    if ~lagrange
       [M, ~] = computeAMG(A_mat, TV0_vec, false, 0);
    else
       [M, ~] = computeRACP(A_mat, TV0_vec, false, 0);
    end
    t_prec = toc(t0);
end