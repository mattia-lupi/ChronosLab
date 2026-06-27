% MEX_sam_adaptive_left.m
function [N, res_norm] = MEX_sam_adaptive_left(Ak, A0, nthread,nstep,step_size,eps)

% Convert the input matrix A to CSC format for efficient processing
[jatk, iak, coefk] = unpack_csc(Ak);
[jat0, ia0, coef0] = unpack_csc(A0);

% Pass directly the csc matrix to use the transposed version of the
% algorithm
[row_N, col_N, val_N] = sam_adaptive_left_mex(...
    jatk, iak, coefk, ...
    jat0, ia0, coef0, ...
    nthread,nstep,step_size,eps);

N = sparse(row_N, col_N, val_N, size(A0,1), size(A0,1));

if nargout > 1
    norm_A0 = norm(A0, 'fro');
    R = N * Ak - A0;
    if norm_A0 > 0
        res_norm = norm(R, 'fro') / norm_A0;
    else
        res_norm = norm(R, 'fro');
    end
end

end