% MEX_sam_adaptive_left.m
function [N, res_norm] = MEX_sam_adaptive_left(Ak, A0, nthread,nstep,step_size,eps)

% Check to not break the algorithm
if nstep*step_size < 1
   N = speye(size(A,1));
   res_norm = normSAM(Ak,A0);
   return;
end

% Convert the input matrix A to CSR format for efficient processing
[jatk, iak, coefk] = unpack_csc(Ak);
% Prepare also the transposed to speed up everything
[~, ~, coefkT]     = unpack_csr(Ak);

% Convert the input matrix A to CSC format for efficient processing
[jat0, ia0, coef0] = unpack_csr(A0);

[row_N, col_N, val_N, res_norm] = sam_adaptive_left_mex(...
    jatk, iak, coefk, coefkT, ...
    jat0, ia0, coef0, ...
    nthread,nstep,step_size,eps);

% Get the trasposed one
N = sparse(col_N,row_N, val_N, size(A0,1), size(A0,1));

end