% MEX_sam_compute_left.m
function [N, res_norm] = MEX_sam_compute_left(Ak, A0, preproc)
% Computes N via C++ MEX
[iatk,jak,coefk] = unpack_csr(Ak);
[iat0,ja0,coef0] = unpack_csr(A0);

% Flat interface bypasses StructArray unpacking in C++
[row_N, col_N, val_N] = sam_compute_left_mex(...
    iatk, jak, coefk, ...
    iat0, ja0, coef0, ...
    int32(preproc.s_ptr), ...
    int32(preproc.s_data), ...
    int32(preproc.r_ptr), ...
    int32(preproc.r_data), ...
    int32(preproc.nnz_total));

N = sparse(row_N, col_N, val_N, preproc.n, preproc.n);

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