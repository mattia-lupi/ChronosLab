function [N, res_norm] = sam_compute_left(Ak, A0, preproc)
%SAM_COMPUTE_LEFT  Compute the LEFT Sparse Approximate Map (SAM) N_k.
%
%   N = SAM_COMPUTE_LEFT(Ak, A0, preproc) solves
%
%       N_k = argmin_{N with sparsity pattern preproc.S}  || N*Ak - A0 ||_F
%
%   by decomposing the global minimisation into n independent row-wise
%   least-squares problems.  Each is solved with MATLAB's backslash.
%
%   [N, res_norm] = SAM_COMPUTE_LEFT(...) also returns the relative
%   residual  || N*Ak - A0 ||_F / || A0 ||_F.
%
%   Inputs
%   ------
%   Ak      : n-by-n sparse source matrix  ("new" matrix A_k)
%   A0      : n-by-n sparse target matrix  (whose preconditioner P0 we recycle)
%   preproc : struct from SAM_PREPROCESS_LEFT(Ak, S).
%             Rebuild if the nonzero PATTERN of Ak changes between steps.
%
%   Output
%   ------
%   N        : n-by-n sparse matrix with pattern preproc.S such that
%              N * Ak ≈ A0  (the left sparse approximate map)
%   res_norm : scalar,  || N*Ak - A0 ||_F / || A0 ||_F  (optional)
%
%   Updated preconditioner
%   ----------------------
%   The updated LEFT preconditioner is  P_k = P0 * N_k.
%   Apply it with SAM_APPLY_LEFT(N, P0, x).
%   In left-preconditioned GMRES solve:  (P_k * Ak) x = P_k * b
%
%   Row-wise LS problem (row i)
%   ----------------------------
%   For row i the problem is:
%       min_z  || Ak(s_i, r_i).' * z - A0(i, r_i).' ||_2
%   where:
%     s_i = preproc.s_idx{i}  (column positions of nonzeros in row i of N)
%     r_i = preproc.r_idx{i}  (relevant column indices)
%   The solution z gives values at (i, s_i) in N.
%
%   This is the transpose of the column problem in SAM_COMPUTE (right),
%   because (N Ak)_{i,:} = n_i^T Ak, and restricting to columns r_i gives
%       Ak(s_i, r_i).' * z ≈ A0(i, r_i).'
%
%   Reference
%   ---------
%   Carr, de Sturler & Gugercin, "Preconditioning parametrized linear
%   systems", SIAM J. Sci. Comput. 43(3), A2242-A2267 (2021). Eq. (2.13).
%
%   See also: SAM_PREPROCESS_LEFT, SAM_APPLY_LEFT, SAM_COMPUTE (right version).

n         = preproc.n;
s_idx     = preproc.s_idx;
r_idx     = preproc.r_idx;
nnz_total = preproc.nnz_total;

% ------------------------------------------------------------------
% Preallocate COO storage
% ------------------------------------------------------------------
row_N = zeros(nnz_total, 1);
col_N = zeros(nnz_total, 1);
val_N = zeros(nnz_total, 1);

cnt = 0;

% ------------------------------------------------------------------
% Main loop: one row LS problem per row i
% ------------------------------------------------------------------
for i = 1:n
    si   = s_idx{i};    % nonzero column positions in row i of N
    ri   = r_idx{i};    % relevant column indices for LS
    nnzi = numel(si);

    if nnzi == 0
        continue;
    end

    if isempty(ri)
        z = zeros(nnzi, 1);
    else
        % Row i of N * Ak restricted to relevant columns ri:
        %   n_i(si).' * Ak(si, ri)  ≈  A0(i, ri)
        %
        % Taking the transpose:
        %   Ak(si, ri).' * z  ≈  A0(i, ri).'
        %
        % where z = n_i(si) (column vector of unknowns).

        Atmp = full(Ak(si, ri)).';   % size: |r_i| x |s_i|
        f    = full(A0(i,  ri)).';   % size: |r_i| x 1

        warning('off', 'MATLAB:rankDeficientMatrix');
        warning('off', 'MATLAB:nearlySingularMatrix');
        z = Atmp \ f;
        warning('on',  'MATLAB:rankDeficientMatrix');
        warning('on',  'MATLAB:nearlySingularMatrix');

        z(~isfinite(z)) = 0;
    end

    % Store: nonzero at (i, si) with value z
    row_N(cnt+1 : cnt+nnzi) = i;
    col_N(cnt+1 : cnt+nnzi) = si;
    val_N(cnt+1 : cnt+nnzi) = z;
    cnt = cnt + nnzi;
end

% ------------------------------------------------------------------
% Assemble sparse matrix in one call
% ------------------------------------------------------------------
N = sparse(row_N(1:cnt), col_N(1:cnt), val_N(1:cnt), n, n);

% ------------------------------------------------------------------
% Optional: relative residual  || N*Ak - A0 ||_F / || A0 ||_F
% ------------------------------------------------------------------
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
