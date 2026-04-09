function [N, res_norm] = sam_compute(Ak, A0, preproc)
%SAM_COMPUTE  Compute the Sparse Approximate Map (SAM) N_k.
%
%   N = SAM_COMPUTE(Ak, A0, preproc) solves
%
%       N_k = argmin_{N with sparsity pattern preproc.S}  || Ak*N - A0 ||_F
%
%   by decomposing the global minimization into n independent small
%   least-squares problems (one per column) which are solved with
%   MATLAB's backslash operator.
%
%   [N, res_norm] = SAM_COMPUTE(...) additionally returns the relative
%   residual  || Ak*N - A0 ||_F / || A0 ||_F.
%
%   Inputs
%   ------
%   Ak      : n-by-n sparse source matrix  ("new" matrix A_k)
%   A0      : n-by-n sparse target matrix  (matrix whose preconditioner P0
%             we wish to recycle)
%   preproc : struct from SAM_PREPROCESS(Ak, S)
%             Must have been built with the SAME sparsity pattern of Ak.
%             Re-run SAM_PREPROCESS if the nonzero pattern of Ak changes.
%
%   Output
%   ------
%   N        : n-by-n sparse matrix with sparsity pattern preproc.S such
%              that  Ak * N ≈ A0  (the sparse approximate map)
%   res_norm : scalar,  || Ak*N - A0 ||_F / || A0 ||_F  (optional)
%
%   Updated preconditioner
%   ----------------------
%   The updated right preconditioner is  P_k = N_k * P0.
%   Apply it with SAM_APPLY(N, P0, x).
%
%   Algorithmic notes
%   -----------------
%   For column k the LS problem is
%       min_z  || Ak(r_k, s_k) * z - A0(r_k, k) ||_2
%   where s_k and r_k come from preproc.  The solution z fills the nonzero
%   values of column k of N.  The result is assembled in COO format and
%   converted to MATLAB sparse at the end (one sparse() call).
%
%   If the system is overdetermined (typical), backslash uses QR.
%   If it is underdetermined (S denser than Ak), the minimum-norm solution
%   is returned (also via backslash / QR with column pivoting).
%
%   Reference
%   ---------
%   Carr, de Sturler & Gugercin, "Preconditioning parametrized linear
%   systems", SIAM J. Sci. Comput. 43(3), A2242-A2267 (2021). Alg. 4.2.

n         = preproc.n;
s_idx     = preproc.s_idx;
r_idx     = preproc.r_idx;
nnz_total = preproc.nnz_total;

% ------------------------------------------------------------------
% Preallocate COO (coordinate) storage
% ------------------------------------------------------------------
row_N = zeros(nnz_total, 1);
col_N = zeros(nnz_total, 1);
val_N = zeros(nnz_total, 1);

cnt = 0;   % running nonzero counter

% ------------------------------------------------------------------
% Main loop: solve one small LS problem per column
% ------------------------------------------------------------------
for k = 1:n
    sk   = s_idx{k};    % row positions of nonzeros in col k of N
    rk   = r_idx{k};    % relevant row indices for this LS problem
    nnzk = numel(sk);

    if nnzk == 0
        continue;        % column k of N is structurally zero – skip
    end

    % Extract small dense submatrix  Atmp = Ak(rk, sk)  and RHS
    if isempty(rk)
        % No relevant rows: set column to zero (cannot improve residual)
        z = zeros(nnzk, 1);
    else
        Atmp = full(Ak(rk, sk));   % size: |r_k| x |s_k|
        f    = full(A0(rk, k));    % size: |r_k| x 1

        % Solve  min_z || Atmp*z - f ||_2
        % Backslash handles over- and under-determined cases via QR.
        % Suppress the rank-deficiency warning for near-singular subproblems.
        warning('off', 'MATLAB:rankDeficientMatrix');
        warning('off', 'MATLAB:nearlySingularMatrix');
        z = Atmp \ f;
        warning('on',  'MATLAB:rankDeficientMatrix');
        warning('on',  'MATLAB:nearlySingularMatrix');

        % Replace any NaN/Inf (degenerate subproblem) with zeros
        z(~isfinite(z)) = 0;
    end

    % Store in COO format (row index = s_k, col index = k, value = z)
    row_N(cnt+1 : cnt+nnzk) = sk;
    col_N(cnt+1 : cnt+nnzk) = k;
    val_N(cnt+1 : cnt+nnzk) = z;
    cnt = cnt + nnzk;
end

% ------------------------------------------------------------------
% Assemble sparse output (single sparse() call for efficiency)
% ------------------------------------------------------------------
N = sparse(row_N(1:cnt), col_N(1:cnt), val_N(1:cnt), n, n);

% ------------------------------------------------------------------
% Optional: relative residual  || Ak*N - A0 ||_F / || A0 ||_F
% ------------------------------------------------------------------
if nargout > 1
    norm_A0 = norm(A0, 'fro');
    R = Ak * N - A0;
    if norm_A0 > 0
        res_norm = norm(R, 'fro') / norm_A0;
    else
        res_norm = norm(R, 'fro');
    end
end
end
