function preproc = sam_preprocess_right(Ak, S)
%SAM_PREPROCESS  Preprocessing for the Sparse Approximate Map (SAM).
%
%   preproc = SAM_PREPROCESS(Ak, S) analyses the sparsity pattern S and
%   the matrix Ak once to build the index arrays for the n small
%   least-squares (LS) problems solved by SAM_COMPUTE.  Re-run only when
%   either the nonzero PATTERN of Ak or S changes; the numerical VALUES of
%   Ak may change freely between calls to SAM_COMPUTE.
%
%   Inputs
%   ------
%   Ak  : n-by-n sparse matrix  (source matrix, used to find relevant rows)
%   S   : n-by-n sparse logical / binary matrix defining the desired
%         sparsity pattern of the SAM N_k.
%         DEFAULT (omitted or empty): S = spones(Ak), i.e. the sparsity
%         pattern of the current matrix, as recommended by Islam, Carr &
%         Jacobs (PAMM 2024) for PDE-derived sequences.
%
%   Output
%   ------
%   preproc : struct with fields
%     .n      - matrix order n
%     .S      - sparsity pattern used (copy of input S)
%     .s_idx  - n-by-1 cell array; s_idx{k} contains the ROW indices of
%               the nonzeros in column k of N_k  (denoted s_k in the paper)
%     .r_idx  - n-by-1 cell array; r_idx{k} contains the row indices
%               relevant to the k-th LS problem  (denoted r_k in the paper)
%     .nnz_total - total number of nonzeros across all columns (= nnz of N)
%
%   Algorithm
%   ---------
%   For column k of N_k:
%     s_k = { i : (i,k) in S }                  % nonzero rows of col k
%     r_k = union_{j in s_k}  find(Ak(:,j))     % relevant rows for LS
%
%   The LS problem for column k is then
%       min_z  || Ak(r_k, s_k) * z - A0(r_k, k) ||_2
%   and z gives the values at positions (s_k, k) of N_k.
%
%   Reference
%   ---------
%   Carr, de Sturler & Gugercin, "Preconditioning parametrized linear
%   systems", SIAM J. Sci. Comput. 43(3), A2242-A2267 (2021). Alg. 4.1.
%
%   Islam, Carr & Jacobs, "Optimization of approximate maps for linear
%   systems arising in discretized PDEs", PAMM (2024).

n = size(Ak, 1);
assert(size(Ak,2) == n, 'sam_preprocess: Ak must be square.');

% Default sparsity pattern: pattern of Ak (recommended for PDE sequences)
if nargin < 2 || isempty(S)
    S = spones(Ak);
end
S = logical(S);   % ensure logical sparse

% ------------------------------------------------------------------
% Build column -> sorted-row-index maps in one pass using accumarray
% ------------------------------------------------------------------
[Si, Sj] = find(S);
[Ai, Aj] = find(Ak);

% S_by_col{k} = sorted row indices where col k of S (= col k of N) is nonzero
if isempty(Sj)
    S_by_col = cell(n, 1);
else
    S_by_col = accumarray(Sj(:), Si(:), [n 1], ...
                          @(v){sort(v)}, {zeros(0,1,'double')});
end

% A_by_col{j} = sorted nonzero row indices of column j of Ak
if isempty(Aj)
    A_by_col = cell(n, 1);
else
    A_by_col = accumarray(Aj(:), Ai(:), [n 1], ...
                          @(v){sort(v)}, {zeros(0,1,'double')});
end

% ------------------------------------------------------------------
% For each column k, collect s_k and r_k
% ------------------------------------------------------------------
s_idx = cell(n, 1);
r_idx = cell(n, 1);

for k = 1:n
    sk = S_by_col{k};          % potential nonzero rows in col k of N
    s_idx{k} = sk;

    if isempty(sk)
        r_idx{k} = zeros(0,1,'double');
        continue;
    end

    % r_k = union of nonzero rows of each column Ak(:, j), j in s_k
    rk = [];
    for idx = 1:numel(sk)
        rk = union(rk, A_by_col{sk(idx)});
    end
    r_idx{k} = rk(:);
end

% ------------------------------------------------------------------
% Pack results
% ------------------------------------------------------------------
preproc.n         = n;
preproc.S         = S;
preproc.s_idx     = s_idx;
preproc.r_idx     = r_idx;
preproc.nnz_total = sum(cellfun(@numel, s_idx));
end
