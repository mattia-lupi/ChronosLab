function preproc = sam_preprocess_left(Ak, S)
%SAM_PREPROCESS_LEFT  Preprocessing for the LEFT Sparse Approximate Map.
%
%   preproc = SAM_PREPROCESS_LEFT(Ak, S) analyses the sparsity pattern S
%   and the matrix Ak once to build the index arrays for the n small
%   least-squares problems solved by SAM_COMPUTE_LEFT.
%
%   The LEFT SAM solves the row-wise minimisation
%
%       N_k = argmin_{N with pattern S}  || N * Ak - A0 ||_F         (*)
%
%   which decomposes into n independent row problems.  The updated LEFT
%   preconditioner is then  P_k = P0 * N_k  (compare the RIGHT SAM where
%   P_k = N_k * P0).
%
%   (*) corresponds to eq. (2.13) in Carr, de Sturler & Gugercin (2021):
%       N_k = argmin_{N in S}  || N Ak - A0 ||_F
%
%   Inputs
%   ------
%   Ak  : n-by-n sparse matrix  (source matrix)
%   S   : n-by-n sparse logical / binary matrix defining the desired
%         sparsity pattern of N_k.
%         DEFAULT (omitted or empty): S = spones(Ak)
%         Following Islam, Carr & Jacobs (PAMM 2024), using the pattern of
%         the source matrix Ak is recommended for PDE-derived sequences.
%
%   Output
%   ------
%   preproc : struct with fields
%     .n         - matrix order n
%     .S         - sparsity pattern used
%     .s_idx     - n-by-1 cell: s_idx{i} = column positions of nonzeros
%                  in ROW i of N_k  (i.e. which columns of Ak row i "sees")
%     .r_idx     - n-by-1 cell: r_idx{i} = relevant COLUMN indices for the
%                  i-th row LS problem (union of nonzero columns across
%                  rows s_idx{i} of Ak)
%     .nnz_total - total nonzeros across all rows
%
%   Row-wise LS problem (row i of N_k * Ak ≈ row i of A0)
%   -------------------------------------------------------
%   Let  s_i = s_idx{i}  (nonzero column positions in row i of N)
%        r_i = r_idx{i}  (relevant column indices)
%
%   Then the LS problem for row i is:
%       min_z  || Ak(s_i, r_i).' * z - A0(i, r_i).' ||_2
%   and the solution z fills positions (i, s_i) of N_k.
%
%   Reference
%   ---------
%   Carr, de Sturler & Gugercin, "Preconditioning parametrized linear
%   systems", SIAM J. Sci. Comput. 43(3), A2242-A2267 (2021). Eq. (2.13).
%
%   Islam, Carr & Jacobs, "Optimization of approximate maps for linear
%   systems arising in discretized PDEs", PAMM (2024).
%
%   See also: SAM_PREPROCESS (right version), SAM_COMPUTE_LEFT, SAM_APPLY_LEFT.

n = size(Ak, 1);
assert(size(Ak, 2) == n, 'sam_preprocess_left: Ak must be square.');

% Default sparsity pattern: pattern of Ak
if nargin < 2 || isempty(S)
    S = spones(Ak);
end
S = logical(S);

% ------------------------------------------------------------------
% Build row -> sorted-column-index maps using accumarray
% ------------------------------------------------------------------
[Si, Sj] = find(S);   % Si = row index, Sj = col index
[Ai, Aj] = find(Ak);  % Ai = row, Aj = col  (of Ak)

% S_by_row{i}  = sorted column indices where row i of S (= row i of N) is nonzero
if isempty(Si)
    S_by_row = cell(n, 1);
else
    S_by_row = accumarray(Si(:), Sj(:), [n 1], ...
                          @(v){sort(v)}, {zeros(0,1,'double')});
end

% A_by_row{i}  = sorted nonzero COLUMN indices in row i of Ak
if isempty(Ai)
    A_by_row = cell(n, 1);
else
    A_by_row = accumarray(Ai(:), Aj(:), [n 1], ...
                          @(v){sort(v)}, {zeros(0,1,'double')});
end

% ------------------------------------------------------------------
% For each row i, collect s_i (nonzero cols of N row i)
%                        r_i (relevant cols = union of nonzero cols
%                             of rows s_i of Ak)
% ------------------------------------------------------------------
s_idx = cell(n, 1);
r_idx = cell(n, 1);

for i = 1:n
    si = S_by_row{i};      % columns where row i of N is nonzero
    s_idx{i} = si;

    if isempty(si)
        r_idx{i} = zeros(0, 1, 'double');
        continue;
    end

    % r_i = union of nonzero column indices of each ROW j of Ak, j in s_i
    ri = [];
    for jdx = 1:numel(si)
        ri = union(ri, A_by_row{si(jdx)});
    end
    r_idx{i} = ri(:);
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
