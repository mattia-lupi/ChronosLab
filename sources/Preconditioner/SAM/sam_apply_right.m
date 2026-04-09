function y = sam_apply(N, P0, x)
%SAM_APPLY  Apply the recycled preconditioner  P_k = N_k * P0  to a vector.
%
%   y = SAM_APPLY(N, P0, x)  computes
%
%       y = N_k * (P0 * x)
%
%   i.e. it applies the updated RIGHT preconditioner P_k = N_k * P0 to
%   the vector (or matrix of vectors) x.
%
%   In a right-preconditioned Krylov solver (e.g. GMRES) the preconditioned
%   system is  A_k * P_k * y = b,  with  x_k = P_k * y_k.  Pass this
%   function as the preconditioner application routine, e.g.
%       gmres(@(v)Ak*v, b, restart, tol, maxit, [], @(v)sam_apply(N,P0,v))
%
%   Inputs
%   ------
%   N  : n-by-n sparse SAM returned by SAM_COMPUTE
%   P0 : base preconditioner for A0.  Accepted formats:
%
%        (a) Function handle  @(v) ...
%            Computes P0*v directly.  Most flexible option.
%            Example:  P0 = @(v) U \ (L \ (perm * v));
%
%        (b) Struct with LU factors (output of ilu / lu)
%            Required fields: .L (lower triangular), .U (upper triangular)
%            Optional fields (mutually exclusive):
%              .P  - permutation MATRIX  (from [L,U,P]=lu(A0) or ilu(A0))
%                    Applies  v  <-  U \ (L \ (P * v))
%              .p  - permutation VECTOR  (row permutation)
%                    Applies  v  <-  U \ (L \ v(p))
%              .Q  - right permutation MATRIX (from 4-output lu)
%                    Applies  v  <-  Q * (U \ (L \ (P * v)))
%            If neither .P nor .p is present:
%                    Applies  v  <-  U \ (L \ v)
%
%        (c) Sparse or dense matrix  M
%            Computes P0*v = M*v  (direct multiplication).
%            Use when the approximate inverse is stored explicitly.
%
%   x  : n-by-1 vector, or n-by-m matrix for m simultaneous vectors
%
%   Output
%   ------
%   y  : n-by-1 (or n-by-m) result  y = N_k * (P0 * x)
%
%   Typical usage with MATLAB's ilu
%   --------------------------------
%       setup.type = 'ilutp';  setup.droptol = 1e-3;  setup.milu = 'off';
%       [L, U, P] = ilu(A0, setup);
%       P0.L = L;  P0.U = U;  P0.P = P;
%       preproc = sam_preprocess(A1);          % S = spones(A1)
%       N       = sam_compute(A1, A0, preproc);
%       % Right-preconditioned GMRES:
%       Afun  = @(v) A1 * v;
%       Pfun  = @(v) sam_apply(N, P0, v);
%       [x, flag] = gmres(Afun, b, 50, 1e-8, 200, [], Pfun);
%
%   Reference
%   ---------
%   Carr, de Sturler & Gugercin, "Preconditioning parametrized linear
%   systems", SIAM J. Sci. Comput. 43(3), A2242-A2267 (2021), eq. (2.4):
%       P_k = N_k * P0

% ------------------------------------------------------------------
% Step 1 – apply the BASE preconditioner P0 to x
% ------------------------------------------------------------------
if isa(P0, 'function_handle')
    % -------------------------------------------------------
    % (a) User-supplied function handle: maximum flexibility
    % -------------------------------------------------------
    v = P0(x);

elseif isstruct(P0)
    % -------------------------------------------------------
    % (b) Struct containing LU factors
    % -------------------------------------------------------
    if ~(isfield(P0,'L') && isfield(P0,'U'))
        error('sam_apply: P0 struct must contain fields .L and .U.');
    end
    L = P0.L;
    U = P0.U;

    if isfield(P0, 'Q')
        % Four-output lu:  P*A0*Q = L*U  =>  A0^{-1} = Q * U^{-1} * L^{-1} * P
        if isfield(P0, 'P')
            v = P0.Q * (U \ (L \ (P0.P * x)));
        elseif isfield(P0, 'p')
            v = P0.Q * (U \ (L \ x(P0.p, :)));
        else
            v = P0.Q * (U \ (L \ x));
        end
    elseif isfield(P0, 'P')
        % Three-output lu / ilu with permutation MATRIX P:
        %   P * A0 ≈ L * U  =>  A0^{-1} x ≈ U \ (L \ (P * x))
        v = U \ (L \ (P0.P * x));
    elseif isfield(P0, 'p')
        % Permutation VECTOR p  (row permutation):
        %   A0(p,:) ≈ L * U  =>  A0^{-1} x ≈ U \ (L \ x(p))
        v = U \ (L \ x(P0.p, :));
    else
        % No permutation stored
        v = U \ (L \ x);
    end

else
    % -------------------------------------------------------
    % (c) Explicit matrix (sparse or dense)
    % -------------------------------------------------------
    v = P0 * x;
end

% ------------------------------------------------------------------
% Step 2 – apply the SAM:  y = N_k * v
% ------------------------------------------------------------------
y = N * v;
end
