function y = sam_apply_left(N, P0, x)
%SAM_APPLY_LEFT  Apply the recycled LEFT preconditioner  P_k = P0 * N_k.
%
%   y = SAM_APPLY_LEFT(N, P0, x)  computes
%
%       y = P0 * (N_k * x)
%
%   i.e. it applies the updated LEFT preconditioner P_k = P0 * N_k to the
%   vector (or matrix of vectors) x.
%
%   Contrast with SAM_APPLY (right version), which computes y = N_k*(P0*x)
%   for the RIGHT preconditioner P_k = N_k * P0.
%
%   In a LEFT-preconditioned Krylov solver the preconditioned system is
%
%       (P_k * Ak) x = P_k * b
%
%   so the solver must apply P_k to arbitrary vectors during the iteration.
%   Pass this function as the LEFT preconditioner in MATLAB's gmres:
%
%       Pfun = @(v) sam_apply_left(N, P0, v);
%       [x, flag] = gmres(@(v)Ak*v, b, restart, tol, maxit, Pfun);
%
%   Note: for left preconditioning, the preconditioned RHS P_k*b must also
%   be computed before calling gmres if the solver does not do it internally.
%   MATLAB's gmres handles this automatically when Pfun is passed as M1.
%
%   Inputs
%   ------
%   N  : n-by-n sparse LEFT SAM returned by SAM_COMPUTE_LEFT
%        (satisfies N * Ak ≈ A0)
%   P0 : base preconditioner for A0.  Accepted formats:
%
%        (a) Function handle  @(v) ...      — computes P0*v directly
%
%        (b) Struct with fields .L and .U (and optionally .P / .p / .Q)
%              .P  permutation MATRIX:  y = U\(L\(P*v))
%              .p  permutation VECTOR:  y = U\(L\v(p))
%              .Q  right permutation:   y = Q*(U\(L\(P*v)))
%              (neither): y = U\(L\v)
%
%        (c) Sparse or dense matrix M — computes P0*v = M*v
%
%   x  : n-by-1 vector, or n-by-m matrix of m simultaneous vectors
%
%   Output
%   ------
%   y  : n-by-1 (or n-by-m) result  y = P0 * (N_k * x)
%
%   Typical usage
%   -------------
%       opts.type    = 'ilutp';
%       opts.droptol = 1e-4;
%       [L, U, P]    = ilu(A0, opts);
%       P0.L = L;  P0.U = U;  P0.P = P;
%
%       preproc = sam_preprocess_left(A1);         % S = spones(A1)
%       N       = sam_compute_left(A1, A0, preproc);
%
%       % Left-preconditioned GMRES:
%       Afun = @(v) A1 * v;
%       Pfun = @(v) sam_apply_left(N, P0, v);
%       [x, flag] = gmres(Afun, b, 50, 1e-8, 200, Pfun);
%
%   Reference
%   ---------
%   Carr, de Sturler & Gugercin, "Preconditioning parametrized linear
%   systems", SIAM J. Sci. Comput. 43(3), A2242-A2267 (2021), eq. (2.13):
%       P_k = P0 * N_k

% ------------------------------------------------------------------
% Step 1 – apply the LEFT SAM N_k to x:  v = N_k * x
% ------------------------------------------------------------------
v = N * x;

% ------------------------------------------------------------------
% Step 2 – apply the BASE preconditioner P0 to v:  y = P0 * v
% ------------------------------------------------------------------
if isa(P0, 'function_handle')
    % -------------------------------------------------------
    % (a) User-supplied function handle
    % -------------------------------------------------------
    y = P0(v);

elseif isstruct(P0)
    % -------------------------------------------------------
    % (b) Struct containing LU factors from ilu / lu
    % -------------------------------------------------------
    if ~(isfield(P0,'L') && isfield(P0,'U'))
        error('sam_apply_left: P0 struct must contain fields .L and .U.');
    end
    L = P0.L;
    U = P0.U;

    if isfield(P0, 'Q')
        if isfield(P0, 'P')
            y = P0.Q * (U \ (L \ (P0.P * v)));
        elseif isfield(P0, 'p')
            y = P0.Q * (U \ (L \ v(P0.p, :)));
        else
            y = P0.Q * (U \ (L \ v));
        end
    elseif isfield(P0, 'P')
        y = U \ (L \ (P0.P * v));
    elseif isfield(P0, 'p')
        y = U \ (L \ v(P0.p, :));
    else
        y = U \ (L \ v);
    end

else
    % -------------------------------------------------------
    % (c) Explicit matrix
    % -------------------------------------------------------
    y = P0 * v;
end
end
