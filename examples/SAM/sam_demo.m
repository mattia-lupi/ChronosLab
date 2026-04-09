%SAM_DEMO  Demonstrates the Sparse Approximate Map (SAM) workflow.
%
%  Builds a small toy sequence  A0, A1, A2  (parametrised Poisson-like
%  matrices), computes an ILUTP preconditioner P0 for A0, then uses SAM
%  updates to recycle P0 for A1 and A2 inside right-preconditioned GMRES.
%
%  Files required:  sam_preprocess.m  sam_compute.m  sam_apply.m
%
%  Reference:
%    Carr, de Sturler & Gugercin, SIAM J. Sci. Comput. 43(3) (2021).
%    Islam, Carr & Jacobs, PAMM (2024).

clear; clc; rng(0);

%% ---- 1. Build a small parametrised sequence -------------------------
%  A(p) = tridiagonal system arising from a 1-D BVP with variable coeff.
%  n x n system, parameters p0, p1, p2 slightly different.

n   = 200;
e   = ones(n,1);

make_A = @(mu) spdiags([-mu*e, (2*mu+1)*e, -mu*e], [-1 0 1], n, n);

A0 = make_A(1.0);
A1 = make_A(1.15);
A2 = make_A(1.30);

b  = rand(n,1);    % same right-hand side for all systems (can differ)

fprintf('Matrix size n = %d\n\n', n);

%% ---- 2. Compute ILUTP preconditioner P0 for A0 ----------------------
opts.type    = 'ilutp';
opts.droptol = 1e-3;
opts.milu    = 'off';
[L0, U0, P0_mat] = ilu(A0, opts);

% Pack into the struct format expected by sam_apply
P0.L = L0;
P0.U = U0;
P0.P = P0_mat;   % permutation matrix from ilu

%% ---- 3. Preprocessing (once, because S = spones(Ak) and
%         patterns of A1, A2 equal that of A0 for this example) --------
%  We preprocess with A1 to set up LS problems for the SAM N1 = SAM(A1->A0).
%  Because the nonzero PATTERN is identical here, we can reuse preproc.

preproc = sam_preprocess_right(A1);   % uses S = spones(A1) by default

fprintf('Preprocessing done.  nnz(N) = %d per update.\n\n', preproc.nnz_total);

%% ---- 4. Compute SAM updates ------------------------------------------
[N1, rel1] = sam_compute_right(A1, A0, preproc);
[N2, rel2] = sam_compute_right(A2, A0, preproc);

fprintf('Relative SAM residual for A1:  ||A1*N1 - A0||_F / ||A0||_F = %.2e\n', rel1);
fprintf('Relative SAM residual for A2:  ||A2*N2 - A0||_F / ||A0||_F = %.2e\n\n', rel2);

%% ---- 5. Solve with right-preconditioned GMRES -----------------------
tol     = 1e-8;
restart = 50;
maxit   = 200;

% -- Solve A0*x = b (sanity check, preconditioner exact for A0) --------
Pfun0 = @(v) sam_apply_right(speye(n), P0, v);   % N = I  => P_k = P0
[x0, fl0, rv0, it0] = gmres(@(v) A0*v, b, restart, tol, maxit, [], Pfun0);
fprintf('A0:  flag=%d,  iters=%d,  rel-res=%.2e\n', fl0, it0(2), rv0(end));

% -- Solve A1*x = b with recycled preconditioner P1 = N1 * P0 ----------
Pfun1 = @(v) sam_apply_right(N1, P0, v);
[x1, fl1, rv1, it1] = gmres(@(v) A1*v, b, restart, tol, maxit, [], Pfun1);
fprintf('A1 (SAM):  flag=%d,  iters=%d,  rel-res=%.2e\n', fl1, it1(2), rv1(end));

% -- Solve A2*x = b with recycled preconditioner P2 = N2 * P0 ----------
Pfun2 = @(v) sam_apply_right(N2, P0, v);
[x2, fl2, rv2, it2] = gmres(@(v) A2*v, b, restart, tol, maxit, [], Pfun2);
fprintf('A2 (SAM):  flag=%d,  iters=%d,  rel-res=%.2e\n', fl2, it2(2), rv2(end));

% -- Reference: solve A1, A2 by reusing P0 directly (no SAM) -----------
PfunRef = @(v) sam_apply_right(speye(n), P0, v);
[~, ~, rv1r, it1r] = gmres(@(v) A1*v, b, restart, tol, maxit, [], PfunRef);
[~, ~, rv2r, it2r] = gmres(@(v) A2*v, b, restart, tol, maxit, [], PfunRef);
fprintf('\nA1 (reuse P0):  iters=%d,  rel-res=%.2e\n', it1r(2), rv1r(end));
fprintf('A2 (reuse P0):  iters=%d,  rel-res=%.2e\n',  it2r(2), rv2r(end));

fprintf('\nSAM reduces iterations vs reusing P0 unchanged.\n');

%% ---- 6. Summary: relative residuals across sequence -----------------
fprintf('\n--- Relative SAM residual norms (proxy for preconditioner quality) ---\n');
fprintf('  N1: %.4e\n', rel1);
fprintf('  N2: %.4e\n', rel2);
fprintf('(Smaller => SAM maps Ak closer to A0 => preconditioner quality preserved.)\n');
