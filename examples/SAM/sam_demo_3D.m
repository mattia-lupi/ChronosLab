%SAM_DEMO_3D_AMG  LEFT SAM preconditioner updates with AMG base preconditioner.
%               3-D variable-coefficient convection-diffusion sequence.
%
%  PDE:  -div( mu(x,y,z,p)*grad u ) + c1(p)*u_x + c2(p)*u_y + c3(p)*u_z = f  on [0,1]^3
%
%  Strategies compared per system:
%    (1) Recompute AMG for every Ak                          [most accurate]
%    (2) Reuse P0 (AMG for A0) unchanged                     [cheapest/worst]
%    (3) Left SAM L1 (MEX):     P_k = P0 * N_k,  sparsity = spones(Ak)
%
%  Files needed:
%    sam_preprocess_left.m      sam_compute_left.m          sam_apply_left.m
%    MEX_sam_preprocess_left.m  MEX_sam_compute_left.m

clear; clc; rng(42);
fprintf('================================================================\n');
fprintf(' SAM demo (LEFT, AMG) — 3-D variable-coeff conv-diffusion\n');
fprintf('================================================================\n\n');

%% ---- Parameters -------------------------------------------------------
m                   = 15;    % grid points per side (n = m^3 = 3375)
Nseq                = 9;     % number of systems in the sequence
p0                  = 0.0;
dp                  = 0.15;

gmres_RIGHT_restart = 100;
gmres_RIGHT_tol     = 1e-8;
gmres_RIGHT_maxit   = 5;

%% ---- Build the sequence -----------------------------------------------
[A_seq, b_seq] = build_sequence(m, Nseq, p0, dp);
A0 = A_seq{1};
n  = size(A0, 1);
fprintf('Grid: %d x %d x %d,  n = %d,  sequence length = %d\n\n', m, m, m, n, Nseq);

%% ---- MEX Verification & Benchmark -------------------------------------
fprintf('--- Verifying MEX implementations against MATLAB ---\n');
Ak_test = A_seq{2};
pat_test = spones(Ak_test);

% 1. Preprocess Check
pre_mat = sam_preprocess_left(Ak_test, pat_test);
pre_mex = MEX_sam_preprocess_left(Ak_test, pat_test);
assert(pre_mat.nnz_total == pre_mex.nnz_total, 'nnz_total mismatch in preprocess!');

% 2. Compute Check
[N_mat, ~] = sam_compute_left(Ak_test, A0, pre_mat);
[N_mex, ~] = MEX_sam_compute_left(Ak_test, A0, pre_mex);

err_N = norm(N_mat - N_mex, 'fro') / norm(N_mat, 'fro');
fprintf('|| N_mat - N_mex ||_F / || N_mat ||_F = %e\n', err_N);
if err_N > 1e-10
    warning('MEX compute output differs significantly from MATLAB formulation.');
else
    fprintf('MEX correctness verification PASSED.\n');
end

% 3. Benchmark (Averaged over 5 runs)
runs = 5;
fprintf('\n--- Benchmarking MEX vs MATLAB (averaged over %d runs) ---\n', runs);
tic; for i=1:runs, sam_preprocess_left(Ak_test, pat_test); end; t_pre_m = toc/runs;
tic; for i=1:runs, MEX_sam_preprocess_left(Ak_test, pat_test); end; t_pre_x = toc/runs;
tic; for i=1:runs, sam_compute_left(Ak_test, A0, pre_mat); end; t_com_m = toc/runs;
tic; for i=1:runs, MEX_sam_compute_left(Ak_test, A0, pre_mex); end; t_com_x = toc/runs;

fprintf('Preprocess: MATLAB = %.4fs | MEX = %.4fs | Speedup = %.2fx\n', t_pre_m, t_pre_x, t_pre_m/t_pre_x);
fprintf('Compute:    MATLAB = %.4fs | MEX = %.4fs | Speedup = %.2fx\n\n', t_com_m, t_com_x, t_com_m/t_com_x);

%% ---- AMG preconditioner P0 for A0 (computed once, time recorded) ------

[AMG_prec,time] = computeAMG(A0,false);
t_amg0    = time;
P0_amg    = @(x) AMG_Vcycle(AMG_prec,A0,b_seq{0});
fprintf('done in %.3f s\n\n', t_amg0);

%% ---- Sparsity pattern preprocessing (once while pattern is fixed) -----
pat_L1 = spones(A0);

t0 = tic; preproc_L1_mex = MEX_sam_preprocess_left(A0, pat_L1); t_pre_L1_mex = toc(t0);

fprintf('SAM preprocessing L1 (MEX): %.3f s\n', t_pre_L1_mex);
fprintf('nnz per row in N: L1 = %.1f\n\n', preproc_L1_mex.nnz_total/n);

%% ---- Storage ----------------------------------------------------------
iter_new       = zeros(1, Nseq);
iter_reuse     = zeros(1, Nseq);
iter_L1_mex    = zeros(1, Nseq);

t_prec_new     = zeros(1, Nseq);
t_prec_reuse   = zeros(1, Nseq);
t_prec_L1_mex  = zeros(1, Nseq);

t_solve_new    = zeros(1, Nseq);
t_solve_reuse  = zeros(1, Nseq);
t_solve_L1_mex = zeros(1, Nseq);

relres_L1_mex  = zeros(1, Nseq);

%% ---- Main loop --------------------------------------------------------
for k = 1:Nseq
    Ak = A_seq{k};
    bk = b_seq{k};
    p  = p0 + (k-1)*dp;
    fprintf('--- System %2d / %2d  (p = %.3f) ---\n', k, Nseq, p);

    if ~isequal(spones(Ak), spones(A0))
        pat_Ak = spones(Ak);
        t0 = tic; preproc_L1_mex = MEX_sam_preprocess_left(Ak, pat_Ak); t_pre_L1_mex = toc(t0);
        fprintf('  [pattern changed] preproc: MEX-L1=%.3fs\n', t_pre_L1_mex);
    end

    % 1. Recompute AMG
    [AMG_prec,time] = computeAMG(Ak,false);
    Pk_amg = @(x) AMG_Vcycle(AMG_prec,Ak,bk);
    t_prec_new(k) = time;

    Pfun_new = @(v) sam_apply_left(speye(n), Pk_amg, v);
    t0 = tic;
    [~, fl, ~, it] = gmres_RIGHT(@(v) Ak*v, bk, gmres_RIGHT_restart, gmres_RIGHT_tol, gmres_RIGHT_maxit, Pfun_new, []);
    t_solve_new(k) = toc(t0);
    iter_new(k)    = extract_iters(it, fl, gmres_RIGHT_restart);

    % 2. Reuse P0
    t_prec_reuse(k) = 0;
    Pfun_reuse = @(v) sam_apply_left(speye(n), P0_amg, v);
    t0 = tic;
    [~, fl, ~, it] = gmres_RIGHT(@(v) Ak*v, bk, gmres_RIGHT_restart, gmres_RIGHT_tol, gmres_RIGHT_maxit, Pfun_reuse, []);
    t_solve_reuse(k) = toc(t0);
    iter_reuse(k)    = extract_iters(it, fl, gmres_RIGHT_restart);

    % 3. Left SAM L1 (MEX)
    if k == 2 || k == 5 || k == 8
      t0 = tic;
      [Nl1_mex, relres_L1_mex(k)] = MEX_sam_compute_left(Ak, A0, preproc_L1_mex);
      t_prec_L1_mex(k) = toc(t0);
    else
       if k == 1, Nl1_mex = speye(n); end
       t_prec_L1_mex(k) = 0;
    end

    Pfun_L1_mex = @(v) sam_apply_left(Nl1_mex, P0_amg, v);
    t0 = tic;
    [~, fl, ~, it] = gmres_RIGHT(@(v) Ak*v, bk, gmres_RIGHT_restart, gmres_RIGHT_tol, gmres_RIGHT_maxit, Pfun_L1_mex, []);
    t_solve_L1_mex(k) = toc(t0);
    iter_L1_mex(k)    = extract_iters(it, fl, gmres_RIGHT_restart);

    % System output
    fprintf('  iters : recomp=%3d | reuse=%3d | L1(MEX)=%3d\n', ...
            iter_new(k), iter_reuse(k), iter_L1_mex(k));
    fprintf('  t_prec: recomp=%5.3fs | reuse=%5.3fs | L1(MEX)=%5.3fs\n', ...
            t_prec_new(k), t_prec_reuse(k), t_prec_L1_mex(k));
    fprintf('  t_solv: recomp=%5.3fs | reuse=%5.3fs | L1(MEX)=%5.3fs\n', ...
            t_solve_new(k), t_solve_reuse(k), t_solve_L1_mex(k));
end

%% ---- Derived totals ---------------------------------------------------
t_total_new     = sum(t_prec_new)     + sum(t_solve_new);
t_total_reuse   = t_amg0              + sum(t_solve_reuse); 
t_total_L1_mex  = t_amg0 + t_pre_L1_mex + sum(t_prec_L1_mex) + sum(t_solve_L1_mex);

%% ---- Summary table ----------------------------------------------------
fprintf('\n================================================================\n');
fprintf(' TIMING BREAKDOWN (seconds)\n');
fprintf('================================================================\n');
fprintf('%-8s  %8s  %8s  %9s\n', '', 'Recomp', 'Reuse', 'SAM-L1(X)');
fprintf('%-8s  %8s  %8s  %9s\n', '', '------', '------', '---------');
fprintf('%-8s  %8.3f  %8.3f  %9.3f\n', 't_prec', ...
        sum(t_prec_new), sum(t_prec_reuse), sum(t_prec_L1_mex));
fprintf('%-8s  %8.3f  %8.3f  %9.3f\n', 't_solve', ...
        sum(t_solve_new), sum(t_solve_reuse), sum(t_solve_L1_mex));
fprintf('%-8s  %8.3f  %8.3f  %9.3f\n', 't_total*', ...
        t_total_new, t_total_reuse, t_total_L1_mex);
fprintf('%-8s  %8d  %8d  %9d\n', 'iters', ...
        sum(iter_new), sum(iter_reuse), sum(iter_L1_mex));
fprintf('================================================================\n\n');

%% ---- Plots ------------------------------------------------------------
ks    = 1:Nseq;
c_new  = [0.1 0.1 0.1];
c_reu  = [0.8 0.1 0.1];
c_L1_x = [0.1 0.7 0.8];

figure('Name','Left SAM with AMG — 3D PDE','Color','w','Position',[100 60 1000 780]);

subplot(3,1,1);
plot(ks, iter_new,   '-o','Color',c_new,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Recompute AMG'); hold on;
plot(ks, iter_reuse, '-s','Color',c_reu,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Reuse P_0');
plot(ks, iter_L1_mex,'-v','Color',c_L1_x,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Left SAM L1 (MEX)');
ylabel('GMRES iterations'); title('GMRES iterations per system');
legend('Location','northwest'); grid on; xlim([0.5, Nseq+0.5]);

subplot(3,1,2);
t_prec_mat = [t_prec_new.' t_prec_reuse.' t_prec_L1_mex.'];
t_solve_mat= [t_solve_new.' t_solve_reuse.' t_solve_L1_mex.'];
b = bar(ks, t_prec_mat + t_solve_mat, 'grouped');
b(1).FaceColor = c_new; b(2).FaceColor = c_reu; b(3).FaceColor = c_L1_x;
hold on;
b2 = bar(ks, t_solve_mat, 'grouped');
b2(1).FaceColor = min(c_new+0.45,1); b2(1).FaceAlpha = 0.6;
b2(2).FaceColor = min(c_reu+0.45,1); b2(2).FaceAlpha = 0.6;
b2(3).FaceColor = min(c_L1_x+0.45,1);b2(3).FaceAlpha = 0.6;
ylabel('Time (s)'); title('Per-system time (dark = t_{prec}, light overlay = t_{solve})');
legend([b(1) b(2) b(3)], {'Recomp','Reuse','SAM-L1(X)'}, 'Location','northwest');
grid on; xlim([0.5, Nseq+0.5]);

subplot(3,1,3);
cum_new    = cumsum(t_prec_new    + t_solve_new);
cum_reuse  = cumsum(t_prec_reuse  + t_solve_reuse)  + t_amg0;
cum_L1_mex = cumsum(t_prec_L1_mex + t_solve_L1_mex) + t_amg0 + t_pre_L1_mex;

plot(ks, cum_new,   '-o','Color',c_new,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Recompute AMG'); hold on;
plot(ks, cum_reuse, '-s','Color',c_reu,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Reuse P_0');
plot(ks, cum_L1_mex,'-v','Color',c_L1_x,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Left SAM L1 (MEX)');
xlabel('System index k'); ylabel('Cumulative time (s)');
title('Cumulative total time (prec setup + solve, incl. AMG(A_0))');
legend('Location','northwest'); grid on; xlim([0.5, Nseq+0.5]);

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function [A_seq, b_seq] = build_sequence(m, Nseq, p0, dp)
    h   = 1/(m+1);
    N   = m^3;
    idx = @(i,j,k) (k-1)*m^2 + (j-1)*m + i;

    A_seq = cell(Nseq, 1);
    b_seq = cell(Nseq, 1);

    for seq = 1:Nseq
        p  = p0 + (seq-1)*dp;
        c1 = 5*p;
        c2 = 3*p;
        c3 = 2*p;

        max_nnz = 7*N;
        ri = zeros(max_nnz,1);
        ci = zeros(max_nnz,1);
        vi = zeros(max_nnz,1);
        cnt = 0;

        for k = 1:m
            zk = k*h;
            for j = 1:m
                yj = j*h;
                for i = 1:m
                    xi  = i*h;
                    row = idx(i,j,k);

                    muE = 1 + 0.8*p*sin(pi*(xi+h/2))*cos(pi*yj)*sin(pi*zk);
                    muW = 1 + 0.8*p*sin(pi*(xi-h/2))*cos(pi*yj)*sin(pi*zk);
                    muN = 1 + 0.8*p*sin(pi*xi)*cos(pi*(yj+h/2))*sin(pi*zk);
                    muS = 1 + 0.8*p*sin(pi*xi)*cos(pi*(yj-h/2))*sin(pi*zk);
                    muF = 1 + 0.8*p*sin(pi*xi)*cos(pi*yj)*sin(pi*(zk+h/2));
                    muB = 1 + 0.8*p*sin(pi*xi)*cos(pi*yj)*sin(pi*(zk-h/2));

                    dE = muE/h^2;  dW = muW/h^2;
                    dN = muN/h^2;  dS = muS/h^2;
                    dF = muF/h^2;  dB = muB/h^2;

                    cxP =  c1/(2*h);  cxM = -c1/(2*h);
                    cyP =  c2/(2*h);  cyM = -c2/(2*h);
                    czP =  c3/(2*h);  czM = -c3/(2*h);

                    cnt=cnt+1; ri(cnt)=row; ci(cnt)=row;
                    vi(cnt) = dE+dW+dN+dS+dF+dB;

                    if i < m
                        cnt=cnt+1; ri(cnt)=row; ci(cnt)=idx(i+1,j,k);
                        vi(cnt) = -dE+cxP;
                    end
                    if i > 1
                        cnt=cnt+1; ri(cnt)=row; ci(cnt)=idx(i-1,j,k);
                        vi(cnt) = -dW+cxM;
                    end
                    if j < m
                        cnt=cnt+1; ri(cnt)=row; ci(cnt)=idx(i,j+1,k);
                        vi(cnt) = -dN+cyP;
                    end
                    if j > 1
                        cnt=cnt+1; ri(cnt)=row; ci(cnt)=idx(i,j-1,k);
                        vi(cnt) = -dS+cyM;
                    end
                    if k < m
                        cnt=cnt+1; ri(cnt)=row; ci(cnt)=idx(i,j,k+1);
                        vi(cnt) = -dF+czP;
                    end
                    if k > 1
                        cnt=cnt+1; ri(cnt)=row; ci(cnt)=idx(i,j,k-1);
                        vi(cnt) = -dB+czM;
                    end
                end
            end
        end

        A_seq{seq} = sparse(ri(1:cnt), ci(1:cnt), vi(1:cnt), N, N);
        [XI, YJ, ZK] = meshgrid((1:m)*h, (1:m)*h, (1:m)*h);
        F          = (3*pi^2*sin(pi*XI).*sin(pi*YJ).*sin(pi*ZK)) * (1 + 0.1*p^2);
        b_seq{seq} = F(:);
    end
end

function total = extract_iters(it_vec, flag, restart)
    if flag == 0
        total = (it_vec(1)-1)*restart + it_vec(2);
    else
        total = 9999;
    end
end
