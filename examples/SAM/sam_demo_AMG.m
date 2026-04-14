%SAM_DEMO_2D_AMG  LEFT SAM preconditioner updates with AMG base preconditioner.
%               2-D variable-coefficient convection-diffusion sequence.
%
%  PDE:  -div( mu(x,y,p)*grad u ) + c1(p)*u_x + c2(p)*u_y = f  on [0,1]^2
%
%  Strategies compared per system:
%    (1) Recompute AMG for every Ak                          [most accurate]
%    (2) Reuse P0 (AMG for A0) unchanged                     [cheapest/worst]
%    (3) Left SAM L1:  P_k = P0 * N_k,  sparsity = spones(Ak)
%    (4) Left SAM L2:  P_k = P0 * N_k,  sparsity = spones(Ak^2)  [denser]
%
%  Timing breakdown per strategy and system:
%    t_prec  — preconditioner setup (AMG Compute or SAM compute)
%    t_solve — GMRES solve time
%    t_total — t_prec + t_solve  (the operationally relevant number)
%
%  Note: for strategy (2) reuse, t_prec = 0 (no setup per system).
%        For strategies (3)/(4), t_prec is the SAM compute time only
%        (the cost of the initial AMG for A0 is amortised and shown
%        separately as t_amg0).
%
%  Files needed:
%    sam_preprocess_left.m  sam_compute_left.m  sam_apply_left.m

clear; clc; rng(42);
fprintf('================================================================\n');
fprintf(' SAM demo (LEFT, AMG) — 2-D variable-coeff conv-diffusion\n');
fprintf('================================================================\n\n');

%% ---- Parameters -------------------------------------------------------
m                   = 40;    % grid points per side  (n = m^2 = 1600)
Nseq                = 9;    % number of systems in the sequence
p0                  = 0.0;
dp                  = 0.15;

gmres_RIGHT_restart = 60;
gmres_RIGHT_tol     = 1e-8;
gmres_RIGHT_maxit   = 400;

%% ---- Build the sequence -----------------------------------------------
[A_seq, b_seq] = build_sequence(m, Nseq, p0, dp);
A0 = A_seq{1};
n  = size(A0, 1);
fprintf('Grid: %d x %d,  n = %d,  sequence length = %d\n\n', m, m, n, Nseq);

%% ---- AMG preconditioner P0 for A0 (computed once, time recorded) ------
fprintf('Computing AMG preconditioner for A0 ... ');
generalsolver                           = [];
generalsolver.simparams                 = [];
generalsolver.simparams.linSolverParams = [];
generalsolver.simparams.relTol          = 1e-8;
generalsolver.domains                   = ones(1,1);
generalsolver.nInterf                   = 0;
generalsolver.nDom                      = 1;

[AMG_prec,time] = computeAMG(A0,false);
P0_amg    = @(r) AMG_Vcycle(AMG_prec,A0,r);
t_amg0    = time; 

fprintf('done in %.3f s\n\n', t_amg0);

%% ---- Sparsity pattern preprocessing (once while pattern is fixed) -----
pat_L1     = spones(A0);
pat_L2     = spones(spones(A0)^2);

t0 = tic; preproc_L1 = sam_preprocess_left(A0, pat_L1); t_pre_L1 = toc(t0);
t0 = tic; preproc_L2 = sam_preprocess_left(A0, pat_L2); t_pre_L2 = toc(t0);

fprintf('SAM preprocessing:  L1 = %.3f s   L2 = %.3f s\n', t_pre_L1, t_pre_L2);
fprintf('nnz per row in N:   L1 = %.1f     L2 = %.1f\n\n', ...
        preproc_L1.nnz_total/n, preproc_L2.nnz_total/n);

%% ---- Storage ----------------------------------------------------------
% iteration counts
iter_new   = zeros(1, Nseq);
iter_reuse = zeros(1, Nseq);
iter_L1    = zeros(1, Nseq);
iter_L2    = zeros(1, Nseq);

% preconditioner setup times (per system)
t_prec_new   = zeros(1, Nseq);   % AMG recompute
t_prec_reuse = zeros(1, Nseq);   % 0 — nothing to do
t_prec_L1    = zeros(1, Nseq);   % SAM L1 compute
t_prec_L2    = zeros(1, Nseq);   % SAM L2 compute

% GMRES solve times (per system)
t_solve_new   = zeros(1, Nseq);
t_solve_reuse = zeros(1, Nseq);
t_solve_L1    = zeros(1, Nseq);
t_solve_L2    = zeros(1, Nseq);

% relative SAM residuals
relres_L1 = zeros(1, Nseq);
relres_L2 = zeros(1, Nseq);

%% ---- Main loop --------------------------------------------------------
for k = 1:Nseq
    Ak = A_seq{k};
    bk = b_seq{k};
    p  = p0 + (k-1)*dp;
    fprintf('--- System %2d / %2d  (p = %.3f) ---\n', k, Nseq, p);

    % Rebuild preprocessing only when the nonzero pattern of Ak changes
    if ~isequal(spones(Ak), spones(A0))
        t0 = tic; preproc_L1 = sam_preprocess_left(Ak, spones(Ak));           t_pre_L1 = toc(t0);
        t0 = tic; preproc_L2 = sam_preprocess_left(Ak, spones(spones(Ak)^2)); t_pre_L2 = toc(t0);
        fprintf('  [pattern changed] preprocessing: L1=%.3fs  L2=%.3fs\n', t_pre_L1, t_pre_L2);
    end

    % ---- 1. Recompute AMG for every Ak ----------------------------------
    [AMG_prec1,time] = computeAMG(Ak,false);
    Pk_amg = @(r) AMG_Vcycle(AMG_prec1,Ak,r);
    t_prec_new(k) = time;

    Pfun_new = @(v) sam_apply_left(speye(n), Pk_amg, v);
    t0 = tic;
    [~, fl, ~, it] = gmres_RIGHT(@(v) Ak*v, bk, ...
                                  gmres_RIGHT_restart, gmres_RIGHT_tol, ...
                                  gmres_RIGHT_maxit, Pfun_new, []);
    t_solve_new(k) = toc(t0);
    iter_new(k)    = extract_iters(it, fl, gmres_RIGHT_restart);

    % ---- 2. Reuse P0 unchanged (no setup cost) --------------------------
    t_prec_reuse(k) = 0;

    Pfun_reuse = @(v) sam_apply_left(speye(n), P0_amg, v);
    t0 = tic;
    [~, fl, ~, it] = gmres_RIGHT(@(v) Ak*v, bk, ...
                                  gmres_RIGHT_restart, gmres_RIGHT_tol, ...
                                  gmres_RIGHT_maxit, Pfun_reuse, []);
    t_solve_reuse(k) = toc(t0);
    iter_reuse(k)    = extract_iters(it, fl, gmres_RIGHT_restart);

    % ---- 3. Left SAM L1:  S = spones(Ak) --------------------------------
    t0 = tic;
    [Nl1, relres_L1(k)] = sam_compute_left(Ak, A0, preproc_L1);
    t_prec_L1(k) = toc(t0);

    Pfun_L1 = @(v) sam_apply_left(Nl1, P0_amg, v);
    t0 = tic;
    [~, fl, ~, it] = gmres_RIGHT(@(v) Ak*v, bk, ...
                                  gmres_RIGHT_restart, gmres_RIGHT_tol, ...
                                  gmres_RIGHT_maxit, Pfun_L1, []);
    t_solve_L1(k) = toc(t0);
    iter_L1(k)    = extract_iters(it, fl, gmres_RIGHT_restart);

    % ---- 4. Left SAM L2:  S = spones(Ak^2) ------------------------------
    t0 = tic;
    [Nl2, relres_L2(k)] = sam_compute_left(Ak, A0, preproc_L2);
    t_prec_L2(k) = toc(t0);

    Pfun_L2 = @(v) sam_apply_left(Nl2, P0_amg, v);
    t0 = tic;
    [~, fl, ~, it] = gmres_RIGHT(@(v) Ak*v, bk, ...
                                  gmres_RIGHT_restart, gmres_RIGHT_tol, ...
                                  gmres_RIGHT_maxit, Pfun_L2, []);
    t_solve_L2(k) = toc(t0);
    iter_L2(k)    = extract_iters(it, fl, gmres_RIGHT_restart);

    % Per-system progress line
    fprintf('  iters :  recomp=%3d | reuse=%3d | SAM-L1=%3d | SAM-L2=%3d\n', ...
            iter_new(k), iter_reuse(k), iter_L1(k), iter_L2(k));
    fprintf('  t_prec:  recomp=%5.3fs | reuse=%5.3fs | SAM-L1=%5.3fs | SAM-L2=%5.3fs\n', ...
            t_prec_new(k), t_prec_reuse(k), t_prec_L1(k), t_prec_L2(k));
    fprintf('  t_solv:  recomp=%5.3fs | reuse=%5.3fs | SAM-L1=%5.3fs | SAM-L2=%5.3fs\n', ...
            t_solve_new(k), t_solve_reuse(k), t_solve_L1(k), t_solve_L2(k));
    fprintf('  t_tot :  recomp=%5.3fs | reuse=%5.3fs | SAM-L1=%5.3fs | SAM-L2=%5.3fs\n', ...
            t_prec_new(k)+t_solve_new(k),   t_prec_reuse(k)+t_solve_reuse(k), ...
            t_prec_L1(k)+t_solve_L1(k),     t_prec_L2(k)+t_solve_L2(k));
    fprintf('  SAM residual:  L1 = %.2e   L2 = %.2e\n', relres_L1(k), relres_L2(k));
end

%% ---- Derived totals ---------------------------------------------------
% For SAM strategies, add the one-off AMG(A0) cost to the total so the
% comparison with recompute (which also uses AMG) is fair.
% Note: recompute pays AMG every step; SAM pays AMG once + SAM per step.
t_total_new   = sum(t_prec_new)   + sum(t_solve_new);
t_total_reuse = t_amg0            + sum(t_solve_reuse);  % AMG(A0) once
t_total_L1    = t_amg0 + t_pre_L1 + sum(t_prec_L1)   + sum(t_solve_L1);
t_total_L2    = t_amg0 + t_pre_L2 + sum(t_prec_L2)   + sum(t_solve_L2);

%% ---- Summary table ----------------------------------------------------
fprintf('\n');
fprintf('================================================================\n');
fprintf(' TIMING BREAKDOWN (seconds)\n');
fprintf('================================================================\n');
fprintf('  AMG(A0) setup (amortised over sequence):  %.3f s\n', t_amg0);
fprintf('  SAM preprocessing L1:                     %.3f s\n', t_pre_L1);
fprintf('  SAM preprocessing L2:                     %.3f s\n\n', t_pre_L2);

fprintf('%-8s  %8s  %8s  %8s  %8s\n', '', 'Recomp', 'Reuse', 'SAM-L1', 'SAM-L2');
fprintf('%-8s  %8s  %8s  %8s  %8s\n', '', '------', '------', '------', '------');
fprintf('%-8s  %8.3f  %8.3f  %8.3f  %8.3f\n', 't_prec', ...
        sum(t_prec_new), sum(t_prec_reuse), sum(t_prec_L1), sum(t_prec_L2));
fprintf('%-8s  %8.3f  %8.3f  %8.3f  %8.3f\n', 't_solve', ...
        sum(t_solve_new), sum(t_solve_reuse), sum(t_solve_L1), sum(t_solve_L2));
fprintf('%-8s  %8.3f  %8.3f  %8.3f  %8.3f\n', 't_total*', ...
        t_total_new, t_total_reuse, t_total_L1, t_total_L2);
fprintf('%-8s  %8d  %8d  %8d  %8d\n', 'iters', ...
        sum(iter_new), sum(iter_reuse), sum(iter_L1), sum(iter_L2));
fprintf('* t_total includes AMG(A0) setup + preprocessing for SAM strategies\n');
fprintf('================================================================\n\n');

fprintf(' k  |  p_k  | t_pr(R) t_sl(R) | t_pr(U) t_sl(U) | t_pr(L1) t_sl(L1) | t_pr(L2) t_sl(L2)\n');
fprintf('    |       |  itr(R)  tot(R)  |  itr(U)  tot(U) |  itr(L1)  tot(L1) |  itr(L2)  tot(L2)\n');
fprintf('----+-------+-----------------|-----------------|-------------------|-------------------\n');
for k = 1:Nseq
    fprintf('%3d | %.3f | %6.3f  %6.3f  | %6.3f  %6.3f  |  %6.3f   %6.3f  |  %6.3f   %6.3f\n', ...
        k, p0+(k-1)*dp, ...
        t_prec_new(k),   t_solve_new(k), ...
        t_prec_reuse(k), t_solve_reuse(k), ...
        t_prec_L1(k),    t_solve_L1(k), ...
        t_prec_L2(k),    t_solve_L2(k));
    fprintf('    |       |  %4d   %6.3f  |  %4d   %6.3f  |   %4d    %6.3f  |   %4d    %6.3f\n', ...
        iter_new(k),   t_prec_new(k)  +t_solve_new(k), ...
        iter_reuse(k), t_prec_reuse(k)+t_solve_reuse(k), ...
        iter_L1(k),    t_prec_L1(k)   +t_solve_L1(k), ...
        iter_L2(k),    t_prec_L2(k)   +t_solve_L2(k));
end
fprintf('================================================================\n');

%% ---- Plots ------------------------------------------------------------
ks    = 1:Nseq;
c_new = [0.1 0.1 0.1];   % near-black
c_reu = [0.8 0.1 0.1];   % red
c_L1  = [0.1 0.3 0.8];   % blue
c_L2  = [0.6 0.1 0.7];   % purple

figure('Name','Left SAM with AMG — 2D PDE','Color','w','Position',[100 60 900 780]);

% ---- subplot 1: GMRES iterations ----------------------------------------
subplot(3,1,1);
plot(ks, iter_new,   '-o','Color',c_new,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Recompute AMG'); hold on;
plot(ks, iter_reuse, '-s','Color',c_reu,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Reuse P_0');
plot(ks, iter_L1,    '-^','Color',c_L1, 'LineWidth',1.5,'MarkerSize',6,'DisplayName','Left SAM L1');
plot(ks, iter_L2,    '-d','Color',c_L2, 'LineWidth',1.5,'MarkerSize',6,'DisplayName','Left SAM L2');
ylabel('GMRES iterations');
title('GMRES iterations per system');
legend('Location','northwest'); grid on;
xlim([0.5, Nseq+0.5]);

% ---- subplot 2: timing breakdown (stacked: prec + solve) ----------------
subplot(3,1,2);
bar_data = [ t_prec_new.'    t_solve_new.' ; ...
             t_prec_reuse.'  t_solve_reuse.' ; ...
             t_prec_L1.'     t_solve_L1.' ; ...
             t_prec_L2.'     t_solve_L2.' ].';
% bar_data rows = [prec; solve], columns = [Recomp, Reuse, L1, L2] per system
% Reshape into a 2 x (4*Nseq) doesn't vectorise cleanly for grouped stacked;
% use a simpler grouped bar with total prec and solve side by side instead.
t_prec_mat = [t_prec_new.' t_prec_reuse.' t_prec_L1.' t_prec_L2.'];
t_solve_mat= [t_solve_new.' t_solve_reuse.' t_solve_L1.' t_solve_L2.'];
b = bar(ks, t_prec_mat + t_solve_mat, 'grouped');
b(1).FaceColor = c_new; b(2).FaceColor = c_reu;
b(3).FaceColor = c_L1;  b(4).FaceColor = c_L2;
hold on;
% Overlay hatched (lighter) portion for the solve time only
b2 = bar(ks, t_solve_mat, 'grouped');
b2(1).FaceColor = min(c_new+0.45,1); b2(1).FaceAlpha = 0.6;
b2(2).FaceColor = min(c_reu+0.45,1); b2(2).FaceAlpha = 0.6;
b2(3).FaceColor = min(c_L1+0.45,1);  b2(3).FaceAlpha = 0.6;
b2(4).FaceColor = min(c_L2+0.45,1);  b2(4).FaceAlpha = 0.6;
ylabel('Time (s)');
title('Per-system time  (dark = t_{prec},  light overlay = t_{solve})');
legend([b(1) b(2) b(3) b(4)], {'Recompute','Reuse','SAM-L1','SAM-L2'}, ...
       'Location','northwest');
grid on; xlim([0.5, Nseq+0.5]);

% ---- subplot 3: cumulative total time -----------------------------------
subplot(3,1,3);
cum_new   = cumsum(t_prec_new   + t_solve_new);
cum_reuse = cumsum(t_prec_reuse + t_solve_reuse) + t_amg0;  % include AMG(A0)
cum_L1    = cumsum(t_prec_L1    + t_solve_L1)    + t_amg0 + t_pre_L1;
cum_L2    = cumsum(t_prec_L2    + t_solve_L2)    + t_amg0 + t_pre_L2;

plot(ks, cum_new,   '-o','Color',c_new,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Recompute AMG'); hold on;
plot(ks, cum_reuse, '-s','Color',c_reu,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Reuse P_0');
plot(ks, cum_L1,    '-^','Color',c_L1, 'LineWidth',1.5,'MarkerSize',6,'DisplayName','Left SAM L1');
plot(ks, cum_L2,    '-d','Color',c_L2, 'LineWidth',1.5,'MarkerSize',6,'DisplayName','Left SAM L2');
xlabel('System index k');
ylabel('Cumulative time (s)');
title('Cumulative total time (prec setup + solve, incl. AMG(A_0) for SAM/Reuse)');
legend('Location','northwest'); grid on;
xlim([0.5, Nseq+0.5]);

% ---- separate figure: SAM residual quality ------------------------------
figure('Name','SAM residual quality','Color','w','Position',[100 60 700 340]);
semilogy(ks, relres_L1, '-^','Color',c_L1,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Left SAM L1'); hold on;
semilogy(ks, relres_L2, '-d','Color',c_L2,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Left SAM L2');
xlabel('System index k');
ylabel('||N_k A_k - A_0||_F / ||A_0||_F');
title('Relative SAM residual — quality of the approximate map');
legend('Location','northwest'); grid on;

%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function [A_seq, b_seq] = build_sequence(m, Nseq, p0, dp)
%BUILD_SEQUENCE  Assemble the PDE matrix sequence.
    h   = 1/(m+1);
    N   = m*m;
    idx = @(i,j) (j-1)*m + i;

    A_seq = cell(Nseq, 1);
    b_seq = cell(Nseq, 1);

    for seq = 1:Nseq
        p  = p0 + (seq-1)*dp;
        c1 = 5*p;
        c2 = 3*p;

        max_nnz = 5*N;
        ri = zeros(max_nnz,1);
        ci = zeros(max_nnz,1);
        vi = zeros(max_nnz,1);
        cnt = 0;

        for j = 1:m
            yj = j*h;
            for i = 1:m
                xi  = i*h;
                row = idx(i,j);

                muE = 1 + 0.8*p*sin(pi*(xi+h/2))*cos(pi*yj);
                muW = 1 + 0.8*p*sin(pi*(xi-h/2))*cos(pi*yj);
                muN = 1 + 0.8*p*sin(pi*xi)*cos(pi*(yj+h/2));
                muS = 1 + 0.8*p*sin(pi*xi)*cos(pi*(yj-h/2));

                dE = muE/h^2;  dW = muW/h^2;
                dN = muN/h^2;  dS = muS/h^2;

                cxP =  c1/(2*h);  cxM = -c1/(2*h);
                cyP =  c2/(2*h);  cyM = -c2/(2*h);

                % diagonal
                cnt=cnt+1; ri(cnt)=row; ci(cnt)=row;
                vi(cnt) = dE+dW+dN+dS;

                if i < m  % east
                    cnt=cnt+1; ri(cnt)=row; ci(cnt)=idx(i+1,j);
                    vi(cnt) = -dE+cxP;
                end
                if i > 1  % west
                    cnt=cnt+1; ri(cnt)=row; ci(cnt)=idx(i-1,j);
                    vi(cnt) = -dW+cxM;
                end
                if j < m  % north
                    cnt=cnt+1; ri(cnt)=row; ci(cnt)=idx(i,j+1);
                    vi(cnt) = -dN+cyP;
                end
                if j > 1  % south
                    cnt=cnt+1; ri(cnt)=row; ci(cnt)=idx(i,j-1);
                    vi(cnt) = -dS+cyM;
                end
            end
        end

        A_seq{seq} = sparse(ri(1:cnt), ci(1:cnt), vi(1:cnt), N, N);
        [XI, YJ]   = meshgrid((1:m)*h, (1:m)*h);
        F          = (2*pi^2*sin(pi*XI).*sin(pi*YJ)) * (1 + 0.1*p^2);
        b_seq{seq} = F(:);
    end
end

function total = extract_iters(it_vec, flag, restart)
    if flag == 0
        total = (it_vec(1)-1)*restart + it_vec(2);
    else
        total = 9999;   % did not converge within maxit
    end
end
