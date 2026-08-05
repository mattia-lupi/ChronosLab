function [iter, t_setup, t_solve] = solveSystem(A, b, M_all, params)
    % SOLVESYSTEM Executes linear solver evaluation for all 3 strategies
    sizeSeq = length(A);
    iter    = params.iter;
    t_setup = params.t_setup;
    t_solve = params.t_solve;

    P_reuse       = [];
    P_sam         = [];
    idx_prec_base = 1; 
    sm            = [];

    for i = 1:sizeSeq
       % --- Strategy 1: Recompute Preconditioner Every Step ---
       M_recomp = M_all{i};
       
       t0 = tic;
       [~, ~, ~, it_recomp] = gmres_RIGHT(@(v) A{i}*v, b{i}, params.gmres_restart, params.gmres_tol, params.gmres_maxit, M_recomp, []);
       t_solve(i, 1) = toc(t0);
       iter(i, 1)    = (it_recomp(1)-1)*params.gmres_restart + it_recomp(2);
       
       % --- Strategy 2: Reuse/Recycle Preconditioner ---
       if i == 1 || ismember(i, params.q)
          P_reuse = M_recomp;
          t_setup(i, 2) = t_setup(i, 1);
       else
          t_setup(i, 2) = 0;
       end
       
       Pfun_reuse = @(v) sam_apply_left(speye(size(A{i}, 1)), P_reuse, v);
       
       t0 = tic;
       [~, ~, ~, it_reuse] = gmres_RIGHT(@(v) A{i}*v, b{i}, params.gmres_restart, params.gmres_tol, params.gmres_maxit, Pfun_reuse, []);
       t_solve(i, 2) = toc(t0);
       iter(i, 2)    = (it_reuse(1)-1)*params.gmres_restart + it_reuse(2);
       
       % --- Strategy 3: SAM Adaptive ---
       if ~isempty(params.tS)
          if i == 1 || ismember(i, params.tP)
             P_sam = M_recomp;
             t_prec_sam = t_setup(i, 1);
             t_setup(i, 3) = t_setup(i, 3) + t_prec_sam;
          
             idx_prec_base = i;  
             sm = speye(size(A{i}, 1)); 
          end
          
          if ismember(i, params.tS)
             if i == idx_prec_base
                sm = speye(size(A{i}, 1));
             else
                t0_sam = tic;
                [sm, ~] = MEX_sam_adaptive_left(A{i}, A{idx_prec_base}, 8, params.nstep, params.step_size, params.epss);
                t_setup(i, 3) = t_setup(i, 3) + toc(t0_sam);
             end
          end
          
          Pfun_sam = @(v) sam_apply_left(sm, P_sam, v);
          
          t0 = tic;
          [~, ~, ~, it_sam] = gmres_RIGHT(@(v) A{i}*v, b{i}, params.gmres_restart, params.gmres_tol, params.gmres_maxit, Pfun_sam, []);
          t_solve(i, 3) = toc(t0);
          iter(i, 3)    = (it_sam(1)-1)*params.gmres_restart + it_sam(2);
       end
    end
end