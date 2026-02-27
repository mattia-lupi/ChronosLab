function [AMG_hrc] = cpt_aspAMG(param,A,TV0)
%-----------------------------------------------------------------------------------------
%
% Function representing the entry point for the recursive computation of the AMG hierarchy
%
% Input
%
% param:     structure containing all the parameters for AMG construction
% A:         system matrix
% TV0:       initial testspace which is eventually padded to the desired size
%
% Output
%
% AMG_hrc:  AMG hierachy, a linked list of all the AMG components at every level
%
%-----------------------------------------------------------------------------------------

global AMG_Total_Mem;
global AMG_Peak_Mem;
AMG_Total_Mem = 0;
AMG_Peak_Mem = 0;

% Init the level counter
level = 0;

% Padd the test space with random vectors if TV0 is not large enough
fprintf('BEGIN: Initializing initial test space\n');
n = size(A,1);
ntv = param.tspace.ntv;
ntv0 = size(TV0,2);
TV = zeros(n,ntv);
TV(:,1:ntv0) = TV0;
TV(:,ntv0+1:ntv) = rand(n,ntv-ntv0);

% Orthonormalize the test space
[TV,~] = qr(TV,0);
fprintf('END: Initializing initial test space\n');

AMG_hrc = REC_cpt_aspAMG(level,param,A,TV);

GB_AMG = whos('AMG_hrc').bytes/(1024^3);
fprintf('*******************************************************************\n');
fprintf('******* Total AMG memory occupation (GB): %10.4f *******\n',GB_AMG);
fprintf('******* Peak memory occupation (GB):      %10.4f *******\n',AMG_Peak_Mem);
fprintf('*******************************************************************\n');

return

%-----------------------------------------------------------------------------------------

function [AMG_hrc] = REC_cpt_aspAMG(level,param,A,TV0)

global AMG_Total_Mem;
global AMG_Peak_Mem;
global DEBINFO;
% Remove some DEBUG file
if DEBINFO.flag
   if level == 1
      ! rm Logfile.*
      ! rm IBDEBUG
      ! rm ROWDEBUG
      ! rm LOG_CORRECT
   end
end

amg_times

% Increase the level counter
level = level + 1;
AMG_hrc.level = level;
fprintf('\n*****   LEVEL: %5d   ***** \n\n',level);

% Get the current size
n = size(A,1);
% Compute operator max eigenvalue
lmax_A = eigs(A,1,'lm','FailureTreatment','keep','Display',0,'Tolerance',0.001,'MaxIterations',10);
AMG_hrc.lmax_A = lmax_A;

% Print operator info
fprintf('Operator # of rows:               %10d\n',n);
fprintf('Operator # of non-zeroes:         %10d\n',nnz(A));
fprintf('Operator # of non-zeroes per row: %10.2f\n',nnz(A)/n);
fprintf('Operator max lambda:              %15.6e\n\n',lmax_A);

% Compute the preconditioner for this level
if n <= param.amg.maxCoarseSZ || level >= param.amg.nLevMax

   %***** This level is the last one *****

   % Compute the Factorization and store it
   if param.symm
      [L,p,S] = chol(A,'lower');
      if p ~= 0
         fprintf('ERROR IN FACTORIZING LAST LEVEL\n');
      end
      AMG_hrc.L = L;
      AMG_hrc.S = S;
   else
      fprintf('Using UNSYMMETRIC factorization\n');
      [L,U,P,Q] = lu(A);
      % P*A*Q = L*U ---> A*x = b ---> P'*LU*Q'*x = b ---> LU*Q'*x = P*b ---> x = Q*(U\(L\(P*b))
      AMG_hrc.L = L;
      AMG_hrc.U = U;
      AMG_hrc.P = P;
      AMG_hrc.Q = Q;
   end
   AMG_hrc.symm = param.symm;
   AMG_hrc.A = A;
   AMG_hrc.next = 0;

else

   %***** Compute this level using adaptive AMG *****

   % Compute the smoother
   time_start = tic;
   fprintf('BEGIN: Computing the smoother\n');
   smootherOp = smoother(A, param.symm, param.smoother);
   fprintf('END: Computing the smoother\n\n');
   T_smoo = T_smoo + toc(time_start);

   %--------------------------------------------------------------------------------------

   % Compute the test space
   tic;
   fprintf('BEGIN: Computing the test space\n');
   if param.symm
      [TV, lambda, res] = tspace(TV0, A, smootherOp, param.tspace);
      fprintf('   i          lambda     res_lam     res_vec\n');
      nl = numel(lambda);
      fprintf('%4d %15.6e %11.2e %11.2e\n',[(1:nl)' lambda(end:-1:1) res(end:-1:1,:)]');
   else
      fprintf('Unsymmetric eigensolver not available yet\n');
      fprintf('The initial test space will be used\n');
      TV = TV0;
   end
   fprintf('END: Computing the test space\n\n');
   T_tspa = T_tspa + toc;

   %--------------------------------------------------------------------------------------

   % Check if aggressive coarsening is required
   use_aggressive = level <= param.coarsen.nl_agg;

   % Compute coarsening
   tic;
   fprintf('BEGIN: Computing the coarse nodes\n');
   if ~use_aggressive
      [fcnode,clist,flist,S,S_patt] = coarsen(param.coarsen,A,smootherOp,TV);
      num_f = numel(flist);
      num_c = numel(clist);
      fclist = 0;
      fflist = 0;
   else
      [fcnode,clist,fclist,fflist,S,S_patt] = agg_coarsen(param.coarsen,A,TV);
      num_f = numel(fclist) + numel(fflist);
      num_c = numel(clist);
   end
   fprintf('END: Computing the coarse nodes\n\n');
   T_coar = T_coar + toc;
   fprintf('Number of FINE nodes:   %10d\n',num_f);
   fprintf('Number of COARSE nodes: %10d\n',num_c);
   fprintf('Coarse node percentage: %10.2f\n\n',100*num_c/numel(fcnode));
   if num_f + num_c ~= numel(fcnode)
      err_msg = 'Number of F/C node is not consistent with level size';
      error(err_msg);
   end

   %--------------------------------------------------------------------------------------

   % Compute prolongation
   tprol = tic;
   fprintf('BEGIN: Computing the prolongation\n');
   [P,clist,fcnode,emin_info] = prolong(use_aggressive,level,param.symm,param,...
                                        clist,fclist,fflist,fcnode,S,S_patt,smootherOp,A,TV);
   fprintf('END: Computing the prolongation\n\n');
   prol_tot_time = toc(tprol);
   fprintf('Prolongation time %10.3f [s]\n',prol_tot_time);
   T_prol = T_prol + prol_tot_time;

   % Store EMIN detailed information
   if ( param.symm && strcmpi(param.prolong.prol_emin,'EMIN') )
      emin_info.time_prolTot  = prol_tot_time;
   end
   AMG_hrc.emin_info = emin_info;

   % Create next level test space (by simple injection of coarse nodes)
   nc = numel(clist);
   TVnext(1:nc,:) = TV(clist,:);

   %--------------------------------------------------------------------------------------

   % Filter prolongation
   tic;
   type = param.filter.filt_type;
   np   = param.filter.np;
   wgt  = param.filter.filt_wgt;
   tol  = param.filter.filt_tol;
   if wgt < 100
      fprintf('Filtering prolongation\n');
      switch upper(type)
         case 'MTLB_FILT'
            % MATLAB filtering
            Pf = Filter_Prol(wgt,tol,TVnext,P);
         case 'MEX_FILT'
            % MEX filtering
            Pf = MEX_FiltProl(np,wgt,tol,TVnext,P);
         otherwise
            err_msg = strcat('This kind of filtering is not supported:  ',type);
            error(err_msg);
      end
   else
      Pf = P;
   end
   T_FilProl = T_FilProl + toc;
   [nn,~] = size(P);
   fprintf('Filtered prolongation non-zeroes per row: %10.2f\n',(nnz(Pf)-nc)/(nn-nc));
   fprintf('Pf density over P: %f\n', nnz(Pf) / nnz(P));
   if size(P,1) < 10000
      fprintf('Conditioning of P: %f\n',condest(Pf'*Pf));
   end

   %--------------------------------------------------------------------------------------

   % Create next level operator
   Anext = Pf'*(A*Pf);
   % Enforce symmetry if a symmetric accelerator is required
   if param.symm
%      Anext = 0.5*(Anext+Anext');
   end

   % Get memory occupation before clean-up
   GB_used = GetWhosMemory;
   AMG_Peak_Mem = max(AMG_Peak_Mem,AMG_Total_Mem+GB_used);

   % Store the components of the current level
   if level > 1
      AMG_hrc.A = A; clear A;
   end
   AMG_hrc.P = P; clear P;
   AMG_hrc.Pf = Pf; clear Pf;
   AMG_hrc.S = S; clear S;
   AMG_hrc.S_patt = S_patt; clear S_patt;
   AMG_hrc.TV = TV; clear TV;
   AMG_hrc.fcnode = fcnode;
   AMG_hrc.omega =smootherOp.omega;
   if strcmp(lower(param.smoother.method),'blk_j')
      AMG_hrc.Snnz = 0;
   elseif strcmp(lower(param.smoother.method),'bafsai')
      AMG_hrc.Snnz = 0;
   elseif strcmp(lower(param.smoother.method),'ddsw')
      AMG_hrc.Snnz = smootherOp.nnz;
   else
      AMG_hrc.Snnz = nnz(smootherOp.right) + nnz(smootherOp.left);
   end
   AMG_hrc.nupre = param.smoother.nupre;
   AMG_hrc.nupost = param.smoother.nupost;
   % Simple smoother
   AMG_hrc.Minv1 = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.left*x));
   AMG_hrc.Minv2 = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.left*x));

   %--------------------------------------------------------------------------------------

   % Filter the next level operator
   tic;
   type = param.filter.filt_type;
   tau = param.filter.filt_tau;
   np = param.filter.np;
   patt_min_flag = param.filter.min_patt;
   fprintf('Avg nnzr Anext before filtering: %f\n',nnz(Anext)/size(Anext,1));
   if tau > 0
      fprintf('Filtering next level operator\n');
      switch upper(type)
         case 'MTLB_FILT'
            % MATLAB filtering
            Anext = Filter_CoarseLev(patt_min_flag,tau,A,Anext,fcnode,Pf,TVnext);
         case 'MEX_FILT'
            % MEX filtering
            Anext = MEX_FiltCLEV(np,patt_min_flag,tau,A,Anext,fcnode,Pf,TVnext);
         otherwise
            err_msg = strcat('This kind of filtering is not supported:  ',type);
            error(err_msg);
      end
   end
   fprintf('Avg nnzr Anext after filtering: %f\n',nnz(Anext)/size(Anext,1));
   T_FilCLev = T_FilCLev + toc;

   %--------------------------------------------------------------------------------------
   GB_used = GetWhosMemory;
   AMG_Total_Mem = AMG_Total_Mem + GB_used;
   AMG_Peak_Mem = max(AMG_Peak_Mem,AMG_Total_Mem);
   fprintf('\n');
   fprintf('************************************************************\n');
   fprintf('******* Level %2d memory occupation (GB):  %10.4f *******\n',level,GB_used);
   fprintf('******* Total memory occupation (GB):     %10.4f *******\n',AMG_Total_Mem);
   fprintf('******* Peak memory occupation (GB):      %10.4f *******\n',AMG_Peak_Mem);
   fprintf('************************************************************\n');
   fprintf('\n');

   % Compute next level of the hierarchy
   AMG_hrc.next = REC_cpt_aspAMG(level,param,Anext,TVnext);

end

return
