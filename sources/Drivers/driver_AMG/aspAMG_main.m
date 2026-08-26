% Adaptive Smoother and Prolongation AMG main driver

clear;
%clc;
%close all;

% Modern random state initialization (rand('state',0) is deprecated)
rng(0);

RCM_flag = true;
RCM_flag = false;

fprintf('EXECUTION BEGIN\n\n');

% Read all lines at once into a cell array
fileIN = fopen('aspAMG.fnames','r');
rawNames = textscan(fileIN, '%s', 'Delimiter', '\n');
fclose(fileIN);

fnames = rawNames{1};

% Helper to normalize slashes for the current OS (Windows '\', Unix '/')
normalizePath = @(p) strrep(strrep(strtrim(p), '/', filesep), '\', filesep);

% Assign and normalize each path
file_AMG     = normalizePath(fnames{1});
file_SMOOTH  = normalizePath(fnames{2});
file_TSPACE  = normalizePath(fnames{3});
file_COARSEN = normalizePath(fnames{4});
file_PROLONG = normalizePath(fnames{5});
file_FILTER  = normalizePath(fnames{6});
file_GENERAL = normalizePath(fnames{7});
file_MATRIX  = normalizePath(fnames{8});
file_TV0     = normalizePath(fnames{9});
file_RHS     = normalizePath(fnames{10});
file_BIN     = normalizePath(fnames{11});

% Read parameters for the AMG hierarchy
param.amg = read_amg(file_AMG);

% Read parameters for the smoother
param.smoother = read_smoother(file_SMOOTH);

% Read parameters for the testspace
param.tspace = read_tspace(file_TSPACE);

% Read parameters for the smoother
param.coarsen = read_coarsen(file_COARSEN);

% Read parameters for the prolongation
param.prolong = read_prolong(file_PROLONG);

% Read parameters for the filtering
param.filter = read_filter(file_FILTER);

% Read general parameters
[ascii_input,rhs_build,sym_flag,solv_method,itmax,tol,restart] =...
        read_general(file_GENERAL);

% Set the symmetry flag
param.symm = sym_flag;

% Input matrix, rhs and test space
if ascii_input

   % Input from ascii file

   % Read the system matrix
   A = read_csr(file_MATRIX);
   %%%%%%%%%%%%%
   % Treat Dir nodes
   lmax = eigs(A,1,'lm','Tolerance',0.01,'MaxIterations',10);
   A(A==1) = lmax/10;
   %%%%%%%%%%%%%

   % Read the initial test space
   TV0 = dlmread(file_TV0,'',1,0);

   if strcmp(lower(rhs_build),'rhs_in')
      rhs = load(file_RHS);
   end

else

   % Input from .mat file
   load(file_BIN);
   %%%%%%%%%%%%%
   % Treat Dir nodes
   lmax = eigs(A,1,'lm','FailureTreatment','keep','Display',0,'Tolerance',0.001,'MaxIterations',3);
   A(A==1) = lmax/10;
   %%%%%%%%%%%%%

   if ~exist('rhs') && exist('b')
      rhs = b;
      TV0 = ones(size(A,1),1);
   end

   if ~exist('rhs') && strcmp(lower(rhs_build),'rhs_in')
      err_msg = 'Missing the expected right-hand side';
   end

end

fprintf('END INPUT\n\n');

%-----------------------------------------------------------------------------------------

% Sort the matrix if equired
if RCM_flag
   perm = symrcm(A);
   A = A(perm,perm);
   TV0 = TV0(perm,:);
   rhs = rhs(perm,:);
end

%-----------------------------------------------------------------------------------------

% Compute right-hand side
switch lower(rhs_build)
   case 'unit_sol'
      rhs = A*ones(size(A,1),1);
   case 'unit_rhs'
      rhs = ones(size(A,1),1);
   case 'rand_sol'
      rand_sol = rand(size(A,1),1);
      rhs = A*rand_sol;
   case 'rand_rhs'
      rhs = rand(size(A,1),1).*10.^(6*(rand(size(A,1),1)-0.5));
   case 'rhs_in'

   otherwise
      err_msg = ['Wrong value for rhs_build in ' filename];
      error(err_msg);
end

%-----------------------------------------------------------------------------------------

global DEBINFO;
% PARTE GENERALE
DEBINFO.flag = true;
DEBINFO.flag = false;
% PARTE PER LA PROLONGATION
DEBINFO.prol = [];
% STAMPARE SI/NO
DEBINFO.prol.prt_flag = false;
% UNITA DI STAMPA
DEBINFO.prol.ofile = 0;
% STAMPA NUMERO DI ITERAZIONI NEL CALCOLO PROL
DEBINFO.prol.it_print = false;
% STAMPA LISTA VICINI NEL CALCOLO PROL
DEBINFO.prol.neigh_print = false;

% PARTE PER IL COARSENING
DEBINFO.coarsen = [];
% STAMPARE SI/NO
DEBINFO.coarsen.draw_dist = false;

%-----------------------------------------------------------------------------------------

% Set timings
amg_times
T_smoo = 0;
T_LowRank = 0;
T_tspa = 0;
T_coar = 0;
T_CR = 0;
T_prol = 0;
T_FilCLev = 0;
T_FilProl = 0;
T_setup = 0;

% Compute the AMG hierarchy
fprintf('BEGIN: Preconditioner computation\n');
time_start = tic;
AMG_prec = cpt_aspAMG(param,A,TV0);
T_setup = toc(time_start);
fprintf('END: Preconditioner computation\n\n');

%-----------------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%sol = AMG_Vcycle(AMG_prec,A,rhs);
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

tic;
fprintf('BEGIN: System solution\n');
switch lower(solv_method)
    case 'stat_amg'

      % Solve the system by AMG stationary iteration
      fprintf('BEGIN: System solution by stationary AMG\n');
      [sol,flag,relres,iter,resvec] = stat_AMG(A,rhs,tol,itmax,AMG_prec);
      fprintf('END: System solution by stationary AMG\n\n');

    case 'pcg'

       % Solve the system by PCG
       fprintf('BEGIN: System solution by PCG\n');
       Mfun = @(r) AMG_Vcycle(AMG_prec,A,r);
       %[sol,flag,relres,iter,resvec] = pcg(A,rhs,tol,itmax,Mfun);
       prodA = @(x) A*x;
       [sol,iter,resvec,flag] = PCG_CJ(Mfun,prodA,rhs,zeros(size(rhs)),itmax,tol);
       relres = resvec(end)/resvec(1);
       fprintf('END: System solution by PCG\n\n');

    case 'bicgstab'

      % Solve the system by BiCGstab
      fprintf('BEGIN: System solution by BiCGstab\n');
      Mfun = @(r) AMG_Vcycle(AMG_prec,A,r);
      [sol,flag,relres,iter,resvec] = bicgstab(A,rhs,tol,itmax,Mfun);
      fprintf('END: System solution by BiCGstab\n\n');

    case 'gmres'

      % Solve the system by GMRES
      fprintf('BEGIN: System solution by GMRES\n');
      Mfun = @(r) AMG_Vcycle(AMG_prec,A,r);
      [sol,flag,relres,iter,resvec] = gmres(A,rhs,restart,tol,itmax,Mfun);
      fprintf('END: System solution by GMRES\n\n');

    case 'gmres_cj'

      % Solve the system by GMRES
      fprintf('BEGIN: System solution by GMRES_CJ\n');
      Mfun = @(r) AMG_Vcycle(AMG_prec,A,r);
      [sol,flag,relres,iter,resvec] = gmres_RIGHT(A,rhs,restart,tol,itmax,Mfun);
      fprintf('END: System solution by GMRES_CJ\n\n');

    case 'sqmr'

      % Solve the system by SQMR
      fprintf('BEGIN: System solution by SQMR\n');
      Afun = @(x) A*x;
      Mfun = @(r) AMG_Vcycle(AMG_prec,A,r);
      IDfun = @(x) x;
      [sol,flag,relres,iter,resvec] = SQMR(Afun,rhs,tol,itmax,Mfun,IDfun);
      fprintf('END: System solution by GMRES_CJ\n\n');


end
fprintf('END: System solution\n');
T_iter = toc;

%-----------------------------------------------------------------------------------------

fprintf('\n');
fprintf('------------------------------------------------------------------------------------------');
fprintf('\n');

% Get AMG hierarchy information
AMG_info = get_AMG_info(AMG_prec,A);

% Print AMG hierarchy information
print_AMG_info(AMG_info);

fprintf('\n');
fprintf('Grid Complexity:     %10.4f\n',AMG_info(1).grid_comp);
fprintf('Operator Complexity: %10.4f\n',AMG_info(1).oper_comp);

% Print results
if flag
   fprintf('Convergence not achieved\n');
end
fprintf('\n');
fprintf('\n');
fprintf('# of iterations:   %15d\n',iter);
fprintf('Relative residual: %15.6e\n',relres);
fprintf('Real residual:     %15.6e\n',norm(rhs-A*sol)/norm(rhs));

fprintf('\n');
fprintf('Smoother time:     %10.2f\n',T_smoo);
fprintf('Low Rank time:     %10.2f\n',T_LowRank);
fprintf('Tspace time:       %10.2f\n',T_tspa);
fprintf('Coarsening time:   %10.2f\n',T_coar);
fprintf('Comp. Relax. time: %10.2f\n',T_CR);
fprintf('Prolongation time: %10.2f\n',T_prol);
fprintf('Filter oper. time: %10.2f\n',T_FilCLev);
fprintf('Filter prol. time: %10.2f\n',T_FilProl);
fprintf('Total Set-Up:      %10.2f\n',T_setup);
fprintf('Iteration time:    %10.2f\n',T_iter);
fprintf('Total time:        %10.2f\n',T_setup+T_iter);

% Print results in a file

ofile = fopen('Out_AMG','w');

fprintf(ofile,'\n');
fprintf(ofile,'\n');

% Print AMG hierarchy information
print_AMG_info(AMG_info,ofile);

fprintf(ofile,'\n');
fprintf(ofile,'Grid Complexity:     %10.4f\n',AMG_info(1).grid_comp);
fprintf(ofile,'Operator Complexity: %10.4f\n',AMG_info(1).oper_comp);

% Print results
if flag
   fprintf(ofile,'Convergence not achieved\n');
end
fprintf(ofile,'\n');
fprintf(ofile,'\n');
fprintf(ofile,'# of iterations:   %15d\n',iter);
fprintf(ofile,'Relative residual: %15.6e\n',relres);
fprintf(ofile,'Real residual:     %15.6e\n',norm(rhs-A*sol)/norm(rhs));

fprintf(ofile,'\n');
fprintf(ofile,'Smoother time:     %10.2f\n',T_smoo);
fprintf(ofile,'Low Rank time:     %10.2f\n',T_LowRank);
fprintf(ofile,'Tspace time:       %10.2f\n',T_tspa);
fprintf(ofile,'Coarsening time:   %10.2f\n',T_coar);
fprintf(ofile,'Comp. Relax. time: %10.2f\n',T_CR);
fprintf(ofile,'Prolongation time: %10.2f\n',T_prol);
fprintf(ofile,'Total Set-Up:      %10.2f\n',T_setup);
fprintf(ofile,'Iteration time:    %10.2f\n',T_iter);
fprintf(ofile,'Total time:        %10.2f\n',T_setup+T_iter);

% Print report on EMIN
if strcmp(upper(param.prolong.prol_emin),'MEX_EMIN');
   prt_EMIN_Repo(AMG_prec);
end

fclose(ofile);
