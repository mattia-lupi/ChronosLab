% Adaptive Smoother and Prolongation AMG main driver

clear;
%clc;
%close all;

rand('state',0);
RCM_flag = false;

fprintf('EXECUTION BEGIN\n\n');

% Read names of the input files
fileIN = fopen('eigdriver.fnames','r');
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_AMG     = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_SMOOTH  = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_TSPACE  = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_COARSEN = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_PROLONG = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_FILTER  = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_GENERAL = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_MATRIX  = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_TV0     = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_RHS     = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_BIN     = D{1};
fclose(fileIN);

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
[ascii_input,precond,method,neig,reslambda_check,itmax,tol] = read_eig(file_GENERAL);

% Set flag for symmetric matrices
param.symm = true;

% Input matrix, rhs and test space
if ascii_input

   % Input from ascii file

   % Read the system matrix
   A = read_csr(file_MATRIX);
   %%%%%%%%%%%%%
   % Treat Dir nodes
   lmax = eigs(A,1,'lm');
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

fprintf('BEGIN: Preconditioner computation\n');
t_init = tic;
switch lower(precond)
   case 'fsai'
      FSAI_prec = smoother(A, param.smoother);
      PREC_apply = @(x) FSAI_prec.right*(FSAI_prec.left*x);
   case 'amg'
      % Compute the AMG hierarchy
      AMG_prec = cpt_aspAMG(param,A,TV0);
      PREC_apply = @(r) AMG_Vcycle(AMG_prec,A,r);
end
T_prec = toc(t_init);
fprintf('END: Preconditioner computation\n\n');
%-----------------------------------------------------------------------------------------

%lmin = 5.1760e+03;
%MAT = @(x) PREC_apply(A*x);
%[V,D] = eigs(MAT,size(A,1),20,'sr','Display',1,'tolerance',1.e-8,'Failuretreatment','keep');
%Minv = PREC_apply(eye(size(A)));
%AA = Minv*(A-lmin*eye(size(A)));

tic;
fprintf('BEGIN: Eigenvalue computation\n');
switch lower(method)
    case 'srqcg'

      fprintf('BEGIN: Eigenvalue with SRQCG\n');
      V0 = rand(size(A,1),neig);
      ProdMat = @(x) A*x;
      ritz_freq = 1;
      [eigV,eigS,iter,resid] = NG_SRQCG(V0,ProdMat,PREC_apply,itmax,tol,ritz_freq);
      fprintf('END: Eigenvalue with SRQCG\n\n');
      tot_iter = iter*neig;

    case 'defl_srqcg'

      fprintf('BEGIN: Eigenvalue with Deflated SRQCG\n');
      V0 = rand(size(A,1),neig);
      ProdMat = @(x) A*x;
      ritz_freq = 1;

      % Load the space for deflation
      if isfile('DEFL_SPACE.mat')
         load DEFL_SPACE;
      else
         fprintf("\nThis function requires the deflation space to work\n");
         return;
      end

      [eigV,eigS,iter,resid] = DEFL_SRQCG(V0,ProdMat,PREC_apply,itmax,tol,ritz_freq,YY);
      fprintf('END: Eigenvalue with Deflated SRQCG\n\n');
      tot_iter = iter*neig;

    case 'lanczos'

      global count;
      global tot_iter;
      count = 0;
      tot_iter = 0;
      fprintf('BEGIN: Eigenvalue with LANCZOS\n');
      % Define application of the inverse
      Ainv = @(x) SolvePCG(A,PREC_apply,x,1000,1.e-8);
      [eigV,eigS] = eigs(Ainv,size(A,1),neig,'lm','Display',1);
      eigS = diag(eigS);
      fprintf('END: Eigenvalue with LANCZOS\n\n');
      iter = count;

    case 'lobpcg'

      largest_flag = false;
      restartControl = 7;
      prodA = @(x) A*x;
      prodB = [];
      fprintf('BEGIN: Eigenvalue with LOBPCG\n');
      if isfile('DEFL_SPACE.mat')
         load DEFL_SPACE;
      else
         YY = [];
      end
         
      X0 = rand(size(A,1),neig);

      % Define application of the inverse
      [iter,eigS,eigV,resnorm_vec,lambda_vec,ierr] =...
             lobpcg(prodA,prodB,PREC_apply,YY,X0,largest_flag,...
                    reslambda_check,itmax,tol,restartControl);
      fprintf('END: Eigenvalue with LOBPCG\n\n');
      tot_iter = iter*neig;

end
fprintf('END: Eigenvalue computation\n');
T_iter = toc;

%-----------------------------------------------------------------------------------------

fprintf('\n');
fprintf('------------------------------------------------------------------------------------------');
fprintf('\n');

fprintf('\n');
fprintf('\n');
fprintf('# of iterations:            %15d\n',iter);
fprintf('# of matvet applications:   %15d\n',tot_iter);
fprintf('\n');
fprintf('Preconditioner time: %10.2f\n',T_prec);
fprintf('EigenSolver time:    %10.2f\n',T_iter);

fprintf('\n');
for i = 1:neig
    fprintf('Eigenvalue (%d): %12.6f\n', i, eigS(i,1));
end

if(isfile('DEFL_SPACE.mat'))
   YY = [YY eigV];
   % Save the eigenvectors to use them as a deflation space
   save("DEFL_SPACE","YY");
end

return
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
fprintf(ofile,'Tspace time:       %10.2f\n',T_tspa);
fprintf(ofile,'Coarsening time:   %10.2f\n',T_coar);
fprintf(ofile,'Comp. Relax. time: %10.2f\n',T_CR);
fprintf(ofile,'Prolongation time: %10.2f\n',T_prol);
fprintf(ofile,'Iteration time:    %10.2f\n',T_iter);

fprintf('\n');
for i = 1:neig
    fprintf(ofile,'Eigenvalue (%d): %12.6f\n', i, D(i,1));
end

% Print report on EMIN
if strcmp(upper(param.prolong.prol_emin),'MEX_EMIN');
   prt_EMIN_Repo(AMG_prec);
end

fclose(ofile);
