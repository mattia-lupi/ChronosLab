clear
clc

simple_flag = false;
DEBUG = false;

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

% Read names of the input files
fileIN = fopen('MCP.fnames','r');
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_AMG     = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_SMOOTH  = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_TSPACE  = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_COARSEN = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_PROLONG = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_FILTER  = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_GENERAL = D{1};
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
[ascii_input,rhs_build,sym_flag,solv_method,itmax_ruiz,tol_ruiz,itmax,tol,restart] =...
        read_general_MCP(file_GENERAL);

% Load the system
load(file_BIN);

fprintf('END INPUT\n\n');

%-----------------------------------------------------------------------------------------


% Split the matrix in 4 blocks

nn = size(A,1);
n11 = size(TV0,1);
n22 = nn - n11;
A11 = A(1:n11,1:n11);
A12 = A(1:n11,n11+1:end);
A21 = A(n11+1:end,1:n11);
A22 = A(n11+1:end,n11+1:end);

% Treat BC in A11
tic;
if true
   D = sum(spones(A11));
   ind_dir_dof = find(D==1);
   ind_col_rem = find(sum(spones(A12))==1);
   [ind_dir_lag,~,~] = find(A12(:,ind_col_rem));
   ind_dir = union(ind_dir_dof,ind_dir_lag);
   A11(:,ind_dir) = 0;
   A11 = A11';
   A11(:,ind_dir) = 0;
   A21(:,ind_dir) = 0;
   A12 = A21';
   fac = max(D);
   D = zeros(n11,1);
   D(ind_dir,1) = fac;
   A11 = A11 + diag(sparse(D));

   % Remove BC columns from A12 (and A21) and A22
   ind_col_retain = setdiff(1:size(A22,1),ind_col_rem)';
   A12 = A12(:,ind_col_retain);
   A21 = A12';
   A22 = A22(ind_col_retain,ind_col_retain);
   rhs = rhs([(1:n11)'; ind_col_retain]);
end

% Replace A
A_orig = [A11 A12; A21 A22];
nn = size(A_orig,1);
n22 = size(A22,1);
rhs_orig = rhs;
time = toc;
fprintf('Tempo BC %f\n',time);

% Apply Ruiz to the entire matriz
tic;
[A_r,D_scal] = ruiz_symmetric(A_orig,itmax_ruiz,tol_ruiz);
A11_scaled = A_r(1:n11,1:n11);
A12_scaled = A_r(1:n11,n11+1:end);
A21_scaled = A_r(n11+1:end,1:n11);
A22_scaled = A_r(n11+1:end,n11+1:end);
D_11 = full(diag(D_scal(1:n11,1:n11)));
D_22 = full(diag(D_scal(n11+1:end,n11+1:end)));
TV0 = diag(sparse(D_11))\TV0;
D_scal = diag(sparse([D_11; D_22]));
time = toc;
fprintf('Tempo scal %f\n',time);

A_scaled = [A11_scaled A12_scaled; A21_scaled A22_scaled];
rhs_scaled = D_scal*rhs_orig;

if DEBUG
   print_SpMat('SCALED_K',A11_scaled);
   print_SpMat('SCALED_B',A12_scaled);
   print_SpMat('SCALED_C',A22_scaled);
end

%-----------------------------------------------------------------------------------------

% Compute global augmentation
lmax_glo = eigs(A11_scaled,1,'lm','Display',1,'Tolerance',1.e-5,...
               'MaxIterations',20,'FailureTreatment','keep');

ADD_scaled = A12_scaled*A21_scaled;
if DEBUG
   print_SpMat('ADD_SCALED',ADD_scaled);
end
D = full(diag(ADD_scaled));
D = D(D>0);
gmean_ADD_s = geomean(D);
fprintf('lmax_glo:            %e\n',lmax_glo);
fprintf('geometric mean BBT:  %e\n',gmean_ADD_s);
fprintf('Augmentation factor: %e\n',lmax_glo/gmean_ADD_s);

% Compute augmentation block
gamma = 1.0;
D22_mat = gamma*(lmax_glo/gmean_ADD_s)*speye(n22);

% Compute new saddle-point system peconditioner
param.symm = sym_flag;
AA = A11_scaled + A12_scaled*D22_mat*A21_scaled; AA = 0.5*(AA+AA');
BB = A12_scaled;
CC = A22_scaled; CC = 0.5*(CC+CC');
if DEBUG
   tmp = [AA BB; BB' CC];
   print_SpMat('AUGMENTED_K',AA);
   print_SpMat('AUGMENTED_B',BB);
   print_SpMat('AUGMENTED_C',CC);
   clear tmp;
end

% Compute preconditioner
MCP_prec = cpt_MCP(param,AA,TV0,BB,CC);

%-----------------------------------------------------------------------------------------
%fprintf('INIZIO PCG\n');
%PREC = @(x) AMG_Vcycle(MCP_prec.AMG_prec,AA,x);
%ProdA = @(x) AA*x;
%rhs11 = ones(size(AA,1),1);
%[x11,iter11,resvec11,flag11] = PCG_CJ(PREC,ProdA,rhs11,zeros(size(AA,1),1),itmax,tol);
%res11 = rhs11-AA*x11;
%D11_scal = D_scal(1:n11,1:n11);
%fprintf('\n\nRESIDUAL NORMS PCG:\n\n');
%fprintf('SCALED RES   %15.6e RHS %15.6e\n',norm(res11),norm(rhs11));
%fprintf('UNSCALED RES %15.6e RHS %15.6e\n',norm(D11_scal\res11),norm(D11_scal\rhs11));
%-----------------------------------------------------------------------------------------

nrest = 100;
M = @(x) apply_MCP(MCP_prec,AA,BB,CC,x);
A_aug = [AA BB; BB' CC];
rhs_aug = rhs_scaled;
[sol_aug,flag,relres,iter,resvec] = gmres_LEFT(A_aug,rhs_aug,nrest,1.e-8,4,M);

fprintf('\n\nRESIDUAL NORMS:\n\n');
fprintf('SCALED RES   %15.6e RHS %15.6e\n',norm(rhs_aug-A_aug*sol_aug),norm(rhs_aug));
fprintf('UNSCALED RES %15.6e RHS %15.6e\n',norm(D_scal\(rhs_aug-A_aug*sol_aug)),...
                                           norm(D_scal\rhs_aug));
fprintf('PRECOND. RES %15.6e RHS %15.6e\n',norm(M(rhs_aug-A_aug*sol_aug)),norm(M(rhs_aug)));
fprintf('\n');
