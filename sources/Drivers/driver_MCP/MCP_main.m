clear
clc

DEBUG = false;
simple_flag = false;
treatBC = true;

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
rawNames = textscan(fileIN, '%s', 'Delimiter', '\n');

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
file_BIN     = normalizePath(fnames{8});

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
if ~exist('A','var')
   A = Amat;
end

fprintf('END INPUT\n\n');

%-----------------------------------------------------------------------------------------

% Check rhs input
if strcmpi(rhs_build,'unit_rhs')
   rhs = ones(size(A,1),1);
elseif strcmpi(rhs_build,'unit_sol')
   xx = ones(size(A,1),1);
   rhs = A*xx;
   clear xx;
elseif strcmpi(rhs_build,'rand_rhs')
   rhs = rand(size(A,1),1);
elseif strcmpi(rhs_build,'rhs_in')
   if ~exist('rhs','var')
      error('rhs input not found');
   end
else
   error('not expected type of rhs, try unit_rhs/unit_sol/rand_rhs/rhs_in');
end

% Split the matrix in 4 blocks

nn = size(A,1);
n11 = size(TV0,1);
n22 = nn - n11;
A11 = A(1:n11,1:n11);
A12 = A(1:n11,n11+1:end);
A21 = A(n11+1:end,1:n11);
A22 = A(n11+1:end,n11+1:end);

% Treat BC in A11
if treatBC
   % Identify target indices
   D = sum(spones(A11));
   ind_dir_dof = find(D==1);
   ind_col_rem = find(sum(spones(A12))==1);
   [ind_dir_lag,~,~] = find(A12(:,ind_col_rem));
   ind_dir = union(ind_dir_dof,ind_dir_lag);
   
   % Native Column Zeroing 
   A11(:, ind_dir) = 0;
   A21(:, ind_dir) = 0;
   
   % A11 Row Zeroing
   A11 = A11.';
   A11(:, ind_dir) = 0;
   A11 = A11.';
   
   % Coupling Matrix Resolution
   if sym_flag == 1
       % Symmetry Exploit: A21 columns are already zeroed. 
       % Transposing it perfectly zeros the corresponding rows for A12.
       A12 = A21.';
   else
       % Independent execution for non-symmetric systems
       A12 = A12.';
       A12(:, ind_dir) = 0;
       A12 = A12.';
   end
   
   % Diagonal Restoration
   fac = max(D);
   D_diag = zeros(n11, 1);
   D_diag(ind_dir, 1) = fac;
   A11 = A11 + spdiags(D_diag, 0, n11, n11);

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


A_aug = [AA BB; BB' CC];
rhs_aug = rhs_scaled;

fprintf('BEGIN: System solution\n');
switch lower(solv_method)

    case 'gmres_cj'

        % Solve the system by SQMR
        fprintf('BEGIN: System solution by GMRES_CJ\n');
        M = @(x) apply_MCP(MCP_prec,AA,BB,CC,x);
        [sol_aug,flag,relres,iter,resvec] = gmres_LEFT(A_aug,rhs_aug,restart,tol,itmax/restart,M);
        fprintf('END: System solution by GMRES_CJ\n\n');

    case 'sqmr'

        % Solve the system by SQMR
        fprintf('BEGIN: System solution by SQMR\n');
        Afun = @(x) A_aug*x;
        M = @(x) apply_MCP(MCP_prec,AA,BB,CC,x);
        IDfun = @(x) x;
        [sol_aug,flag,relres,iter,resvec] = SQMR(Afun,rhs_aug,tol,itmax,M,IDfun);
        fprintf('END: System solution by SQMR\n\n');

end
fprintf('END: System solution\n');

% Print results
if flag
   fprintf('Convergence not achieved\n');
end
fprintf('\n');
fprintf('\n');
fprintf('# of iterations:   %15d\n',iter);

fprintf('\n\nRESIDUAL NORMS:\n\n');
fprintf('SCALED RES   %15.6e    RHS %15.6e\n',norm(rhs_aug-A_aug*sol_aug),norm(rhs_aug));
fprintf('UNSCALED RES %15.6e    RHS %15.6e\n',norm(D_scal\(rhs_aug-A_aug*sol_aug)),...
                                              norm(D_scal\rhs_aug));
fprintf('PRECOND. RES %15.6e    RHS %15.6e\n',norm(M(rhs_aug-A_aug*sol_aug)),norm(M(rhs_aug)));
