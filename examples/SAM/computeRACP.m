function [M,T_setup] = computeRACP(A,TV0,sym_flag)

DEBUG = false;
simple_flag = false;
treatBC = true;
itmax_ruiz = 0;
tol_ruiz = 1;

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


fileIN = fopen('aspAMG.fnames','r');
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_AMG     = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_SMOOTH  = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_TSPACE  = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_COARSEN = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_PROLONG = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_FILTER  = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_GENERAL = D{1};
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
[ascii_input,rhs_build,~,solv_method,itmax,tol,restart] =...
        read_general(file_GENERAL);

% Set the symmetry flag
param.symm = sym_flag;

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
end

% Replace A
A_orig = [A11 A12; A21 A22];
nn = size(A_orig,1);
n22 = size(A22,1);
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

if DEBUG
   print_SpMat('SCALED_K',A11_scaled);
   print_SpMat('SCALED_B',A12_scaled);
   print_SpMat('SCALED_C',A22_scaled);
end

%-----------------------------------------------------------------------------------------

GLOBAL = false;
gamma = 1.0;

if GLOBAL

   % Compute global augmentation

   lmax_glo = eigs(A11_scaled,1,'lm','Display',1,'Tolerance',1.e-5,...
                  'MaxIterations',20,'FailureTreatment','keep');
   fprintf('Eigmax_glo: %15.6e\n',lmax_glo);

   ADD_scaled = A12_scaled*A21_scaled;
   if DEBUG
      print_SpMat('ADD_SCALED',ADD_scaled);
   end
   D = full(diag(ADD_scaled));
   D = D(D>0);
   gmean_ADD_s = geomean(D);
   fprintf('gmean_ADD_s: %15.6e\n',gmean_ADD_s);

   ind = find(full(diag(ADD_scaled)) > 0);
   X = A11_scaled(ind,ind);
   fprintf('Gloabl Augmentation factor: %e\n',lmax_glo/gmean_ADD_s);

   inv_D22 = gamma*(lmax_glo/gmean_ADD_s)*speye(n22);

   param.symm = sym_flag;
   A11_aug = A11_scaled + A12_scaled*inv_D22*A21_scaled; A11_aug = 0.5*(A11_aug+A11_aug');

else

   % Compute local augmentation
   
   % Compute augmentation
   AUG_FILE = fopen('AUG_FACT','w');
   AA_list = {};
   BB_list = {};
   aug = zeros(size(A22_scaled,1),1);
   D_11 = full(diag(A11_scaled));
   mean_diag_A = mean(D_11);
   D_22 = full(diag(A22_scaled));
   A21_scaled_T = A21_scaled';
   for icol = 1:n22
      v12 = A12_scaled(:,icol);
      v21 = A21_scaled_T(:,icol);
      [ii_12,~,bb_12] = find(v12);
      [ii_21,~,bb_21] = find(v21);
      if (numel(ii_12)+numel(ii_21) > 0)
         BB = bb_12*bb_21';
         if simple_flag 
            m_a= max(D_11(ii_12));
            m_b = max(diag(BB));
         else
            AA = A11_scaled(ii_12,ii_21);
            AA_list{icol} = AA;
            BB_list{icol} = BB;
            m_a = max(eig(full(AA)));
            m_b = max(eig(full(BB)));
         end
         m_a_sav = m_a;
         if m_a == 0
            m_a = mean_diag_A;
         end
         alpha = m_a / m_b;
         fprintf(AUG_FILE,'%15.6e %15.6e %15.6e %15.6e\n',m_a,m_b,alpha,m_a_sav);
         aug(icol) = 1 / alpha;
      end
   end
   fclose(AUG_FILE);

   aug_mat = diag(sparse(aug));
   A22_aug = A22_scaled - gamma*aug_mat;

   % Compute augmented 11 block
   param.symm = sym_flag;
   inv_D22 = -inv(diag(diag(A22_aug)));
   ADD = A12_scaled*inv_D22*A21_scaled; 

   % Symmetrise strongly if the matrix is supposed to be symmetric
   if sym_flag == 1
      ADD = 0.5*(ADD+ADD');
   end
   A11_aug = A11_scaled+ADD;

end

% Compute reverse augmented peconditioner
time_start = tic;
AMG_prec = cpt_aspAMG(param,A11_aug,TV0,true);
T_setup = toc(time_start);

M = @(x) apply_RevAug(AMG_prec,A11_aug,A12_scaled,inv_D22,x);

end
