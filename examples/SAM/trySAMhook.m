clc
clear 

%%
load mats/hypEl_nHook_cube.mat

%% Permute the matrix
perm = symrcm(A);

A_perm = A(perm,perm);
b_perm = b(perm);
TV0_perm = TV0(perm,:);

%% Take Away dirichlet nodes
n = 1864;
A_new = A_perm;
b_new = b_perm;
TV0_new = TV0_perm;

%% Scale and shift
D = diag(A_new);
A_scaled = sparse(diag((1./D))*A_new);
b_scaled = diag(1./D)*b_new;

alpha = [1e-5 1e-4 1e-3 5e-3 1e-2 5e-2 1e-1 5e-1].*mean(diag(A_scaled));

nn = size(A_scaled,1);
ID = speye(nn);

%% Compute Sam
clc
clear it
pre_mex = MEX_sam_preprocess_left(A_scaled);

restart = 400;
nn = size(A,1);

for i = 1:size(alpha,2)

   AA = A_scaled + alpha(i)*ID;

   % v = sort(eig(full(AA)));

   sam = MEX_sam_compute_left(AA,A_scaled,pre_mex);
   simpleNorm(i) = 2*norm(AA - A_scaled,'f')/(norm(AA,'f')+norm(A_scaled,'f'));
   samNorm(i) = 2*norm(sam*AA - A_scaled,'f')/(norm(sam*AA,'f')+norm(A_scaled,'f'));
   
   [M,time] = computeAMG(AA,TV0,false);

   S = @(x) sam_apply_left(speye(nn),M,x);

   [x, ~, ~, it_recomp] = gmres_RIGHT(@(v) AA*v, ones(size(AA,1),1), restart, 1e-6, 10, S, []);
   it(i) = (it_recomp(1)-1)*restart + it_recomp(2);

end


%%
options.type = 'ilutp';
options.droptol = 5e-3;
[L,U] = ilu(A,options);

%%
[x, ~, ~, it_LU] = gmres_RIGHT(@(v) A*v, b, restart, 1e-6, 100, L, U);
it1 = (it_LU(1)-1)*restart + it_LU(2);