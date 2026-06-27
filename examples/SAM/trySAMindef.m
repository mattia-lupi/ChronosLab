clc
clear 

load mats/moonContact2.mat

A = Amat;
rhs = ones(size(A,1),1);
% TV0 = diag(1./diag(A));
% A = diag(1./diag(A))*A;
% p = symrcm(A);
% A = A(p,p);
% rhs = rhs(p);
% TV0 = TV0(p,:);

% d = full(diag(A));
% idx = (d == 1);
% A = A(idx==0,idx==0);
% rhs = rhs(idx==0);
% TV0 = TV0(idx==0,:);


%% Check Eigs

lambda = sort(eig(full(A)));


%% Shifted problem

alpha = [0 2e-3 3e-3 5e-3 7e-3 2.5e-1 3e-1 4e-1];%2
D = diag(diag(A));
lagrange = true;

gmres_restart = 100;
gmres_tol     = 1e-6;
gmres_maxit   = 20;

for i = 1:length(alpha)
   A_mod = A + alpha(i)*avg(D)*ones(size(A,1,1));

   lam = sort(eig(full(A_mod)));

   normLam(:,i) = abs(lam-lambda)./abs(lambda);

   if ~lagrange
      [M,time(1)] = computeAMG(A_mod,TV0,false);
   else
      [M,time(1)] = computeRACP(A_mod,TV0,false);
   end

   normm(i) = norm(A-A_mod,'f')/norm(A,'f');
   
   t_amg0    = time;
   P0_amg    = M;
   fprintf('done in %.3f s\n\n', t_amg0);

   t0 = tic;
   [~, ~, relres(i), it_recomp] = gmres_RIGHT(@(v) A_mod*v, rhs, gmres_restart, gmres_tol, gmres_maxit, M, []);
   t_solve_recomp(i) = toc(t0);
   it(i) = (it_recomp(1)-1)*gmres_restart + it_recomp(2);
end



%% Choose shift as beta = 3e-3

beta = 4e-1;

B = A + beta*D;

% % Define your parameters
% N = size(A,1);  % Size of the square matrix (N x N)
% n = 30;   % Number of upper and lower diagonals to fill with 1s
% 
% % 1. Define the diagonal offsets we want to target
% % For n=2, this creates the vector: [-2, -1, 0, 1, 2]
% diags_to_fill = -n:n;
% num_diags = length(diags_to_fill);
% 
% % 2. Create a columns-matrix of all ones for the values
% % spdiags expects a matrix where each column represents the values for a diagonal
% values = ones(N, num_diags);
% 
% % 3. Construct the sparse matrix
% S = spdiags(values, diags_to_fill, N, N);

tic
pre_mex = MEX_sam_preprocess_left(B);
time_prepr = toc;
tic

tic
sam = MEX_sam_compute_left(B,A,pre_mex);
time_compute = toc;
normmSam = 2 * norm(A-sam*B,'f')/(norm(A,'f')+norm(sam*B,'f'));

if ~lagrange
   [M,time(1)] = computeAMG(B,TV0,false);
else
   [M,time(1)] = computeRACP(B,TV0,false);
end

t0 = tic;
MM = @(v) sam_apply_left(sam, M, v);
[~, ~, relres_shift, it_recomp] = gmres_RIGHT(@(v) A*v, rhs, gmres_restart, gmres_tol, gmres_maxit, MM, []);
t_solve_shift = toc(t0);
it_shift = (it_recomp(1)-1)*gmres_restart + it_recomp(2);
