
clc;
clear;

sizeSeq = 4;

% Load Matrices and rhs
for i = 1:sizeSeq
   fname = "mats/flowNonConforming_" + num2str(i-1) + ".mat";
   prob = load(fname);

   A{i} = prob.Amat;
   b{i} = prob.b;
   TV0{i} = ones(1435,1);
end

clear prob fname i


%% See how much they differ in norm

normm = zeros(sizeSeq,1);

for i = 1:sizeSeq
   normm(i) = 2 * norm(A{1}-A{i},'f')/(norm(A{1},'f')+norm(A{i},'f'));
end

%% Compute Sam

pre_mex1 = MEX_sam_preprocess_left(A{1});
pre_mex2 = MEX_sam_preprocess_left(A{1},A{1}^2);
for i = 1:sizeSeq
   [sam1{i},normmSam1(i)] = MEX_sam_compute_left(A{i},A{1},pre_mex1);
   normm1(i) = 2 * norm(A{1}-sam1{i}*A{i},'f')/(norm(A{1},'f')+norm(sam1{i}*A{i},'f'));

   [sam2{i},normmSam2(i)] = MEX_sam_compute_left(A{i},A{1},pre_mex2);
   normm2(i) = 2 * norm(A{1}-sam2{i}*A{i},'f')/(norm(A{1},'f')+norm(sam2{i}*A{i},'f'));
end

%% Compute the Preconditioner

gmres_restart = 100;
gmres_tol     = 1e-6;
gmres_maxit   = 10;

[M,time] = computeRACP(A{1},TV0{1},false);
t_amg0    = time;
P0_amg    = M;
fprintf('done in %.3f s\n\n', t_amg0);

%% Solve the preconditioner

for i = 1:sizeSeq


   Pfun_new = @(v) sam_apply_left(speye(size(A{i},1)), P0_amg, v);
   t0 = tic;
   [~, ~, ~, it_reuse(i,:)] = gmres_RIGHT(@(v) A{i}*v, b{i}, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
   t_solve_new(1) = toc(t0);
   it(i,1) = (it_reuse(i,1)-1)*gmres_restart + it_reuse(i,2);
   
   Pfun_new = @(v) sam_apply_left(sam1{i}, P0_amg, v);
   t0 = tic;
   [~, ~, ~, it_sam1(i,:)] = gmres_RIGHT(@(v) A{i}*v, b{i}, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
   t_solve_new(1) = toc(t0);
   it(i,2) = (it_sam1(i,1)-1)*gmres_restart + it_sam1(i,2);
   
   
   Pfun_new = @(v) sam_apply_left(sam2{i}, P0_amg, v);
   t0 = tic;
   [~, ~, ~, it_sam2(i,:)] = gmres_RIGHT(@(v) A{i}*v, b{i}, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
   t_solve_new(1) = toc(t0);
   it(i,3) = (it_sam2(i,1)-1)*gmres_restart + it_sam2(i,2);
end

