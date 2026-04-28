
clc;
clear;

%%

ffname = "moonContact";
files = dir('mats/'+ffname+'*.mat');
fileNames = {files.name};

% Extract trailing numbers to enforce numeric sort (prevents 1, 10, 2 ordering)
tokens = regexp(fileNames, ffname+'(\d+)\.mat', 'tokens', 'once');
nums = cellfun(@(x) str2double(x{1}), tokens);
[~, sortIdx] = sort(nums);
fileNames = fileNames(sortIdx);

sizeSeq = length(fileNames);
A = cell(sizeSeq, 1);
b = cell(sizeSeq, 1);
TV0 = cell(sizeSeq, 1);

for i = 1:sizeSeq
    % Load into struct to prevent namespace collision
    fname = "mats/"+fileNames{i};
    data = load(fname);
    
    % Replace 'mat_name' and 'rhs_name' with actual internal variable names
    A{i} = data.Amat;
    TV0{i} = data.TV0;
    b{i} = data.b;
end




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

[M,time] = computeRACP(A{2},TV0{1},true);
t_amg0    = time;
P0_amg    = M;
fprintf('done in %.3f s\n\n', t_amg0);

%% Solve the preconditioner

for i = 2:2

   Pfun_new = @(v) sam_apply_left(speye(size(A{i},1)), P0_amg, v);
   t0 = tic;
   [~, ~, relres(1,i), it_reuse(i,:),resvec(:,1)] = gmres_RIGHT(@(v) A{i}*v, b{i}, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
   t_solve_new(1) = toc(t0);
   it(i,1) = (it_reuse(i,1)-1)*gmres_restart + it_reuse(i,2);
   
   % Pfun_new = @(v) sam_apply_left(sam1{i}, P0_amg, v);
   % t0 = tic;
   % [~, ~, relres(2,i), it_sam1(i,:),resvec(:,2)] = gmres_RIGHT(@(v) A{i}*v, b{i}, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
   % t_solve_new(1) = toc(t0);
   % it(i,2) = (it_sam1(i,1)-1)*gmres_restart + it_sam1(i,2);
   % 
   % 
   % Pfun_new = @(v) sam_apply_left(sam2{i}, P0_amg, v);
   % t0 = tic;
   % [~, ~, relres(3,i), it_sam2(i,:),resvec(:,3)] = gmres_RIGHT(@(v) A{i}*v, b{i}, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
   % t_solve_new(1) = toc(t0);
   % it(i,3) = (it_sam2(i,1)-1)*gmres_restart + it_sam2(i,2);

end
relres = relres';

