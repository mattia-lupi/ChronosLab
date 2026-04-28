
clc;
clear;

files = dir('mats/richards_*.mat');
fileNames = {files.name};

% Extract trailing numbers to enforce numeric sort (prevents 1, 10, 2 ordering)
tokens = regexp(fileNames, 'richards_(\d+)\.mat', 'tokens', 'once');
nums = cellfun(@(x) str2double(x{1}), tokens);
[~, sortIdx] = sort(nums);
fileNames = fileNames(sortIdx);

sizeSeq = length(fileNames);
A = cell(sizeSeq, 1);
b = cell(sizeSeq, 1);

for i = 1:sizeSeq
    % Load into struct to prevent namespace collision
    fname = "mats/"+fileNames{i};
    data = load(fname);
    
    % Replace 'mat_name' and 'rhs_name' with actual internal variable names
    A{i} = data.Amat;
    b{i} = data.b;
end




%% See how much they differ in norm

normm = zeros(sizeSeq,1);

for i = 1:sizeSeq
   normm(i) = 2 * norm(A{1}-A{i},'f')/(norm(A{1},'f')+norm(A{i},'f'));
end

%% See if pattern has changed

normPatt = zeros(sizeSeq,sizeSeq);

for i = 1:sizeSeq
   for j = i:sizeSeq
      normPatt(i,j) = norm(logical(A{j})-logical(A{i}),'f');
   end
end
%% Compute Sam

tic
pre_mex1 = MEX_sam_preprocess_left(A{1});
time_prepr(1) = toc;
tic
pre_mex2 = MEX_sam_preprocess_left(A{1},A{1}^2);
time_prepr(2) = toc;

for i = 1:sizeSeq
   tic
   sam1{i} = MEX_sam_compute_left(A{i},A{1},pre_mex1);
   time_compute1(i) = toc;
   normmSam1(i) = 2 * norm(A{1}-sam1{i}*A{i},'f')/(norm(A{1},'f')+norm(sam1{i}*A{i},'f'));

   tic
   sam2{i} = MEX_sam_compute_left(A{i},A{1},pre_mex2);
   time_compute2(i) = toc;
   normmSam2(i) = 2 * norm(A{1}-sam2{i}*A{i},'f')/(norm(A{1},'f')+norm(sam2{i}*A{i},'f'));
end

%% Compute the Preconditioner

gmres_restart = 100;
gmres_tol     = 1e-6;
gmres_maxit   = 10;

[M,time] = computeAMG(A{1},false);
t_amg0    = time;
P0_amg    = M;
fprintf('done in %.3f s\n\n', t_amg0);

%% Solve the preconditioner

for i = 1:sizeSeq

   uni = ones(size(A{i},1),1);

   if i > 1
      [MM,time] = computeAMG(A{i},false,0);
   else
      MM = P0_amg;
   end
   t0 = tic;
   [~, ~, ~, it_recomp(i,:)] = gmres_RIGHT(@(v) A{i}*v, uni, gmres_restart, gmres_tol, gmres_maxit, MM, []);
   t_solve_recomp(i) = toc(t0);
   it(i,1) = (it_recomp(i,1)-1)*gmres_restart + it_recomp(i,2);

   Pfun_new = @(v) sam_apply_left(speye(size(A{i},1)), P0_amg, v);
   t0 = tic;
   [~, ~, ~, it_reuse(i,:)] = gmres_RIGHT(@(v) A{i}*v, uni, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
   t_solve_reuse(i) = toc(t0);
   it(i,2) = (it_reuse(i,1)-1)*gmres_restart + it_reuse(i,2);
   
   Pfun_new = @(v) sam_apply_left(sam1{i}, P0_amg, v);
   t0 = tic;
   [~, ~, ~, it_sam1(i,:)] = gmres_RIGHT(@(v) A{i}*v, uni, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
   t_solve_sam1(i) = toc(t0);
   it(i,3) = (it_sam1(i,1)-1)*gmres_restart + it_sam1(i,2);
   
   
   Pfun_new = @(v) sam_apply_left(sam2{i}, P0_amg, v);
   t0 = tic;
   [~, ~, ~, it_sam2(i,:)] = gmres_RIGHT(@(v) A{i}*v, uni, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
   t_solve_sam2(i) = toc(t0);
   it(i,4) = (it_sam2(i,1)-1)*gmres_restart + it_sam2(i,2);
end


%%

hold on
for i = 1:4
   plot(1:sizeSeq, it(:,1), 'b-o', 1:sizeSeq, it(:,2), 'g-s', 1:sizeSeq, it(:,3), 'm-d', 1:sizeSeq, it(:,4), 'c-^');
   xlabel('Index');
   ylabel('Iterations');
   grid on;
   ylim([0 150])
end

legend('RecompPrec','ReusePrec1','SAM1+P1','SAM2+P1','Location','best')





























%% Check shifted fsai

alpha = 1.e-1;
alpha = alpha*mean(diag(A{1}));
cheap = false;
nn = size(A{1},1);
ID = speye(nn);
AA = A{1} + alpha*ID;


normWOSAM = 2*norm(AA-A{1},'f')/(norm(A{1},'f')+norm(AA,'f'));

[samShift1,~] = MEX_sam_compute_left(AA,A{1},pre_mex1);
[samShift2,~] = MEX_sam_compute_left(AA,A{1},pre_mex2);

normSamShift1 = 2*norm(samShift1*AA-A{1},'f')/(norm(A{1},'f')+norm(samShift1*AA,'f'));
normSamShift2 = 2*norm(samShift2*AA-A{1},'f')/(norm(A{1},'f')+norm(samShift2*AA,'f'));

% semilogy(0:2,[normWOSAM,normSamShift1,normSamShift2],'r*-')

B = [AA ID; ID ID/alpha];

[GG,GG1] = NSY_rfsai_cpp(15,1,1e-3,B);
[GG2,GG3] = NSY_rfsai_cpp(15,1,1e-3,A{1});
R = [ID sparse(nn,nn)];
P = [ID; sparse(nn,nn)];
MAfun = @(x) R*(GG1*(GG*(P*(A{1}*x))));
lambda = eigs(MAfun,size(AA,1),7,'lm','Tolerance',1.e-5,'Display',1,...
              'FailureTreatment','keep');

lambda1 = eigs(@(x) A{1}*x,size(A{1},1),7,'lm','Tolerance',1.e-5,'Display',1,...
              'FailureTreatment','keep');

ratio(1) = max(lambda1)/min(lambda1);
ratio(2) = max(lambda)/min(lambda);

plot(lambda,'r*')
hold on
plot(lambda1,'k*')

Pfun_new = @(v) sam_apply_left(samShift1,@(x) R*(GG1*(GG*(P*x))),v);
Pfun_new2 = @(v) sam_apply_left(samShift2,@(x) R*(GG1*(GG*(P*x))),v);
Pfun_new3 = @(x) GG3*(GG2*x);

[~, ~, ~, it_sa(1,:)] = gmres_RIGHT(@(v) A{i}*v, ones(size(AA,1),1), 100, 1e-6, 10, Pfun_new, []);
[~, ~, ~, it_sa(2,:)] = gmres_RIGHT(@(v) A{i}*v, ones(size(AA,1),1), 100, 1e-6, 10, Pfun_new2, []);
[~, ~, ~, it_sa(3,:)] = gmres_RIGHT(@(v) A{i}*v, ones(size(AA,1),1), 100, 1e-6, 10, Pfun_new3, []);