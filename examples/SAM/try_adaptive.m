
clc;
clear;
%%

load mats/StickSlipOpenBigFull/StickSlipOpenBigFull_1.mat
A0 = speye(size(Amat,1));
A0 = Amat;
load mats/StickSlipOpenBigFull/StickSlipOpenBigFull_22.mat
% Amat = A;

nstep = 5;
stp_size = 1;
eps = 1e-4;

N = [];
% tic
% [N,avg_resnormL1] = sam_adaptive_left(Amat',A0',nstep,stp_size,eps);
% t1 = toc;

tic;
[N1, res_norm] = MEX_sam_adaptive_left(Amat,A0,8,nstep,stp_size,eps);
t2 = toc;

M = [];
tic
pre_mex1 = MEX_sam_preprocess_left(Amat);
[M,avg_resnormS] = MEX_sam_compute_left(Amat,A0,pre_mex1);
t3 = toc;


%%
normSAM(Amat,A0)
res_norm
normSAM(M*Amat,A0)

% figure
% spy(N')
figure
spy(N1)

figure
spy(M)


%%
gmres_restart = 100;
gmres_tol     = 1e-6;
gmres_maxit   = 20;
clear itt;

[amg,time(1)] = computeRACP(A0,TV0,false);
MM = @(v) sam_apply_left(N1, amg, v);
[~, ~, ~, it] = gmres_RIGHT(@(v) Amat*v, ones(size(Amat,1),1), gmres_restart, gmres_tol, gmres_maxit, MM, []);
itt(1) = (it(1)-1)*gmres_restart + it(2);

if ~isempty(N)
   MM = @(v) sam_apply_left(N', amg, v);
   [~, ~, ~, it] = gmres_RIGHT(@(v) Amat*v, ones(size(Amat,1),1), gmres_restart, gmres_tol, gmres_maxit, MM, []);
   itt(end+1) = (it(1)-1)*gmres_restart + it(2);
end

if ~isempty(M)
   MM1 = @(v) sam_apply_left(M, amg, v);
   [~, ~, ~, it] = gmres_RIGHT(@(v) Amat*v, ones(size(Amat,1),1), gmres_restart, gmres_tol, gmres_maxit, MM1, []);
   itt(end+1) = (it(1)-1)*gmres_restart + it(2);
end

% riciclo
% MM2 = @(v) sam_apply_left(speye(size(A0,1)), amg, v);
% [~, ~, ~, it] = gmres_RIGHT(@(v) Amat*v, ones(size(Amat,1),1), gmres_restart, gmres_tol, gmres_maxit, MM2, []);
% itt(end+1) = (it(1)-1)*gmres_restart + it(2);

[amg_cor,time(1)] = computeRACP(Amat,TV0,false);
[~, ~, ~, it] = gmres_RIGHT(@(v) Amat*v, ones(size(Amat,1),1), gmres_restart, gmres_tol, gmres_maxit, amg_cor, []);
itt(end+1) = (it(1)-1)*gmres_restart + it(2);


%%



[iat,ja,coef] = unpack_csc(Amat);
iat = iat - 1;
ja = ja - 1;

writematrix(iat,"iat.dat")
writematrix(ja,"ja.dat")
writematrix(coef,"coef.dat")

[~,~,coefT] = unpack_csr(Amat);
writematrix(coefT,"coefT.dat")

[iat1,ja1,coef1] = unpack_csr(A0);
iat1 = iat1 - 1;
ja1 = ja1 - 1;

writematrix(iat1,"iat2.dat")
writematrix(ja1,"ja2.dat")
writematrix(coef1,"coef2.dat")

%%
tic
[~,~,coef] = unpack_csc(Amat);
[~,~,coef] = unpack_csr(Amat);
toc
