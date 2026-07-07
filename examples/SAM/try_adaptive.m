
clc;
clear;

load mats/richards/richards_42.mat
A0 = speye(size(Amat,1));
% A0 = Amat;
load mats/richards/richards_73.mat

% [~,~,Aaug0] = computeRACP(A0,TV0,false);
% [~,~,Aaug1] = computeRACP(Amat,TV0,false);

nstep = 25;
stp_size = 1;
eps = 1e-5;

% profile off
% profile on
tic
[N,avg_resnormL] = sam_adaptive_left(Amat,A0,nstep,stp_size,eps);
% [N,avg_resnormL] = sam_adaptive_left(Aaug1',Aaug0',nstep,stp_size,eps);
t1 = toc;
Nleft = N';

tic
[N1, res_norm] = MEX_sam_adaptive_left(Amat,A0,1,nstep,stp_size,eps);
% [N,avg_resnormL] = sam_adaptive_left(Aaug1',Aaug0',nstep,stp_size,eps);
t2 = toc;
% Nleft1 = N1';
% tic
% [Nright,avg_resnormR] = sam_adaptive_left(Amat,A0,nstep,stp_size,eps);
% % [Nright,avg_resnormR] = sam_adaptive_left(Aaug1,Aaug0,nstep,stp_size,eps);
% t2 = toc;


% tic
% pre_mex1 = MEX_sam_preprocess_left(Amat);
% [M,avg_resnormS] = MEX_sam_compute_left(Amat,A0,pre_mex1);
% % pre_mex1 = MEX_sam_preprocess_left(Aaug1);
% % [M,avg_resnormS] = MEX_sam_compute_left(Aaug1,Aaug0,pre_mex1);
% t3 = toc;

% profile viewer

%%
2*norm(Amat-A0,'f')/(norm(A0,'f')+norm(Amat,'f'))
2*norm(Nleft*Amat-A0,'f')/(norm(A0,'f')+norm(Nleft*Amat,'f'))
2*norm(Amat*Nright-A0,'f')/(norm(A0,'f')+norm(Amat*Nright,'f'))
2*norm(M*Amat-A0,'f')/(norm(A0,'f')+norm(M*Amat,'f'))
% 2*norm(Aaug1-Aaug0,'f')/(norm(Aaug0,'f')+norm(Aaug1,'f'))
% 2*norm(Nleft*Aaug1-Aaug0,'f')/(norm(Aaug0,'f')+norm(Nleft*Aaug1,'f'))
% 2*norm(Aaug1*Nright-Aaug0,'f')/(norm(Aaug0,'f')+norm(Aaug1*Nright,'f'))
% 2*norm(M*Aaug1-Aaug0,'f')/(norm(Aaug0,'f')+norm(M*Aaug1,'f'))

figure
spy(Nleft)
figure
spy(Nright)
figure
spy(M)



%%
gmres_restart = 100;
gmres_tol     = 1e-6;
gmres_maxit   = 20;

[amg,time(1)] = computeAMG(A0,ones(size(A0,1),1),false);
MM = @(v) sam_apply_left(Nleft, amg, v);
[~, ~, ~, it] = gmres_RIGHT(@(v) Amat*v, ones(size(Amat,1),1), gmres_restart, gmres_tol, gmres_maxit, MM, []);
itt(1) = (it(1)-1)*gmres_restart + it(2);

MM1 = @(v) sam_apply_right(Nright, amg, v);
[~, ~, ~, it] = gmres_RIGHT(@(v) Amat*v, ones(size(Amat,1),1), gmres_restart, gmres_tol, gmres_maxit, MM1, []);
itt(2) = (it(1)-1)*gmres_restart + it(2);

% riciclo
MM2 = @(v) sam_apply_left(speye(size(A0,1)), amg, v);
[~, ~, ~, it] = gmres_RIGHT(@(v) Amat*v, ones(size(Amat,1),1), gmres_restart, gmres_tol, gmres_maxit, MM2, []);
itt(3) = (it(1)-1)*gmres_restart + it(2);

[amg_cor,time(1)] = computeAMG(Amat,ones(size(A0,1),1),false);
[~, ~, ~, it] = gmres_RIGHT(@(v) Amat*v, ones(size(Amat,1),1), gmres_restart, gmres_tol, gmres_maxit, amg_cor, []);
itt(4) = (it(1)-1)*gmres_restart + it(2);

%%


[iat,ja,coef] = unpack_csr(Amat);
iat = iat - 1;
ja = ja - 1;

writematrix(iat,"iat.dat")
writematrix(ja,"ja.dat")
writematrix(coef,"coef.dat")