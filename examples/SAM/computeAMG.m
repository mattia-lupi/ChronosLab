function [M,T_setup] = computeAMG(A,TV0,sym_flag,verb)

% Read the parameters inside the input files
param = readDefaultParams();

% Set the symmetry flag
param.symm = sym_flag;

%%%%%%%%%%%%%
% Treat Boundary conditions 
warning('off', 'MATLAB:eigs:NotAllEigsConvKeep');
lmax = eigs(A,1,'lm','FailureTreatment','keep','Display',0,'Tolerance',0.001,'MaxIterations',3);
d = diag(A);
idx = (d == 1);
d(idx) = lmax/10;
A = spdiags(d, 0, A);
%%%%%%%%%%%%%

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
% Set timings
amg_times
T_smoo = 0;
T_LowRank = 0;
T_tspa = 0;
T_coar = 0;
T_CR = 0;
T_prol = 0;
T_FilCLev = 0;
T_FilProl = 0;
T_setup = 0;

% Compute the AMG hierarchy
time_start = tic;
AMG_prec = cpt_aspAMG(param,A,TV0,true);
T_setup = toc(time_start);

M = @(x) AMG_Vcycle(AMG_prec,A,x);

end
