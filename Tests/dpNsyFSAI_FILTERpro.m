
% Set reference number of iterations
refIter = 46;

% Load the matrix
name = fullfile("..","examples","matrices/","drackPrag.mat");
load(name);
A = Amat;
rhs = b;

% Read the default params
param = readDefaultParams();

% Modify eventual params to have the correct test
param.prolong.prol_emin = 'SMOOTH';
param.filter.filt_wgt = 85;
param.filter.filt_tol = 0.01;

% Check the symmetry of the matrix
param.symm = (norm(A-A','f')/norm(A,'f')<1e-14);

% Compute preconditioner for the test case with the correct params
[M, ~] = computeAMG(A, TV0, param, 0);

% Set Solver parameters
tol = 1e-6; % Tolerance for the solver
itmax = 100; % Maximum number of iterations
restart = 50; % Restart parameter for GMRES

% Solve the system
if param.symm
   [~,~,~,niter,~] = SQMR(@(v) A*v,rhs,tol,itmax,M,@(v) speye(size(A))*v,[],0);
else
   [~, ~, ~, it] = gmres_RIGHT(@(v) A*v, rhs, restart, tol, itmax/restart, M,[],[],0);
   niter = (it(1)-1)*restart + it(2);
end

if abs(refIter / niter - 1) > 0.05 && abs(refIter-niter) - 2 > 0
   fprintf('Converged in %d iterations\n', niter)
   error('Iterations do not coincide')
else
   fprintf('Test passed, converged in %d iterations\n', niter);
end