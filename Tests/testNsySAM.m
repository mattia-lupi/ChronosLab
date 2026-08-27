
% Set reference number of iterations
refIter = 111;

% Load the matrix
name = fullfile("..","examples","matrices","stickSlipFast","stickSlipFast_1.mat");
load(name);
A0 = Amat;
TV00 = TV0;
name = fullfile("..","examples","matrices","stickSlipFast","stickSlipFast_2.mat");
load(name);
A = Amat;
rhs = b;

% Read the default params
param = readDefaultParams();

% Modify params to suite fluid dynamics
param.prolong.prol_emin = 'EMIN';

% Check the symmetry of the matrix
param.symm = (norm(A-A','f')/norm(A,'f')<1e-14);

% Compute preconditioner for the test case with the correct params
[M, ~] = computeRACP(A0, TV00, param, 0);

% Set Solver parameters
tol = 1e-6; % Tolerance for the solver
itmax = 300; % Maximum number of iterations
restart = 50; % Restart parameter for GMRES

% Compute the Sparse approximate mapping
[sam, ~] = MEX_sam_adaptive_left(A, A0, maxNumCompThreads, 5, 1, 1e-3);

% To see the effect substitute sam with speye(size(A,1)) in the line
% below. The system will converge after at least 2000 iterations.
M_sam = @(v) sam_apply_left(sam, M, v);

% Solve the system
if param.symm
   [~,~,~,niter,~] = SQMR(@(v) A*v,rhs,tol,itmax,M_sam,@(v) speye(size(A))*v,[],0);
else
   [~, ~, ~, it] = gmres_RIGHT(@(v) A*v, rhs, restart, tol, itmax/restart, M_sam,[],[],0);
   niter = (it(1)-1)*restart + it(2);
end

if abs(refIter / niter - 1) > 0.1 && abs(refIter-niter) - 2 > 0
   fprintf('Converged in %d iterations\n', niter)
   error('Iterations do not coincide')
else
   fprintf('Test passed, converged in %d iterations\n', niter);
end