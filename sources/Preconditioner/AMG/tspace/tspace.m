% Compute the near null kernel of the preconditioned matrix F*A*F', starting from V0
% input: A, F, V0, param
% output: V, lambda, res

function [V, lambda, res] = tspace(V0_input, A, smootherOp, param, verb)

if nargin < 5
   verb = 1;
end

% Mandatory parameters
ntv         = param.ntv;
itmax       = param.itmax;
tol         = param.tol;
method      = param.method;
init_approx = param.init_approx;

% Optional parameters
if (isfield(param, 'ritz_freq'))
   ritz_freq = param.ritz_freq;
else
   ritz_freq = 1;
end
if (isfield(param, 'dual_orth'))
   dual_orth = param.dual_orth;
else
   dual_orth = true;
end

% Initial test space
n = size(A,1);
if (isempty(V0_input))
   rand('state', 0);
   V0 = rand(n, ntv);
else
   n1 = size(V0_input, 1);
   n2 = size(V0_input, 2);
   if (n1 ~= n)
       error('Initial near null kernel: wrong dimensions');
   end
   if (n2 < ntv)
      V0 = [V0_input, rand(n, ntv-n2)];
   elseif (n1 > ntv)
      V0 = V0_input(:, 1:ntv);
   else
      V0 = V0_input;
   end
end

switch lower(method)
   case 'none'  % Non-Generalized SRQCG
      % Just orthonormalize initial V0
      [V,~] = qr(V0,0);
      lambda = zeros(ntv,1);
      res    = zeros(ntv,1);
   case 'smoothing'  % Simple smoothing of the input V0
      % Just orthonormalize initial V0
      smooth = @(x) x - smootherOp.omega*smootherOp.right*(smootherOp.left*(A*x));
      [V, iter, res] = simple_smoothing(V0,smooth,itmax,tol,verb);
      [V,~] = qr(V,0);
      lambda = zeros(ntv,1);
   case 'ng-srqcg'  % Non-Generalized SRQCG
      % Define the matrix vector operation
      ProdMat = @(x) A*x;
      % Define the preconditioner applicarion
      Prec = @(x) smootherOp.right*(smootherOp.left*x);
      % Apply 1 step of smoothing to remove Dirichlet nodes
      V0 = V0 - Prec(ProdMat(V0));
      % Orthonormalize
      [V0,~] = qr(V0,0);
      % Call SRQCG
      [V, lambda, iter, res] = NG_SRQCG(V0, ProdMat, Prec, itmax, tol, ritz_freq, verb);
      % NON DOVREBBE SERVIRE--->  % Orthonormalize
      %                     --->  [V,~] = qr(V,0);
      V = normc(V);
   case 'srqcg'
      % Define preconditioned matrix vector operation
      ProdMat = @(x) smootherOp.left*(A*(smootherOp.right*x));
      % Define initial test vector space applying inv(F') or its approximation
      V0 = cpt_initApp(init_approx,smootherOp.right,V0,verb);
      % Apply 1 step of smoother to remove Dirichlet nodes
      V0 = V0 - ProdMat(V0);
      % Orthonormalize
      [V0,~] = qr(V0,0);
      % Call SRQCG
      [V, lambda, iter, res] = SRQCG(V0, ProdMat, itmax, tol, ritz_freq, verb);
      % Apply F'
      V = smootherOp.right*V;
      % Orthonormalize
      [V,~] = qr(V,0);
   case 'ng-lanczos'
      % Define preconditioned matrix vector operation
      ProdMat = @(x) A*x;
      % Apply 1 step of smoother to remove Dirichlet nodes
      V0 = V0 - smootherOp.right*(smootherOp.left*ProdMat(V0));
      % Orthonormalize
      [V0,~] = qr(V0,0);
      % Call Lanczos
      [V, lambda, iter, res] = Lanczos(V0, ProdMat, itmax, ntv, verb);
   case 'lanczos'
      % Define preconditioned matrix vector operation
      ProdMat = @(x) smootherOp.left*(A*(smootherOp.right*x));
      % Define initial test vector space applying inv(F') or its approximation
      V0 = cpt_initApp(init_approx,smootherOp.right,V0,verb);
      % Apply 1 step of smoother to remove Dirichlet nodes
      V0 = V0 - ProdMat(V0);
      % Orthonormalize
      [V0,~] = qr(V0,0);
      % Call Lanczos
      [V, lambda, iter, res] = Lanczos(V0, ProdMat, itmax, ntv, verb);
      % Apply F'
      V = smootherOp.right*V;
      % Orthonormalize
      [V,~] = qr(V,0);
   case 'arnoldi'
      % Define preconditioned matrix vector operation
      ProdMat = @(x) smootherOp.right*(smootherOp.left*(A*x));
      % Apply 1 step of smoother to remove Dirichlet nodes
      V0 = V0 - ProdMat(V0);
      % Orthonormalize
      [V0,~] = qr(V0,0);
      % Call Arnoldi
      [V, lambda, iter, res] = Arnoldi(V0, ProdMat, itmax, ntv, dual_orth, verb);
      % Orthonormalize
      [V,~] = qr(V,0);

   case 'block_arnoldi'
      ProdMat = @(x) smootherOp.right*(smootherOp.left*(A*x));
      % Apply 1 step of smoother to remove Dirichlet nodes
      V0 = V0 - ProdMat(V0);
      % Orthonormalize
      [V0,~] = qr(V0,0);
      % Call Arnoldi
      [iter, lambda, V, res_norm_X, res_norm_D, ~] = ...
         block_arnoldi(ProdMat, ntv, [], V0, 'smallest', itmax, tol, verb);
      % [V, lambda] = eigs(ProdMat, size(A,1), ntv,'smallestreal','Tolerance',1e-15, ...
      %    'MaxIterations',itmax, 'Display',1);
      res = [res_norm_D,res_norm_X];%zeros(ntv,2);%
      % lambda = diag(lambda);
      % V = real(V);
   otherwise
      error('Not existing method');
end

return
