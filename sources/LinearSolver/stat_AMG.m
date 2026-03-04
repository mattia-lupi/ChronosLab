function [xk,flag,relres,iter,resvec] = stat_AMG(A,rhs,tol,itmax,AMG_prec)
%-----------------------------------------------------------------------------------------
%
%  Stationary AMG solver
%
%-----------------------------------------------------------------------------------------

% Initialize iteration
nn = size(A,1);
iter = 0;
xk = zeros(nn,1);
bnorm = norm(rhs);
if bnorm ~= 0
   relres = 1;
else
   relres = 0;
   resvec(1) = 0;
end
res = rhs;

% Stationary iteration loop
while iter < itmax && relres > tol

   % Increase iteration counter
   iter = iter + 1;

   % Apply the AMG V-cycle to update solution
   xk = xk + AMG_Vcycle(AMG_prec,A,res);

   % Compute residual and relative residual norm
   res = rhs - A*xk;
   relres = norm(res)/bnorm;

   % Store residual
   resvec(iter) = relres;

end

% Check convergence
if relres <= tol
   flag = 0;
else
   flag = 1;
end
