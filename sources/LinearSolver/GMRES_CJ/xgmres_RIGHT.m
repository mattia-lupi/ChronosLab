%-----------------------------------------------------------------------------------------
%
% Simplified version of RIGHT preconditioned GMRES without restarting
%
%-----------------------------------------------------------------------------------------

function [x,iter,resvec] = xgmres_LEFT(A,rhs,tol,itmax,Mfun);

nn = size(A,1);

iter = 0;
resvec = zeros(itmax,1);
x_k = zeros(nn,1);
H = zeros(itmax+1,itmax);
V = zeros(nn,itmax);

r_k = rhs-A*x_k;
beta = norm(r_k);
v_k = r_k / beta;
V(:,1) = v_k;
bnorm = beta;

%norm(r_k)

while iter < itmax && norm(r_k) > tol*bnorm

   iter = iter + 1;

   w = A*Mfun(V(:,iter));

   for k = 1:iter
      H(k,iter) = w'*V(:,k);
      w = w - H(k,iter)*V(:,k);
   end
   
   H(iter+1,iter) = norm(w);
   w = w / H(k+1,iter);
   V(:,iter+1) = w;

   e_1 = zeros(iter+1,1);
   e_1(1) = 1;
   %%%%%%%%%%%%%%%%%%%%
   %size( H(1:iter+1,1:iter))
   %size(e_1)
   %beta
   %%%%%%%%%%%%%%%%%%%%
   y_k = H(1:iter+1,1:iter) \ (beta*e_1);
   x_k = Mfun(V(:,1:iter)*y_k);
   r_k = rhs-A*x_k;
   resvec(iter) = norm(r_k) / bnorm;

end

resvec = resvec(1:iter);
x = x_k;

end
