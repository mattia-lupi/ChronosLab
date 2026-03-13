function [x,iter,resvec,flag] = PCG_CJ(prodM,prodA,b,x0,itmax,tol);

nn = size(b,1);

% Init PCG
iter = 0;
flag = 0;
resvec = zeros(itmax,1);

% Compute initial residual
bnorm = norm(b);
if bnorm == 0.0
   x = zeros(nn,1);
   resvec = resvec(1);
   return;
end

% Alllocate some vectors
res = zeros(nn,1);
pres = zeros(nn,1);
p = zeros(nn,1);
axp = zeros(nn,1);

% Set initial solution
x = x0;

% Compute initial residual
res = b - prodA(x);
resiter = norm(res) / bnorm;
resvec(iter+1) = resiter;
exit_test = itmax <= 0;

% PCG loop
while ~exit_test

   iter = iter + 1;

   pres = prodM(res);

   % Compute beta
   if iter == 1
      p = pres;
   else
      beta = -(pres'*axp) / ptap;
      p = pres + beta*p;
   end

   % Compute alpha
   axp = prodA(p);
   ptap = p'*axp;
   alpha = (p'*res) / ptap;

   % Update solution and residual
   x = x + alpha*p;
   res = res - alpha*axp;

   % Check convergence
   resiter = norm(res) / bnorm;
   resvec(iter+1) = resiter;
   if mod(iter,10) == 0
      fprintf('%4d %15.6e\n',iter,resiter);
   end

   exit_test = resiter < tol || iter == itmax;

end
resvec = resvec(1:iter+1);

if resiter > tol
   flag = 1;
end

end
