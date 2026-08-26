function [pnew,info,relres,niter,resvec] = SQMR(Afun,rhs,tol,itmax,lprec,rprec,x0,verb)

%%%%%%%%%%%%%%%%%%%%%
%% VETTORI SCRATCH
nn   = size(rhs,1);
d    = zeros(nn,1);
%%%%%%%%%%%%%%%%%%%%%

if nargin >= 7
   pnew = x0;
else
   pnew = zeros(nn,1);
end

if nargin < 8
   verb = true;
end

clear resvec;
resvec = zeros(1,itmax);
info=0;
exit_test = false;
niter=0;
bnorm = norm(rhs);

% Initial residual: r
r = rhs - Afun(pnew);

% Preconditioned residual: q
vscr = lprec(r);
q = rprec(vscr);

% Auxiliary variables
tau    = norm(vscr);
rho0   = r'*q;
theta0 = 0;

% Main loop
while ~exit_test && niter < itmax

   niter = niter + 1;

   % alpha = rho_0 / (q^T A q)
   vscr = Afun(q);
   sigma = q'*vscr;
   if abs(sigma) < 1e-20*rho0
      fprintf('Small sigma %e\n',sigma)
      info=10;
      break;
   end
   alpha=rho0/sigma;

   % r <-- r - alpha A q
   r = r - alpha*vscr;

   % Update theta, gamma e tau
   vscr   = lprec(r);
   theta1 = norm(vscr)/tau;
   gamma  = 1 / (sqrt(1+theta1*theta1));
   tau    = tau*theta1*gamma;

   % d <-- gamma^2 theta_0^2 d + gamma^2 alpha q
   d = d*gamma*gamma*theta0*theta0 + gamma*gamma*alpha*q;
   theta0=theta1;

   %  pnew <-- pnew + d
   pnew = pnew + d;

   % Check convergence
   resnorm = norm(r)/bnorm;
   if resnorm < tol
      % Exit main loop
      exit_test = true;
    else

      u = rprec(vscr);
      rho1 = r'*u;
      % beta = r^T L^-T L^-1 r / rho_0
      if abs(rho0) < 1e-20*rho1
         fprintf('Small rho0 %e\n',rho0)
         info=20;
         break
      end
      beta = rho1/rho0;
      rho0 = rho1;

      % q <-- u + beta q
      q = u + beta*q;
   end

   % Store relative residual
   resvec(niter) = resnorm;
   if mod(niter,10) == 0
      if verb
         fprintf('iter/res: %d %e\n',niter,resnorm);
      end
   end

end

if niter >= itmax
   info = 1;
end

% Compute real relative residual
relres = norm(rhs - Afun(pnew))/bnorm;

return
