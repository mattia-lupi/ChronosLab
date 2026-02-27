%-----------------------------------------------------------------------------------------
% Simultaneous Rayleigh quotient minimization through Conjugate Gradient
%
% This function computes the smallest eigenpairs of the eigenproblem:
%
% A*x = lambda*x
%
% Input variables:
%
% V0         initial approximation of the smallest eigenvectors
% ProdMat    function handle representing the matrix by vector product (A*x)
% itmax      maximum number of iteration
% tol        exit tolerance on the maximum relative eigenvalue residual
% ritz_freq  frequency of the Rayleigh Ritz orthonormalization
%
% Output variables:
%
% Zk         set of approximate smallest eignevectors (same number of columns as V0)
% qk         set of approximate smallest eignevalues
% iter       number of iterations performed
% resid      relative residual on eigenvalue and eigenvectors
%
%-----------------------------------------------------------------------------------------

function [Zk,qk,iter,resid] = SRQCG(V0,ProdMat,itmax,tol,ritz_freq)

% Get dimensions
[nn,mm] = size(V0);

% Orthogonalize initial test space
[Zk,~] = qr(V0,0);

% Init GCRay
iter = 0;
NRM_res = 2.*tol;

% Preallocate some array
a = zeros(mm,1);
b = zeros(mm,1);
c = zeros(mm,1);
d = zeros(mm,1);
e = zeros(mm,1);
f = zeros(mm,1);
beta = zeros(mm,1);
qk = zeros(mm,1);
Rk = zeros(nn,mm);
Pk = zeros(nn,mm);
Wk = zeros(nn,mm);
resid_vec = zeros(mm,1);
resid_lam = zeros(mm,1);

% Initialize Eig_CG
AZk = ProdMat(Zk);
for i = 1:mm
   f(i) = Zk(:,i)'*Zk(:,i);
   qk(i) = Zk(:,i)'*AZk(:,i) / f(i);
   Rk(:,i) = AZk(:,i) - qk(i)*Zk(:,i);
   Pk(:,i) = 2*Rk(:,i) / f(i);
end

while NRM_res > tol && iter < itmax

   iter = iter + 1;
   if mod(iter,floor(itmax/10)) == 0
      fprintf('SRQCG iter %5d out of %5d | NRM_res: %15.6e\n',iter,itmax,NRM_res)
   end
   APk = ProdMat(Pk);
   for i = 1:mm
      a(i) = Pk(:,i)'*AZk(:,i);
      b(i) = Pk(:,i)'*APk(:,i);
      c(i) = Pk(:,i)'*Zk(:,i);
      d(i) = Pk(:,i)'*Pk(:,i);
      e(i) = Zk(:,i)'*AZk(:,i);
   end

   Delta = (d.*e-b.*f).^2 - 4*(b.*c-a.*d).*(a.*f-c.*e);
   % Search minimum eigenvalue
   alpha = 0.5*( (d.*e-b.*f) + sqrt(Delta) ) ./ (b.*c-a.*d);

   % New approximation
   for i = 1:mm
      Zk(:,i) = Zk(:,i) + alpha(i)*Pk(:,i);
   end

   % Normalization step
   if mod(iter,ritz_freq) == 0
      % Search on Ritz space
      A_ritz = Zk'*ProdMat(Zk);
      M_ritz = Zk'*Zk;
      [V_r,~] = eigs(A_ritz,M_ritz,mm,'SM');
      Zk = Zk*V_r;
      for i = 1:mm
          Zk(:,i) = Zk(:,i)/norm(Zk(:,i));
      end
      f(:) = 1;
   else
      % Simply Normalize
      for i = 1:mm
         f(i) = Zk(:,i)'*Zk(:,i);
      end
   end

   AZk = ProdMat(Zk);
   qkm1 = qk;
   for i = 1:mm
      qk(i) = Zk(:,i)'*AZk(:,i) / f(i);
      Rk(:,i) = AZk(:,i) - qk(i)*Zk(:,i);
      Wk(:,i) = 2*Rk(:,i) / f(i);
      beta(i) = -Wk(:,i)'*APk(:,i) / b(i);
      Pk(:,i) = Wk(:,i) + beta(i)*Pk(:,i);
      resid_vec(i) = norm(Rk(:,i)) / norm(AZk(:,i));
      resid_lam(i) = abs(qkm1(i)-qk(i)) / qk(i);
   end
   NRM_res = max(resid_lam);
end

resid = [resid_lam,resid_vec];

return
