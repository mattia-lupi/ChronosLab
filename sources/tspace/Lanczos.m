% Lanczos with no reorthogonalization
%
% Input:
%
% ProdMat  : matrix by vector product
% m        : dimension of the Lanczos space
% startvec : starting vector
%
% Returns:
%
% val      : estimated eigenvalues
% vec      : estimated eigenvectors
% id_conv  : convergece estimator
%

function [vec,val,iter,resid] = Lanczos(V0,ProdMat,m,ntv)

% Compress input test space
startvec = sum(V0,2);

% Retrieve problem dimension
n = length(startvec);

% Normalize starting vector
v = normc(startvec);

% Start Lanczos iteration, m steps
a = zeros(m,1);
b = zeros(m-1,1);
V = zeros(n,m);
for k = 1:m
   % Store k-th vector
   V(:,k) = v;
   if k == 1
      r = ProdMat(v);
   else
      r = ProdMat(v) - b(k-1)*v2;
   end
   a(k) = v'*r;
   r = r - a(k)*v;
   b(k) = norm(r);
   v2 = v;
   v = 1/b(k)*r;
   % estimate |A|_2 by |T|_1
   if k == 1
      anorm = abs( a(1)+b(1) );
   else
      anorm = max( anorm, b(k-1)+abs(a(k))+b(k) );
   end;
end

% Compute eigenpairs of the tridiagonal matrix
TT = diag(a) + diag(b(1:end-1),1) + diag(b(1:end-1),-1);
[ritz_vec,ritz_val] = eig(TT);

% Sort ritz values in increasing order
ritz_val = diag(ritz_val);
[~,perm] = sort(ritz_val,'descend');
ritz_val = ritz_val(perm);
ritz_vec = ritz_vec(:,perm);

% Estimate eigenvectors
vec = V*ritz_vec;
val = ritz_val;

% Estimate residual
res = ProdMat(vec)-vec*diag(val);
lres = zeros(m,1);
for i = 1:m
   lres(i) = norm(res(:,i));
end
lres_rel = lres./val;
resid = lres_rel;

% Sort according to minimal val
[val, perm] = sort(val, 'ascend');
vec = vec(:, perm);
resid = resid(perm);

vec = vec(:, 1:ntv);
val = val(1:ntv);
resid = [nan(ntv,1), resid(1:ntv)];

% Just for compatibility
iter = m;

return
