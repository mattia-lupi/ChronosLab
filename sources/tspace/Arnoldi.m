function [vec, val, iter, resid] = Arnoldi(V0,ProdMat,m,ntv,dual_orth)

n = size(V0,1);
H = zeros(m);
V = zeros(n,m);

% Compress input test space
startvec = sum(V0,2);
V(:,1) = normc(startvec);

for i = 1:m

   w = ProdMat(V(:,i));
   % First orthogonalization
   for j = 1:i
      H(j,i) = w'*V(:,j);
      w = w - H(j,i)*V(:,j);
   end
   % Second orthogonalization
   if dual_orth
      for j = 1:i
         fac = w'*V(:,j);
         H(j,i) = H(j,i) + fac;
         w = w - fac*V(:,j);
      end
   end
   if i == m; break; end
   H(i+1,i) = norm(w);
   V(:,i+1) = w/H(i+1,i);

end

% Compute eigenpairs of H
[Vr,lambda] = eig(H);
val = diag(lambda);
vec = V*Vr;

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

end
