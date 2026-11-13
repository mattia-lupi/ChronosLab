function [V, iter, res] = simple_smoothing(V0,smooth,itmax,tol,verb)

iter = 0;
res = 2*tol;
nrm_V = vecnorm(V0, 2, 1);
V = V0 ./ nrm_V;
nrm_Zold = nrm_V;
while (iter < itmax && max(res) > tol)
   iter = iter + 1;
   Z = smooth(V);
   nrm_Z = vecnorm(Z, 2, 1);
   res = abs(nrm_Z-nrm_Zold) ./ nrm_Zold;
   V = Z ./ nrm_Z;
   nrm_Zold = nrm_Z;
end
res = nrm_Z';

end
