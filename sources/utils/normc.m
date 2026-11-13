function V = normc(V)

nrm = vecnorm(V, 2, 1);
nrm(nrm == 0) = 1;
V = V ./ nrm;

return
