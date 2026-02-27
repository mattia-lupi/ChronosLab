function Af = Filter_Mat(perc,A)

[nn,mm] = size(A);

% Filter SoC
ncut = floor(nnz(A)*perc/100);
[ii,jj,ss] = find(A);
[~,perm] = sort(abs(ss),'descend');
ii = ii(perm(1:ncut));
jj = jj(perm(1:ncut));
sf = ss(perm(1:ncut));
Af = sparse(ii,jj,sf,nn,mm);

end
