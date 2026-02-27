function patt_out = drop_patt(nnzr,patt_in,patt_0);

[nrows,ncols] = size(patt_in);

% Use patt_0 to avoid dropping original entries
fac = 1.000*max(max(abs(patt_in)));
[ii,jj,pp] = find(patt_0);
pp = fac;
patt_0 = sparse(ii,jj,pp,nrows,ncols);

patt_in = patt_in + patt_0;
[iat,ja,coef] = unpack_csr(patt_in);

ntmax = nnzr*nrows;
ii_out = zeros(ntmax,1);
jj_out = zeros(ntmax,1);
pp_out = zeros(ntmax,1);
ind = 1;
iend = iat(1);
for irow = 1:nrows
   istart = iend;
   iend = iat(irow+1);
   if nnzr < (iend-istart)
      % Select the nnzr largest entries
      col_indices = ja(istart:iend-1);
      entries = coef(istart:iend-1);
      [~,perm] = sort(abs(entries),'descend');
      perm = perm(1:nnzr);
      ii_out(ind:ind-1+nnzr) = irow;
      jj_out(ind:ind-1+nnzr) = col_indices(perm);
      pp_out(ind:ind-1+nnzr) = entries(perm);
      ind = ind + nnzr;
   else
      % Leave the row unchanged
      k_out = iend-istart;
      ii_out(ind:ind-1+k_out) = irow;
      jj_out(ind:ind-1+k_out) = ja(istart:iend-1);
      pp_out(ind:ind-1+k_out) = coef(istart:iend-1);
      ind = ind + k_out;
   end
end

ii_out = ii_out(1:ind-1);
jj_out = jj_out(1:ind-1);
pp_out = pp_out(1:ind-1);
patt_out = sparse(ii_out,jj_out,pp_out,nrows,ncols);

% Remove patt_0
patt_out = patt_out - patt_0;

end
