function [iat,ja,coef] = unpack_csr(matrix)

[irow,ja,coef] = find(matrix);
[irow,perm] = sort(irow);
ja = ja(perm);
coef = coef(perm);

nn = size(matrix,1);
nt = nnz(matrix);
iat = ones(nn+1,1);
irow_old = 0;
for k = 1:nt
   irow_new = irow(k);
   if  irow_new > irow_old
      for j = irow_old+1:irow_new
         iat(j) = k;
      end
      irow_old = irow_new;
   end
end
k = nt+1;
for j = irow_old+1:nn+1
   iat(j) = k;
end

end
