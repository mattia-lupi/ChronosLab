function MEX_DumpChronos(mat_A,TV0)

nn_A = size(mat_A,1);
nt_A = nnz(mat_A);
[irow,jcol,coef] = find(mat_A);

ntv = size(TV0,2);
TV0 = TV0(:);

nn_A = int64(nn_A);
nt_A = int64(nt_A);
ntv = int64(ntv);
irow = int64(irow);
jcol = int64(jcol);

DumpChronos(nn_A,nt_A,ntv,irow,jcol,coef,TV0);

end
