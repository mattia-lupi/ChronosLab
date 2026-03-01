function mat_I = cpt_ProlMEX_Classical(params,nn_C,fcnodes,mat_S,mat_A)
%------------------------------------------------------------------------------------------
%
% Computes the prolongation using the Classical algorithm
%
% MEX Wrapper function
%
% Input
%
% params:   structure for CLASSICAL .np        --> number of omp threads
% nn_C:     number of COARSE nodes
% fcnodes:  FINE/COARSE partition
% mat_S:    SoC matrix
% mat_S:    system matrix
%
%------------------------------------------------------------------------------------------

% Extract necessary parameters from params
np        = params.np;

% Unpack the SoC matrix
nn_S = size(mat_S,1);
% nt_S = nnz(mat_S);
[~,~,coef_S] = unpack_csr(mat_S);

% Unpack the system matrix
nn_A = size(mat_A,1);
nt_A = nnz(mat_A);
[iat_A,ja_A,coef_A] = unpack_csr(mat_A);

% Create vecstart
vecstart = zeros(np+1,1);
bsize = floor(nn_S/np);
resto = mod(nn_S,np);
for i = 2:resto+1
   vecstart(i) = vecstart(i-1) + bsize + 1;
end
for i = resto+2:np+1
   vecstart(i) = vecstart(i-1) + bsize;
end

% Adjust numbering of COARSE nodes in fcnodes
cnt = 0;
for i = 1:nn_S
   if fcnodes(i) > 0
      fcnodes(i) = cnt;
      cnt = cnt + 1;
   else
      fcnodes(i) = -1;
   end
end

% Create some extra parameter
nn_I = nn_S;
nc_I = nn_C;

% Convert reals into integers
level = int32(0);
np = int32(np);
vecstart = int32(vecstart);
nn_A = int32(nn_A);
nt_A = int32(nt_A);
iat_A = int32(iat_A) - 1;
ja_A = int32(ja_A) - 1;
fcnodes = int32(fcnodes);
coef_S = int32(coef_S);
nn_I = int32(nn_I);
nc_I = int32(nc_I);

[nt_I,iat_I,ja_I,coef_I] = cpt_Prolongation_Classical(level,np,...
                           vecstart,nn_A,nt_A,iat_A,ja_A,coef_A,coef_S,...
                           fcnodes,nn_I,nc_I);

% Create a sparse matrix for the prolongation
irow_I = zeros(nt_I,1);
iend = iat_I(1)-1;
for i = 1:nn_I
    istart = iend + 1;
    iend = iat_I(i+1)-1;
    irow_I(istart:iend) = i;
end
nn_I = double(nn_I);
nc_I = double(nc_I);
mat_I = sparse(irow_I,ja_I,coef_I,nn_I,nc_I);

return
