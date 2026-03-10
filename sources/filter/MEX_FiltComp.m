function mat_AC = MEX_FiltComp(mat_A,pattern_min,np,tau,TV);
%------------------------------------------------------------------------------------------
%
% Filters the coarse level operator and compensate deleted entries using C++ functions
%
% MEX Wrapper function
%
%------------------------------------------------------------------------------------------

% Unpack mat_A
nn_A = size(mat_A,1);
[iat_A,ja_A,coef_A] = unpack_csr(mat_A);

% Unpack minimal pattern
nt_patt = nnz(pattern_min);
[iat_patt,ja_patt,~] = unpack_csr(pattern_min);

% Create some extra parameter
ntv = size(TV,2);

% Convert reals into integers
iat_A = int32(iat_A) - 1;
ja_A = int32(ja_A) - 1;
iat_patt = int32(iat_patt) - 1;
ja_patt = int32(ja_patt) - 1;

% Convert TV in proper 1-D array
TV = TV';
TV = TV(:);

[nt_AC,iat_AC,ja_AC,coef_AC] = FilterComp_wrap(np,tau,nn_A,iat_A,ja_A,coef_A,...
                                               nt_patt,iat_patt,ja_patt,ntv,TV);

% Create a sparse matrix for the filtered operator
nn_AC = nn_A;
irow_AC = zeros(nt_AC,1);
iend = iat_AC(1)-1;
for i = 1:nn_AC
    istart = iend + 1;
    iend = iat_AC(i+1)-1;
    irow_AC(istart:iend) = i;
end
nn_AC = double(nn_AC);
mat_AC = sparse(irow_AC,ja_AC,coef_AC,nn_AC,nn_AC);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('max irow: %d %d\n',max(irow_I))
%fprintf('size(mat_I): %d %d\n',size(mat_I))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
