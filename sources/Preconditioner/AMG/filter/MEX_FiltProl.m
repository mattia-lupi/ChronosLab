%-----------------------------------------------------------------------------------------
%
% Filters the prolongation operator ensuring the same action on the test space, that is
% PF*TV = P*TV, where TV is the test space in the next level.
%
% Each row is filtered until the absolute norm of the remaining entries equals perc
% percentage of the original absolute norm.
%
% If the final norm of the row differes more that tol from the original row, more entries
% are allowed in the row.
%
% MEX Wrapper function
%
%-----------------------------------------------------------------------------------------

function mat_PF = Filter_Prol(np,perc,tol,TV,mat_P)

% Unpack mat_P
[nn_P,nc_P] = size(mat_P);
[iat_P,ja_P,coef_P] = unpack_csr(mat_P);

% Create some extra parameter
[nr_TV,ntv] = size(TV);

% Convert reals into integers
np = int32(np);
nn_P = int32(nn_P);
iat_P = int32(iat_P) - 1;
ja_P = int32(ja_P) - 1;
ntv = int32(ntv);

% Convert TV in proper 1-D array
TV = TV';
TV = TV(:);

% Call the wrapper function
[nt_PF,iat_PF,ja_PF,coef_PF] = FilterProl_wrap(np,perc,tol,...
                                               nn_P,iat_P,ja_P,coef_P,nr_TV,ntv,TV);

% Create a sparse matrix for the filtered operator
nn_PF = nn_P;
irow_PF = zeros(nt_PF,1);
iend = iat_PF(1)-1;
for i = 1:nn_PF
    istart = iend + 1;
    iend = iat_PF(i+1)-1;
    irow_PF(istart:iend) = i;
end
nn_PF = double(nn_PF);
nc_PF = nc_P;
mat_PF = sparse(irow_PF,ja_PF,coef_PF,nn_PF,nc_PF);

end
