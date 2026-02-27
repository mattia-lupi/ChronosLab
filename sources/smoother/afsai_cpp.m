function F = afsai_cpp(M,nthread,nstep,step_size,eps)

%-----------------------------------------------------------------------------------------
%
% Computes the FSAI
%
% Input
%
% M         : gathered matrix used to compute FSAI
% nthread   : # of threads
% nstep     : # of iterative step
% step_size : # of terms added per step
% eps       : tolerance
%
% Output
%
% F : FSAI matrix
%
%-----------------------------------------------------------------------------------------


% Unpack the Master matrix
nrows_M = size(M,1);
nterm_M = nnz(M);
[iat_M,ja_M,coef_M] = unpack_csr(M);


% Switch from column to row
iat_M  = iat_M';
ja_M   = ja_M';
coef_M = coef_M';


% Switch from matlab type to c++ type
nthread   = int32(nthread); 
nstep     = int32(nstep); 
step_size = int32(step_size);
iat_M     = int64(iat_M);
ja_M      = int32(ja_M);


% Initialization
nrows_M = size(iat_M,2) - 1;
nterm_M = size(ja_M,2);
nrows = nrows_M;
tau = 0.;


% Compute FSAI ---------------------------------------------------------------------------
[nterm_G,iat_G,ja_G,coef_G] = compute_local_fsai(nthread,nstep,step_size,tau,eps,nrows, ...
                                                 nrows_M,nterm_M,iat_M,ja_M,coef_M);
 
% Create a sparse matrix for the FSAI
irow_G = zeros(nterm_G,1);
iend   = iat_G(1)-1;
for i = 1:nrows_M;
     istart = iend + 1;
     iend = iat_G(i+1)-1;
     irow_G(istart:iend) = i;
end
ja_G   = double(ja_G);
F      = sparse(irow_G,ja_G,coef_G);


% end function
return


