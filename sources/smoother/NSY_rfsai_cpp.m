function [FL,FU] = NSY_rfsai_cpp(nstep,step_size,epsilon,A)

%-----------------------------------------------------------------------------------------
%
% Computes the Adaptive NSY FSAI using a residual-based approach
%
% Input
%
% nstep      : number of steps of the adaptive procedure
% step_size  : number of entries added for each step
% epsilon    : exit tolerance from the adaptive procedure
% A          : matrix used to compute FSAI
%
% Output
%
% FL : Lower FSAI factor
% FU : Upper FSAI factor
%
%-----------------------------------------------------------------------------------------

% Unpack the input matrix
nn_A = size(A,1);
nt_A = nnz(A);
[iat_A,ja_A,coef_A] = unpack_csr(A);

% Switch from column to row
iat_A  = iat_A';
ja_A   = ja_A';
coef_A = coef_A';

% Switch from matlab type to c++ type
nstep     = int32(nstep); 
step_size = int32(step_size); 
nn_A      = int32(nn_A);
nt_A      = int32(nt_A);
iat_A     = int32(iat_A) - 1;
ja_A      = int32(ja_A) - 1;

% Compute NSY r-adaptive FSAI ------------------------------------------------------------
[iat_FL,ja_FL,coef_FL,iat_FU,ja_FU,coef_FU] = ...
      NSY_rFSAI_compute(nstep,step_size,epsilon,nn_A,iat_A,ja_A,coef_A);

% Create a sparse matrices for FL
nt_FL = size(ja_FL,2);
irow_FL = zeros(nt_FL,1);
iend   = iat_FL(1)-1;
for i = 1:nn_A
   istart = iend + 1;
   iend = iat_FL(i+1)-1;
   irow_FL(istart:iend) = i;
end
ja_FL = double(ja_FL);
FL = sparse(irow_FL,ja_FL,coef_FL);

% Create a sparse matrices for FU
nt_FU = size(ja_FU,2);
irow_FU = zeros(nt_FU,1);
iend   = iat_FU(1)-1;
for i = 1:nn_A
   istart = iend + 1;
   iend = iat_FU(i+1)-1;
   irow_FU(istart:iend) = i;
end
ja_FU = double(ja_FU);
FU = sparse(irow_FU,ja_FU,coef_FU);

% end function
return
