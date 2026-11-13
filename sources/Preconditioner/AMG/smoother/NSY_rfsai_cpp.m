function [FL,FU] = NSY_rfsai_cpp(nstep,step_size,epsilon,A,nthreads)

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
% nthreads   : number of threads to be used for openmp
%
% Output
%
% FL : Lower FSAI factor
% FU : Upper FSAI factor
%
%-----------------------------------------------------------------------------------------

% Unpack the input matrix
nn_A = size(A,1);
% nt_A = nnz(A);
[iat_A,ja_A,coef_A] = unpack_csr(A);

% Switch from column to row
iat_A  = iat_A';
ja_A   = ja_A';
coef_A = coef_A';

% Switch from matlab type to c++ type
iat_A     = int32(iat_A) - 1;
ja_A      = int32(ja_A) - 1;

% Compute NSY r-adaptive FSAI ------------------------------------------------------------
[iat_FL,ja_FL,coef_FL,iat_FU,ja_FU,coef_FU] = ...
      NSY_rFSAI_compute(nstep,step_size,epsilon,nn_A,iat_A,ja_A,coef_A,nthreads);

% Create a sparse matrices for FL
irow_FL = repelem((1:nn_A)', diff(iat_FL));
ja_FL = double(ja_FL);
FL = sparse(irow_FL, ja_FL, coef_FL, nn_A, nn_A);

% Create a sparse matrices for FU
irow_FU = repelem((1:nn_A)', diff(iat_FU));
ja_FU = double(ja_FU);
FU = sparse(irow_FU, ja_FU, coef_FU, nn_A, nn_A);

% end function
return
