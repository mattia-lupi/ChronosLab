function [Pout,info] = MEX_EMIN_enhance_prol(level,param,A,Ppatt,Pin,TV,fcnode)

%-----------------------------------------------------------------------------------------
%
% Improves an input prolongation using Energy Minimization, enforcing the sparsity pattern
% of P_patt and ensuring representation of the test space
%
% Input
%
% level:     : level number (starting from 1)
% param      : MATLAB structure storing EMIN parameters
%              np         ==> number of threads
%              itmax_emin ==> number of PCG iterations for Energy Minimization
%              prec_emin  ==> preconditioner for energy minimization
%              solv_emin  ==> solution algorithm (Restr PCG, nullspace, matrix-free)
%              min_lfil   ==> minimum fill-in for Cholesky factorization
%              max_lfil   ==> maximum fill-in for Cholesky factorization
%              D_lfil     ==> fill-in increase for Cholesky factorization
% A          : matrix used to compute FSAI
% Ppatt      : prolongation pattern to be enforced
% Pin        : initial prolongation
% TV         : test space
% fcnode     : fine/coarse indicator
%
% Output
%
% Pout       : improved prolongation
% info       : array with timings and information on EMIN processin
%
%-----------------------------------------------------------------------------------------

% Unpack the input matrix
nn = size(A,1);
nt_A = nnz(A);
[iat_A,ja_A,coef_A] = unpack_csr(A);

% Unpack the input prolongation
nn_C = size(Pin,2);
nt_P = nnz(Pin);
[iat_Pin,ja_Pin,coef_Pin] = unpack_csr(Pin);

% Unpack the input pattern
nt_patt = nnz(Ppatt);
[iat_patt,ja_patt,~] = unpack_csr(Ppatt);

% Switch from column to row
iat_A     = iat_A';
ja_A      = ja_A';
coef_A    = coef_A';
iat_Pin   = iat_Pin';
ja_Pin    = ja_Pin';
coef_Pin  = coef_Pin';
iat_patt  = iat_patt';
ja_patt   = ja_patt';

% Convert TV in proper 1-D array
ntv = size(TV,2);
TV = TV';
TV = TV(:);

% Adpat fcnode to C++
fcnode(fcnode>0) =fcnode(fcnode>0) - 1;

% Switch from matlab type to c++ type
level     = int32(level);
np        = int32(param.np); 
itmax     = int32(param.itmax_emin); 
energ_tol = param.energ_tol; 
condmax   = param.condmax_emin; 
maxwgt    = param.maxwgt_emin; 
prec      = int32(param.prec_emin); 
sol_type  = int32(param.solv_emin); 
min_lfil  = int32(param.min_lfil); 
max_lfil  = int32(param.max_lfil); 
D_lfil    = int32(param.D_lfil); 
nn        = int32(nn);
nn_C      = int32(nn_C);
ntv       = int32(ntv);
nt_A      = int32(nt_A);
nt_P      = int32(nt_P);
nt_patt   = int32(nt_patt);
fcnode    = int32(fcnode);
iat_A     = int32(iat_A) - 1;
ja_A      = int32(ja_A) - 1;
iat_Pin   = int32(iat_Pin) - 1;
ja_Pin    = int32(ja_Pin) - 1;
iat_patt  = int32(iat_patt) - 1;
ja_patt   = int32(ja_patt) - 1;

% Compute Energy Min prolongation --------------------------------------------------------
[iat_Pout,ja_Pout,coef_Pout,info] = ...
          EMIN_Prolong_compute(level,np,itmax,energ_tol,condmax,maxwgt,prec,sol_type,...
                               min_lfil,max_lfil,D_lfil,nn,nn_C,ntv,nt_A,nt_P,nt_patt,...
                               fcnode,iat_A,ja_A,coef_A,iat_Pin,ja_Pin,coef_Pin,...
                               iat_patt,ja_patt,TV);

% Create a sparse matrices for Pout
nt_Pout = size(ja_Pout,2);
irow_Pout = zeros(nt_Pout,1);
iend   = iat_Pout(1)-1;
for i = 1:nn
   istart = iend + 1;
   iend = iat_Pout(i+1)-1;
   irow_Pout(istart:iend) = i;
end
ja_Pout = double(ja_Pout);
Pout = sparse(irow_Pout,ja_Pout,coef_Pout,nn,nn_C);

return
