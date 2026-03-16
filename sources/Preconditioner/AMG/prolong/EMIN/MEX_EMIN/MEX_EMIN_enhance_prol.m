function [Pout,info] = MEX_EMIN_enhance_prol(level,param,A,Ppatt,Pin,TV,fcnode,verb)

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

PRINT_DEBUG = false;
if PRINT_DEBUG
   print_SpMat(strcat('matrix_',num2str(level),'.csr'),A);
   print_SpMat(strcat('pattern_',num2str(level),'.csr'),Ppatt);
   print_SpMat(strcat('prolong_',num2str(level),'.csr'),Pin);
   fid = fopen(strcat('fcnode_',num2str(level)),'w');
   for i = 1 : length(fcnode)
      fprintf(fid,'%10d\n',fcnode(i));
   end
   fclose(fid);
   fid = fopen(strcat('TV_',num2str(level)),'w');
   fprintf(fid,'%10d %10d\n',size(TV,1),size(TV,2));
   for i = 1 : size(TV,1)
      for j = 1 : size(TV,2)
         fprintf(fid,'%25.15e',TV(i,j));
      end
      fprintf(fid,'\n');
   end
   fclose(fid);
end

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
np        = param.np; 
itmax     = param.itmax_emin; 
energ_tol = param.energ_tol; 
condmax   = param.condmax_emin; 
prec      = param.prec_emin; 
sol_type  = param.solv_emin; 
nn        = nn;
nn_C      = nn_C;
ntv       = ntv;
nt_A      = nt_A;
nt_P      = nt_P;
nt_patt   = nt_patt;
fcnode    = int32(fcnode);
iat_A     = int32(iat_A) - 1;
ja_A      = int32(ja_A) - 1;
iat_Pin   = int32(iat_Pin) - 1;
ja_Pin    = int32(ja_Pin) - 1;
iat_patt  = int32(iat_patt) - 1;
ja_patt   = int32(ja_patt) - 1;

% Compute Energy Min prolongation --------------------------------------------------------
[iat_Pout,ja_Pout,coef_Pout,info] = ...
          EMIN_Prolong_compute(level,np,itmax,energ_tol,condmax,prec,sol_type,...
                               nn,nn_C,ntv,nt_A,nt_P,nt_patt,fcnode,iat_A,ja_A,coef_A,...
                               iat_Pin,ja_Pin,coef_Pin,iat_patt,ja_patt,TV);

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
