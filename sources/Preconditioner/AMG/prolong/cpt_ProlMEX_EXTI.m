function mat_I = cpt_ProlMEX_EXTI(params,nn_C,fcnodes,mat_S,mat_A)
%------------------------------------------------------------------------------------------
%
% Computes the prolongation using the Extended + I algorithm (see De Sterk et al.,
% "Distance-two interpolation for parallel algebraic multigrid", Numerical Linear Algebra
% with applications,15, pp. 115-139, 2008).
%
% MEX Wrapper function
%
% Input
%
% params:     structure for EXTI .np        --> number of omp threads
% nn_C:       number of COARSE nodes
% fcnodes:    FINE/COARSE partition
% mat_S:      SoC matrix
% mat_A:      system matrix
%
%------------------------------------------------------------------------------------------

np   = params.np;
nn_S = size(mat_S,1);
[~,~,coef_S] = unpack_csr(mat_S);

nn_A = size(mat_A,1);
nt_A = nnz(mat_A);
[iat_A,ja_A,coef_A] = unpack_csr(mat_A);

assert(nnz(mat_S) == nt_A && nn_S == nn_A, 'EXTI:patternMismatch', ...
    'mat_S and mat_A must have identical size and sparsity pattern.');

c_mask = (fcnodes > 0);
fcnodes(c_mask) = (0:sum(c_mask)-1)';

nn_Fc = sum(fcnodes == -1);
nn_I = nn_C + nn_Fc;

mat_I = cpt_Prolongation_EXTI(0,np,nn_A,nt_A,int32(iat_A)-1,int32(ja_A)-1,...
                              coef_A,int32(coef_S),int32(fcnodes),nn_I,nn_C);
end
