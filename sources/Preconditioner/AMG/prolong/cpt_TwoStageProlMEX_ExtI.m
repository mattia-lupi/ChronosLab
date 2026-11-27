function [P] = cpt_TwoStageProlMEX_ExtI(param,Clist,Fclist,Fflist,A,S_patt)

nC  = numel(Clist);
nFc = numel(Fclist);
nn  = size(A, 1);
nt_A = nnz(A);
np  = param.np;

% Unpack S_patt and A once across both stages
assert(nnz(S_patt) == nt_A && size(S_patt,1) == nn, 'EXTI:patternMismatch', ...
    'S_patt and A must have identical size and sparsity pattern.');
[~,~,coef_S] = unpack_csr(S_patt);
[iat_A,ja_A,coef_A] = unpack_csr(A);
iat_A  = int32(iat_A) - 1;
ja_A   = int32(ja_A) - 1;
coef_S = int32(coef_S);

% Compute P1: interpolation from C+Fc to Ff
fcnode = -2 * ones(nn, 1);
fcnode(Clist) = 1;
fcnode(Fclist) = 1;
fcnode(Fflist) = -1;
c_mask = (fcnode > 0);
fcnode(c_mask) = (0:sum(c_mask)-1)';
nn_I1 = nC + nFc + sum(fcnode == -1);
P1 = cpt_Prolongation_EXTI(0,np,nn,nt_A,iat_A,ja_A,coef_A,coef_S,int32(fcnode),nn_I1,nC+nFc);

% Compute P2: interpolation from C to Fc
fcnode = -2 * ones(nn, 1);
fcnode(Clist) = 1;
fcnode(Fclist) = -1;
c_mask = (fcnode > 0);
fcnode(c_mask) = (0:sum(c_mask)-1)';
nn_I2 = nC + sum(fcnode == -1);
P2 = cpt_Prolongation_EXTI(0,np,nn,nt_A,iat_A,ja_A,coef_A,coef_S,int32(fcnode),nn_I2,nC);

P = P1 * P2;
end
