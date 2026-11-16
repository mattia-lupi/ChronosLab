function [P] = cpt_TwoStageProlMEX_ExtI(param,Clist,Fclist,Fflist,A,S_patt,A_comp)

% Count the entries in the lists
nC = numel(Clist);
nFc = numel(Fclist);
nFf = numel(Fflist);
nn = size(A, 1);

% Compute P1: interpolation from C+Fc to Ff
fcnode = -2 * ones(nn, 1);
fcnode(Clist) = 1;
fcnode(Fclist) = 1;
fcnode(Fflist) = -1;
[P1] = cpt_ProlMEX_EXTI(param,nC+nFc,fcnode,S_patt,A,A_comp);

% Compute P2: interpolation from C to Fc
fcnode = -2 * ones(nn, 1);
fcnode(Clist) = 1;
fcnode(Fclist) = -1;
fcnode(Fflist) = -2;
[P2] = cpt_ProlMEX_EXTI(param,nC,fcnode,S_patt,A,A_comp);

% Compute final interpolation
P = P1*P2;

return
