function [fcnode,Clist,Fclist,Fflist,S,S_patt] = agg_coarsen(param,A,TV)
%-----------------------------------------------------------------------------------------
%
% Function to create the coarse space as a Fine/Coarse indicator and the lists of fine and
% coarse nodes
%
% Input
%
% param:    parameters controlling coarsening  .SoC_type ---> type of SoC (affinity,
%                                                             classical, classical NSY,
%                                                             diagonal dominance,
%                                                             algebraic distance
%                                              .tau ---> drop tolerance for strong
%                                                        connections or avg nonzero per
%
% A:           original system matrix
% smootherOp:  smoother operator
% TV:          test space to match
%
% Output
%
% fcnode:   Fine/Coarse indicator (-1 == fine, 1 == coarse)
% Clist:    list of coarse nodes
% Fclist:   list of fine/coarse nodes (coarse in the first stage)
% Fflist:   list of fine nodes (fine from the first stage)
% S:        strong connections matrix
% S_patt:   pattern of the strong connections
%
%-----------------------------------------------------------------------------------------

% Extract parameters
SoC_type = param.SoC_type;
tau      = param.tau;

% Initialize all nodes as free
CFc_Ff_node = zeros(size(A,1),1);

% Find Dirichlet nodes
[CFc_Ff_node,ndir] = mark_isolated_nodes(A,CFc_Ff_node);

switch upper(SoC_type)
   case 'AFF'
      % Compute strength of connections (use the pattern of A to evaluate affinities)
      [S,S_patt] = cpt_SoC_aff(tau,A,A,TV);
   case 'CLA'
      % Compute classicla strength of connections
      [S,S_patt] = cpt_SoC_cla(tau,A);
   case 'CLANSY'
      % Compute classicla strength of connections (non-symmetrize)
      [S,S_patt] = cpt_SoC_nsyCla(tau,A);
   case 'DOM'
      % Compute strength of connections based on diagonal dominance
      [S,S_patt] = cpt_SoC_dom(tau,A);
   case 'ALG'
      % Compute strength of connections based on diagonal dominance
      [S,S_patt] = cpt_SoC_AlgDist(tau,smootherOp,A);
   otherwise
      err_msg = [SoC_type ' is not a valid key for SoC'];
      error(err_msg);
end
fprintf('NNZR of original matrix:     %10.2f\n',nnz(A)/size(A,1)-1);
Stmp = S;
Stmp = Stmp - diag(diag(Stmp));
fprintf('NNZR of filtered SoC matrix: %10.2f\n',nnz(Stmp)/size(Stmp,1));

% Compute the Maximum Independent Set

% First stage (same of the usual coarsening)
[CFc_Ff_node, CFclist, Fflist] = cpt_PMIS(S, CFc_Ff_node);

% Condense the SoC graph for second stage of coarsening
S_CFc = mk_SoC_FcC(S, CFclist,CFc_Ff_node);

% Second stage (aggressive coarsening)
% Initialize all nodes as free
Fc_C_node = zeros(size(S_CFc,1),1);

% Compute PMIS
[Fc_C_node, Clist_loc, Fclist_loc] = cpt_PMIS(S_CFc, Fc_C_node);
Clist = CFclist(Clist_loc);
Fclist = CFclist(Fclist_loc);

% Compute fcnode
nC = numel(Clist);
nFc = numel(Fclist);
nFf = numel(Fflist);
fcnode = zeros(nC+nFc+nFf,1);
fcnode(Clist) = 1;
fcnode(fcnode>0) = 1:nC;
fcnode(fcnode==0) = -(1:nFc+nFf);

fprintf('# of Dirichlet nodes: %10d\n',ndir);

return
