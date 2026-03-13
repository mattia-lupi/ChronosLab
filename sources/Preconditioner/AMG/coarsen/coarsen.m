function [fcnode,clist,flist,S,S_patt] = coarsen(param,A,smootherOp,TV,verb)
%-----------------------------------------------------------------------------------------
%
% Function to create the coarse space as a Fine/Coarse indicator and the lists of fine and
% coarse nodes
%
% Input
%
% param:     parameters controlling coarsening .SoC_type ---> type of SoC (affinity,
%                                                             classical, classical NSY,
%                                                             diagonal dominance,
%                                                             algebraic distance
%                                              .tau ---> drop tolerance for strong
%                                                        connections or avg nonzero per
%                                                        row
% A:           original system matrix
% smootherOp:  smoother operator
% TV:          test space to match
%
% Output
%
% fcnode:   Fine/Coarse indicator (-1 == fine, 1 == coarse)
% clist:    list of coarse nodes
% flist:    list of fine nodes
% S:        strong connections matrix
% S_patt:   pattern of the strong connections
%
%-----------------------------------------------------------------------------------------

% Extract parameters
SoC_type = param.SoC_type;
tau      = param.tau;

% Initialize all nodes as free
fcnode = zeros(size(A,1),1);

% Find Dirichlet nodes
[fcnode,ndir] = mark_isolated_nodes(A,fcnode);

switch upper(SoC_type)
   case 'AFF'
      % Compute strength of connections (use the pattern of A to evaluate affinities)
      [S,S_patt] = cpt_SoC_aff(tau,A,A,TV,verb);
   case 'CLA'
      % Compute classicla strength of connections
      [S,S_patt] = cpt_SoC_cla(tau,A,verb);
   case 'CLANSY'
      % Compute classicla strength of connections (non-symmetrize)
      [S,S_patt] = cpt_SoC_nsyCla(tau,A,verb);
   case 'DOM'
      % Compute strength of connections based on diagonal dominance
      [S,S_patt] = cpt_SoC_dom(tau,A,verb);
   case 'ALG'
      % Compute strength of connections based on diagonal dominance
      [S,S_patt] = cpt_SoC_AlgDist(tau,smootherOp,A,verb);
   otherwise
      err_msg = [SoC_type ' is not a valid key for SoC'];
      error(err_msg);
end
if verb
   fprintf('NNZR of original matrix:     %10.2f\n',nnz(A)/size(A,1)-1);
end
Stmp = S;
Stmp = Stmp - diag(diag(Stmp));
if verb
   fprintf('NNZR of filtered SoC matrix: %10.2f\n',nnz(Stmp)/size(Stmp,1));
end

% Compute the Maximum Independent Set
[fcnode, clist, flist] = cpt_PMIS(S, fcnode,verb);

if verb
   fprintf('# of Dirichlet nodes: %10d\n',ndir);
end

return
