function [S,S_patt] = cpt_SoC_aff(tau,A,S_patt,TV)
%-----------------------------------------------------------------------------------------
%
% Function to compute affinity-based strength of connection. This version also computes
% negative affinities to check that they are not taken into account.
%
% Input
%
% tau        >= 0 --->  drop tolerance to remove weak connections
%            <  0 --->  avg number of connections per row to keep in filtered SoC
% A          system matrix
% S_patt     sparsity pattern containing every possible connection
% TV         test space used
%
% Output
%
% S       SoC matrix
% S_patt  SoC matrix with 1 indicating strong connections and zero elsewhere
%
%-----------------------------------------------------------------------------------------

global DEBINFO

S_patt = tril(S_patt,-1);
nn = size(S_patt,1);
nt = nnz(S_patt);
TV = TV';
% Compute row affinity
[ii,jj,aa] = find(S_patt);
ss = zeros(nt,1);
for ind = 1:nt

   if mod(ind,floor(nt/10)) == 0
      fprintf('Processed %9i entries out of %9i\n',ind,nt)
   end

   v = TV(:,ii(ind));
   w = TV(:,jj(ind));
   %ss(ind) = cpt_sign_aff(v,w);
   ss(ind) = cpt_aff(v,w);

end
S = sparse(ii,jj,ss,nn,nn);

% Complete the matrix with a unitary diagonal
S = S + S' + speye(nn);

% Filter SoC sig retaining a given percentage of A off-diagonals
% Since S contains the diagonal, we retain nn more

% Filter SoC
[ii,jj,ss] = find(S);
fprintf('Filtering SoC\n');
if tau >= 0
   sf = ss;
   sf(sf <= tau)=0;
   S = sparse(ii,jj,sf,nn,nn);
   sp = ss;
   sp(ss <= tau)=-1;
   sp(ss >  tau)=1;
   S_patt = sparse(ii,jj,sp,nn,nn);
else
   ncut = ceil(-nn*tau);
   if ncut<nnz(A)
      [ss,perm] = sort(ss,'descend');
      lim_val = ss(ncut);
      ncut = find(ss<lim_val,1,'first')-1;
      ii_f = ii(perm(1:ncut));
      jj_f = jj(perm(1:ncut));
      sf = ss(1:ncut);
      S = sparse(ii_f,jj_f,sf,nn,nn);
      ss(:) = 1;
      ss(ncut+1:end) = -1;
      S_patt = sparse(ii,jj,ss,nn,nn);
   else
      S_patt = S;
      S_patt(S_patt>0)=1;
   end
end
fprintf('End Filtering SoC\n')
fprintf('NNZR(Sf) %10.2f\n',nnz(S)/nn);
fprintf('S density over A:    %10.2f\n',(nnz(S)-nn)/(nnz(A)-nn));

return
