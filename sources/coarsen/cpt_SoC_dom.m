function [S,S_patt] = cpt_SoC_cla(tau,A,verb)
%-----------------------------------------------------------------------------------------
%
% Function to compute strength of connection based on diagonal dominance
%
% Input
%
% tau        >= 0 --->  drop tolerance to remove weak connections
%            <  0 --->  avg number of connections per row to keep in filtered SoC
% A          system matrix
%
% Output
%
% S          SoC matrix with 1 indicating strong connections and zero elsewhere
% S_patt     SoC matrix with 1 indicating strong connections and zero elsewhere
%
%-----------------------------------------------------------------------------------------

global DEBINFO

if verb
   fprintf('Computing SoC\n')
end

% Get the minimum value of each column
dd = full(diag(A));
dd = dd(:);

% Compute SoC
[ii,jj,aa] = find(A);
ss = abs(aa) ./ sqrt(dd(ii).*dd(jj));

if verb
   fprintf('End Computing SoC\n')

   fprintf('Filtering SoC\n')
   fprintf('Average NNZR before filtering: %f\n',nnz(A)/size(A,1))
end
if tau >= 0
   sf = ss;
   sf(sf <= tau)=0;
   S = sparse(ii,jj,sf);
   sp = ss;
   sp(ss <= tau)=-1;
   sp(ss >  tau)=1;
   S_patt = sparse(ii,jj,sp);
else
   nn = size(A,1);
   ncut = ceil(-nn*tau);
   if ncut<nnz(A)
       [ss,perm] = sort(ss,'descend');
       ii_f = ii(perm(1:ncut));
       jj_f = jj(perm(1:ncut));
       sf = ss(1:ncut);
       S = sparse(ii_f,jj_f,sf);
       ss = 1;
       ss(ncut+1:end) = -1;
       S_patt = sparse(ii,jj,ss);
   else
       S = sparse(ii,jj,ss);
       S_patt = S;
       S_patt(S_patt>0)=1;
   end
end
if verb
   fprintf('Average NNZR after filtering: %f\n',nnz(S)/size(S,1))
   fprintf('End Filtering SoC\n')
%save('SSS','S','S_patt');
end

nn = size(S,1);
if verb
   fprintf('# of S off-diagonal entries:    %10d\n',nnz(S)-nn);
   fprintf('S density over A:    %10.2f\n',(nnz(S)-nn)/(nnz(A)-nn));
end

%if DEBINFO.coarsen.draw_dist
%   % Draw SoC distribution
%   figure; plot(ss(end:-1:1),'-b'); title('SoC distribution');
%   figure; plot(sf(end:-1:1),'-r'); title('Filtered SoC distribution');
%   drawnow
%end

return
