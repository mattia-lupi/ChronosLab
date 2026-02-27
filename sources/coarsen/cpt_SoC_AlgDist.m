function [S,S_patt] = cpt_SoC_AlgDist(tau,smootherOp,A)
%-----------------------------------------------------------------------------------------
%
% Function to compute strength of connection based on algebraic distance
%
% Input
%
% tau        >= 0 --->  drop tolerance to remove weak connections
%            <  0 --->  avg number of connections per row to keep in filtered SoC
% smootherOp smoother used to unveil strength of connection
% A          system matrix
%
% Output
%
% S          SoC matrix with 1 indicating strong connections and zero elsewhere
% S_patt     SoC matrix with 1 indicating strong connections and zero elsewhere
%
%-----------------------------------------------------------------------------------------

global DEBINFO

fprintf('Computing SoC\n')

% Generare a random vector
nn = size(A,1);
x = rand(nn,1);

% Perform k steps of relaxation
k = 50;
Minv = @(x) smootherOp.omega*(smootherOp.right*(smootherOp.left*x));
nrm_0 = norm(x);
nrm_old = nrm_0;
fprintf('Computing Algebraic Distance SoC\n');
fprintf('%4s   %12s\n','iter','res_norm');
for i = 1:k
   x = x - Minv(A*x);
   fprintf('%3d   %12.3e %12.3e\n',i,norm(x)/nrm_0,norm(x)/nrm_old);
   nrm_old = norm(x);
end

% Compute SoC
[ii,jj,ss] = find(A);
ss = 1 - abs( x(ii) - x(jj) ) ./ max( abs(x(ii)),abs(x(jj)) );

fprintf('End Computing SoC\n')

fprintf('Filtering SoC\n')
fprintf('Average NNZR before filtering: %f\n',nnz(A)/size(A,1))
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
   [ss,perm] = sort(ss,'descend');
   ii_f = ii(perm(1:ncut));
   jj_f = jj(perm(1:ncut));
   sf = ss(1:ncut);
   S = sparse(ii_f,jj_f,sf);
   ss = 1;
   ss(ncut+1:end) = -1;
   S_patt = sparse(ii,jj,ss);
end
fprintf('Average NNZR after filtering: %f\n',nnz(S)/size(S,1))
fprintf('End Filtering SoC\n')
%save('SSS','S','S_patt');

nn = size(S,1);
fprintf('# of S off-diagonal entries:    %10d\n',nnz(S)-nn);
fprintf('S density over A:    %10.2f\n',(nnz(S)-nn)/(nnz(A)-nn));

%if DEBINFO.coarsen.draw_dist
%   % Draw SoC distribution
%   figure; plot(ss(end:-1:1),'-b'); title('SoC distribution');
%   figure; plot(sf(end:-1:1),'-r'); title('Filtered SoC distribution');
%   drawnow
%end

return
