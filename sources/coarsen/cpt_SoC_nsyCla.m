function [S,S_patt] = cpt_SoC_nsyCla(tau,A)
%-----------------------------------------------------------------------------------------
%
% Function to compute classical strength of connection.
%
% Input
%
% tau        drop tolerance to remove weak connections
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

% Get the minimum value of each row (which should be negative)
min_row = full(min(A-diag(diag(A))))';

% Treat possible Dirichlet (isolated) nodes
min_row(min_row==0) = 1;

% Compute SoC
[ii,jj,aa] = find(A);
ss = aa ./ min_row(ii);
S = sparse(ii,jj,ss);
fprintf('End Computing SoC\n')

fprintf('Filtering SoC\n')
sf = ss;
sf(sf <= tau) = 0;
S = sparse(ii,jj,sf);
sp = ss;
sp(ss <= tau)=-1;
sp(ss >  tau)=1;
S_patt = sparse(ii,jj,sp);
fprintf('End Filtering SoC\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save('SSS.mat','S','S_patt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Complete the matrix with a unitary diagonal
nn = size(S,1);
S = S + speye(nn);
S_patt = S_patt - diag(diag(S_patt)) + speye(nn);

fprintf('# of S off-diagonal entries:    %10d\n',nnz(S)-nn);
fprintf('S density over A:    %10.2f\n',(nnz(S)-nn)/(nnz(A)-nn));

if DEBINFO.coarsen.draw_dist
   % Draw SoC distribution
   figure; plot(ss(end:-1:1),'-b'); title('SoC distribution');
   figure; plot(sf(end:-1:1),'-r'); title('Filtered SoC distribution');
   drawnow
end

return
