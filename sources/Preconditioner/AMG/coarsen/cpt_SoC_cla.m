function [S,S_patt] = cpt_SoC_cla(tau,A,verb)
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

if verb
   fprintf('Computing SoC\n')
end

[nRows, nCols] = size(A);

% Get off-diagonal part of A
A_off = A - diag(diag(A));

% Get the minimum value of each column
min_col = full(min(A_off))';

% Check if all column minimums are 0 (no negative off-diagonal entries)
if all(min_col >= 0)
   if verb
      fprintf('All off-diagonal column minimums >= 0: falling back to absolute value SoC.\n');
   end

   % Fallback: compute SoC based on absolute value magnitudes
   % Maximum off-diagonal magnitude per column
   max_abs_col = full(max(abs(A_off)))';
   max_abs_col(max_abs_col == 0) = 1; % Handle isolated/Dirichlet nodes

   [ii,jj,aa] = find(A);
   ss = abs(aa) ./ max(max_abs_col(ii), max_abs_col(jj));
else
   % Standard Classical M-matrix SoC:
   % Treat possible Dirichlet (isolated) nodes
   min_col(min_col==0) = 1;

   % Compute SoC
   [ii,jj,aa] = find(A);
   ss = aa ./ max(min_col(ii),min_col(jj));
end

if verb
   fprintf('End Computing SoC\n')

   fprintf('Filtering SoC\n')
end

sf = ss;
sf(sf <= tau) = 0;
S = sparse(ii,jj,sf,nRows,nCols);
sp = ss;
sp(ss <= tau)=-1;
sp(ss >  tau)=1;
S_patt = sparse(ii,jj,sp,nRows,nCols);
if verb
   fprintf('Average NNZR after filtering: %f\n',nnz(S)/size(S,1))
   fprintf('End Filtering SoC\n')
end

% Complete the matrix with a unitary diagonal
nn = size(S,1);
S = S + speye(nn);
S_patt = S_patt - diag(diag(S_patt)) + speye(nn);

if verb
   fprintf('# of S off-diagonal entries:    %10d\n',nnz(S)-nn);
   fprintf('S density over A:    %10.2f\n',(nnz(S)-nn)/(nnz(A)-nn));
end

if DEBINFO.coarsen.draw_dist
   % Draw SoC distribution
   figure; plot(ss(end:-1:1),'-b'); title('SoC distribution');
   figure; plot(sf(end:-1:1),'-r'); title('Filtered SoC distribution');
   drawnow
end

return
