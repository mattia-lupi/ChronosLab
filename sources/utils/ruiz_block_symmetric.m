function [Abal, D] = ruiz_block_symmetric(A, maxit, tol, verb)
%BALANCE_BLOCK  Iterative diagonal balancing for an NxN block sparse matrix.
%   A    : NxN cell array of sparse blocks (some entries may be empty)
%   D    : Nx1 cell array of column vectors containing diagonal scaling factors
%   Abal : NxN cell array of balanced blocks
if nargin < 2, maxit = 10;   end
if nargin < 3, tol   = 1e-2; end
if nargin < 4, verb  = 1;    end

nb = size(A, 1);

% --- Infer block sizes from first non-empty block in each row ---
blksize = zeros(nb, 1);
for i = 1:nb
   for j = 1:nb
      if ~isempty(A{i,j})
         blksize(i) = size(A{i,j}, 1);
         break;
      end
   end
   if blksize(i) == 0
      error('balance_block: cannot determine size of block-row %d (entire row is empty)', i);
   end
end

% --- Initialize D as cell array of dense column vectors of ones ---
D = cell(nb, 1);
for i = 1:nb
   D{i} = ones(blksize(i), 1);
end

Abal = A;   % COW: no copy until a block is modified

if maxit < 1
   if verb, fprintf('Made 0 iterations\n'); end
   return;
end

for k = 1:maxit
   % 1. Row norms: sum over all non-empty blocks in each block-row
   r = cell(nb, 1);
   for i = 1:nb
      ri = zeros(blksize(i), 1);
      for j = 1:nb
         if ~isempty(Abal{i,j})
            ri = ri + sum(abs(Abal{i,j}).^2, 2);
         end
      end
      r{i} = max(sqrt(ri), 1e-15);
   end

   % 2. Scaling factor calculation
   s = cell(nb, 1);
   for i = 1:nb
      s{i} = 1 ./ sqrt(r{i});
   end

   conv = max(abs(vertcat(s{:}) - 1));
   if verb
      fprintf('ITER: %d | Max(abs(s-1)): %e\n', k, conv);
   end
   if conv < tol, break; end

   % 3. Accumulate D elementwise (O(N) operation)
   for i = 1:nb
      D{i} = D{i} .* s{i};
   end

   % 4. Scale non-empty blocks using implicit expansion (broadcasting)
   for j = 1:nb
      sj_trans = s{j}.';  % Transpose once per block-column
      for i = 1:nb
         if ~isempty(Abal{i,j})
            % (s{i} .* Abal{i,j}) scales rows; .* sj_trans scales columns
            Abal{i,j} = (s{i} .* Abal{i,j}) .* sj_trans;
         end
      end
   end
end

if verb
   fprintf('Made %d iterations\n', k);
   fprintf('Error %e\n', conv);
end
end