function [jat,ia,coef] = unpack_csc(matrix)
   % 1. Use find() directly on the CSC matrix (MATLAB default column-major order)
   % This iterates over columns sequentially.
   [ia, jcol, coef] = find(matrix);

   % 2. Vectorized construction of column pointers (jat)
   mm = size(matrix, 2);

   % Count non-zeros per column.
   col_counts = accumarray(jcol, 1, [mm, 1]);

   % Cumulative sum to generate pointers.
   jat = [1; cumsum(col_counts) + 1];
end
