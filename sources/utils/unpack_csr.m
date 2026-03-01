function [iat,ja,coef] = unpack_csr(matrix)

   % 1. Transpose to switch from CSC (MATLAB default) to CSR ordering
   % find() on the transpose iterates over rows of the original matrix sequentially.
   [ja, irow, coef] = find(matrix.');

   % 2. Vectorized construction of row pointers (iat)
   nn = size(matrix, 1);

   % Count non-zeros per row. accumarray handles empty rows automatically (returns 0).
   row_counts = accumarray(irow, 1, [nn, 1]);

   % Cumulative sum to generate pointers.
   % Prepend 1 because CSR pointers are 1-based in this context.
   iat = [1; cumsum(row_counts) + 1];

end
