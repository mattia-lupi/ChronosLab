function [iat, ja, coef] = unpack_csr(matrix)
   % 1. Transpose once to convert row iteration into CSC column order
   Mt = matrix.';
   
   % 2. Suppress row index output (~) to avoid allocating an NNZ-length array
   [ja, ~, coef] = find(Mt);
   
   % 3. Count non-zeros per row directly from Mt's internal structure in O(M) time
   row_counts = full(sum(Mt ~= 0, 1)).';
   
   % 4. Pre-pend 1 and compute cumulative sum in a single operation
   iat = cumsum([1; row_counts]);
end