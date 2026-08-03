function [jat, ia, coef] = unpack_csc(matrix)
   % 1. Suppress column index output (~) to avoid allocating an NNZ-length array
   [ia, ~, coef] = find(matrix);
   
   % 2. Get non-zeros per column directly from sparse structure in O(N) instead of O(NNZ)
   col_counts = full(sum(matrix ~= 0, 1)).';
   
   % 3. Pre-pend 1 and compute cumulative sum in a single operation
   jat = cumsum([1; col_counts]);
end