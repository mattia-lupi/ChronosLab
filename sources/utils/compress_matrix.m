function B = compress_matrix(A,stride)
   [ii,jj,aa] = find(A);
   aa(:) = 1;
   ii = ceil(ii/stride);
   jj = ceil(jj/stride);
   B = sparse(ii,jj,aa);
end
