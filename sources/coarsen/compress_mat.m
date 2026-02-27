function A = compress_mat(A, blkSize)

    [iA, jA, cA] = find(A);
    iA = floor((iA-1)/blkSize) + 1;
    jA = floor((jA-1)/blkSize) + 1;
    A = sparse(iA, jA, cA.^2);
    A = sqrt(A);

end
