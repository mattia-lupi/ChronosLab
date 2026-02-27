function [fcnodes,ndir] = mark_isolated_nodes(A,fcnodes)

n = size(A,1);
ndir = 0;
for irow = 1:n
    if (nnz(A(:, irow)) == 1)
        ndir = ndir + 1;
        % Mark it as Dirichlet node
        fcnodes(irow) = -2;
    end
end

return
