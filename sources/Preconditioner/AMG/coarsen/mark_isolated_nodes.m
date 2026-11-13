function [fcnodes,ndir] = mark_isolated_nodes(A,fcnodes)

col_counts = full(sum(A ~= 0, 1))';
iso = find(col_counts == 1);
fcnodes(iso) = -2;
ndir = numel(iso);

return
