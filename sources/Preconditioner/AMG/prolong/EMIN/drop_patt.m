function patt_out = drop_patt(nnzr,patt_in,patt_0)

[nrows,ncols] = size(patt_in);

mask_0 = spones(patt_0);
row_counts_0 = full(sum(mask_0, 2));
budget = max(0, nnzr - row_counts_0);

extra = patt_in - patt_in .* mask_0;
[jj_e, ii_e, vv_e] = find(extra.');
if ~isempty(ii_e)
   [~, ord] = sortrows([ii_e, -abs(vv_e)]);
   ii_s = ii_e(ord);
   jj_s = jj_e(ord);
   vv_s = vv_e(ord);
   row_start = find([true; diff(ii_s) ~= 0]);
   row_len   = diff([row_start; numel(ii_s) + 1]);
   ranks     = (1:numel(ii_s)).' - repelem(row_start - 1, row_len);
   keep      = ranks <= budget(ii_s);
   extra_kept = sparse(ii_s(keep), jj_s(keep), vv_s(keep), nrows, ncols);
else
   extra_kept = sparse(nrows, ncols);
end

patt_out = mask_0 + extra_kept;

end
