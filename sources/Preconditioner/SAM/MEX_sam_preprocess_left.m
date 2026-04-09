function preproc = MEX_sam_preprocess_left(Ak, S)
n = size(Ak, 1);

if nargin < 2 || isempty(S)
    S = spones(Ak);
else
    S = spones(S);
end

[iatk, jak, ~] = unpack_csr(Ak);
[iats, jas, ~] = unpack_csr(S);

[s_ptr, s_data, r_ptr, r_data, nnz_total] = ...
    sam_preprocess_left_mex(iatk, jak, iats, jas);

preproc.n         = n;
preproc.S         = logical(S);
preproc.s_ptr     = s_ptr;
preproc.s_data    = s_data;
preproc.r_ptr     = r_ptr;
preproc.r_data    = r_data;
preproc.nnz_total = nnz_total;
end
