function preproc = MEX_sam_preprocess_left(Ak, S)
n = size(Ak, 1);

if nargin < 2 || isempty(S)
    S = spones(Ak);
else
    % Force conversion to double sparse matrix to interface seamlessly with C++ TypedArray<double>
    S = spones(S); 
end

% Execute C++ MEX 
[iatk,jak,coefk] = unpack_csr(Ak);
[s_idx, r_idx, nnz_total] = sam_preprocess_left_mex(iatk,jak,coefk, S);

% Pack struct
preproc.n         = n;
preproc.S         = logical(S);
preproc.s_idx     = s_idx;
preproc.r_idx     = r_idx;
preproc.nnz_total = nnz_total;
end
