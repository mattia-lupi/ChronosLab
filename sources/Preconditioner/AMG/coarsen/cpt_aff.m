function affinity = cpt_aff(v,w)

nrm_v = norm(v);
nrm_w = norm(w);

if nrm_v < 1.e-10 || nrm_w < 1.e-10
   % Set affinity to zero if one of the two nodes has zero norm
   affinity = 0;
else
   affinity = abs(v'*w)/(nrm_v*nrm_w);
end

return
