function [mat_I,c_mark] = cpt_ProlMEX_BAMG(level,params,nn_C,fcnodes,mat_S,TV)
%------------------------------------------------------------------------------------------
%
% Computes the prolongation using the BAMG algorithm (
%
% MEX Wrapper function
%
% Input
%
% level:    level number (starting from 1)
% params:   structure for DPLS .np        --> number of omp threads
%                              .itmax_vol --> max iteration for maxVol
%                              .dist_min  --> min distance for interpolation
%                              .dist_max  --> max distance for interpolation
%                              .mmax      --> max number of interpolatory points for
%                                             each F node (NOT USED YET)
%                              .maxcond   --> maximum conditioning bound
%                              .maxrownrm --> max norm in each mat_I row
%                              .tol_vol   --> exit tolerance for max vol iteration
%                              .eps_prol  --> accuracy in the Least Square interpolation
% nn_C:     number of COARSE nodes
% mat_S:    SoC matrix
% fcnodes:  FINE/COARSE partition
% TV:       test vector space
%
% Output
% 
% mat_I:    prolongation matrix
% c_mark:   marker for nodes to be promoted
%
%------------------------------------------------------------------------------------------

% Extract necessary parameters from params
np        = params.np;
itmax_vol = params.itmax_Vol;
dist_min  = params.dist_min;
dist_max  = params.dist_max;
mmax      = size(TV,2);
maxcond   = params.maxcond;
maxrownrm = params.maxrownrm;
tol_vol   = params.tol_Vol;
eps_prol  = params.eps_prol; 

% Unpack the SoC matrix
nn_S = size(mat_S,1);
nt_S = nnz(mat_S);
[iat_S,ja_S,coef_S] = unpack_csr(mat_S);

% Adjust numbering of COARSE nodes in fcnodes
cnt = 0;
for i = 1:nn_S
   if fcnodes(i) > 0
      fcnodes(i) = cnt;
      cnt = cnt + 1;
   else
      fcnodes(i) = -1;
   end
end

% Create some extra parameter
nn_I = nn_S;
nc_I = nn_C;
ntv = size(TV,2);

% Convert reals into integers
level = int32(level);
np = int32(np);
mmax = int32(mmax);
dist_max = int32(dist_max);
nn_S = int32(nn_S);
nt_S = int32(nt_S);
iat_S = int32(iat_S) - 1;
ja_S = int32(ja_S) - 1;
ntv = int32(ntv);
fcnodes = int32(fcnodes);
coef_S = int32(coef_S);
nn_I = int32(nn_I);
nc_I = int32(nc_I);

% Convert TV in proper 1-D array
TV = TV';
TV = TV(:);

[nt_I,iat_I,ja_I,coef_I,c_mark] = ...
                        cpt_Prolongation_BAMG(level,np,itmax_vol,dist_min,dist_max,...
                        mmax,maxcond,maxrownrm,tol_vol,eps_prol,nn_S,nt_S,...
                        iat_S,ja_S,coef_S,ntv,fcnodes,TV,nn_I,nc_I);

% Create a sparse matrix for the prolongation
irow_I = zeros(nt_I,1);
iend = iat_I(1)-1;
for i = 1:nn_I
    istart = iend + 1;
    iend = iat_I(i+1)-1;
    irow_I(istart:iend) = i;
end
nn_I = double(nn_I);
nc_I = double(nc_I);
mat_I = sparse(irow_I,ja_I,coef_I,nn_I,nc_I);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('max irow: %d %d\n',max(irow_I))
%fprintf('size(mat_I): %d %d\n',size(mat_I))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
