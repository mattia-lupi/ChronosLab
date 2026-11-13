function [mat_I,c_mark] = cpt_ProlMEX_BAMG(level,params,nn_C,fcnodes,mat_S,TV,verb)
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
is_c = fcnodes > 0;
fcnodes(is_c) = int32(0:sum(is_c)-1);
fcnodes(~is_c) = int32(-1);

% Create some extra parameter
nn_I = nn_S;
nc_I = nn_C;
ntv = size(TV,2);

% Convert reals into integers
iat_S = int32(iat_S) - 1;
ja_S = int32(ja_S) - 1;
fcnodes = int32(fcnodes);
coef_S = int32(coef_S);
verb = double(verb);

% Convert TV in proper 1-D array
TV = TV';
TV = TV(:);

[nt_I,iat_I,ja_I,coef_I,c_mark] = ...
                        cpt_Prolongation_BAMG(level,np,itmax_vol,dist_min,dist_max,...
                        mmax,maxcond,maxrownrm,tol_vol,eps_prol,nn_S,nt_S,...
                        iat_S,ja_S,coef_S,ntv,fcnodes,TV,nn_I,nc_I,verb);

% Create a sparse matrix for the prolongation
irow_I = repelem((1:nn_I)', diff(iat_I));
nn_I = double(nn_I);
nc_I = double(nc_I);
mat_I = sparse(irow_I,ja_I,coef_I,nn_I,nc_I);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('max irow: %d %d\n',max(irow_I))
%fprintf('size(mat_I): %d %d\n',size(mat_I))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
