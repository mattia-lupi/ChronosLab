%-----------------------------------------------------------------------------------------
% Simple program to generate a Porous Flow Matrix starting from a grid of tetrahedrons
% and a list of 3D permeabilities.
%-----------------------------------------------------------------------------------------

clear
clc
close all

%file_coord = 'MESH/coord.txt';
%file_tetra = 'MESH/tetra.txt';
file_coord = 'MESH/Cubo_591.coor';
file_tetra = 'MESH/Cubo_591.tetra';
%file_coord = 'MESH/Cubo_4820.coor';
%file_tetra = 'MESH/Cubo_4820.tetra';
%file_coord = 'MESH/Cubo_35199.coor';
%file_tetra = 'MESH/Cubo_35199.tetra';

%file_k_perm = 'MESH/perm_aniso.txt';
file_k_perm = 'MESH/perm_iso.txt';

% Load coordinates
coord = dlmread(file_coord,'',1,0);
coord = coord(:,2:4);

% Load tetrahedrons
tetra  = dlmread(file_tetra,'',1,0);
material = tetra(:,end);
tetra = tetra(:,2:5);

% Load permeabilities (kxx, kyy, kzz, kxy, kxz, kyz)
k_perm = dlmread(file_k_perm,'',0,0);

% Create the list of permeability tensors
K_list = mk_perm_tensor(k_perm);

% Assemble the stiffness matrix
tic;
nnod = size(coord,1);
nele = size(tetra,1);
nmax = 16*nele;
irow = zeros(nmax,1);
jcol = zeros(nmax,1);
coef = zeros(nmax,1);
ii = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];
jj = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
ind = 1;
for iel = 1:size(tetra,1)
   % Retrieve permeability tensor
   K = K_list{material(iel)};
   % Compute the local matrix
   Aloc = cpt_loc_ANISOMAT(tetra(iel,:),coord,K);
   % Add the entries of the local matrix to the list for assembly
   irow(ind:ind+15) = tetra(iel,ii);
   jcol(ind:ind+15) = tetra(iel,jj);
   coef(ind:ind+15) = Aloc(:);
   ind = ind + 16;
end
% Create the sparse matrix
A = sparse(irow,jcol,coef,nnod,nnod);
% Symmtrize the matrix
A = 0.5*(A+A');
toc

print_SpMat('Mat.out',A);

return
