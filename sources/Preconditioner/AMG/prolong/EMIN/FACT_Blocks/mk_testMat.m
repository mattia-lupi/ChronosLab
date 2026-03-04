clear

load BLOCKS_591;
%load BLOCKS_35199;

cat_size = 50;

dim_blk = pt_ind(2:end) - pt_ind(1:end-1);
n_blk = numel(dim_blk);
n_max = max(dim_blk);
n_min = min(dim_blk);

ncat = ceil((n_max-n_min) / cat_size) + 1;
ind_cat = zeros(n_blk,1);

% Assign each block to the proper category
ind_cat = floor( (dim_blk-1) / cat_size ) + 1;

% Create the permutation to sort blocks in increasing category
[~,perm] = sort(ind_cat);

% Copy the set of blocks in a sorted matrix BB
nt = nnz(AA);
ii_B = zeros(nt,1);
jj_B = zeros(nt,1);
aa_B = zeros(nt,1);
pt_cat = ones(ncat+1,1);
pt_BB = ones(n_blk+1,1);
k_cat = ind_cat(perm(1));
% Create pointers
for i = 1:n_blk
   if ~mod(i,100)
      fprintf('Block %d of %d\n',i,n_blk);
   end
   iblk = perm(i);
   istart_A = pt_ind(iblk);
   iend_A = pt_ind(iblk+1)-1;
   dim = iend_A - istart_A + 1;
   pt_BB(i+1) = pt_BB(i) + dim;
   if ind_cat(iblk) > k_cat
      k_cat = k_cat + 1;
      pt_cat(k_cat) = i+1;
   end
end
pt_cat(ncat+1) = n_blk + 1;

% Copy entries in BB
ind_B = 1;
offset = 0;
for i = 1:n_blk
   if ~mod(i,100)
      fprintf('Block %d of %d\n',i,n_blk);
   end
   iblk = perm(i);
   istart_A = pt_ind(iblk);
   iend_A = pt_ind(iblk+1)-1;
   [ii_A,jj_A,aa_A] = find(AA(istart_A:iend_A,istart_A:iend_A));
   dim = numel(ii_A);
   iend_B = ind_B + dim - 1;
   ii_B(ind_B:iend_B) = ii_A + offset;
   jj_B(ind_B:iend_B) = jj_A + offset;
   aa_B(ind_B:iend_B) = aa_A;
   ind_B = ind_B + dim;
   offset = offset + iend_A - istart_A + 1;
end
BB = sparse(ii_B,jj_B,aa_B,size(AA,1),size(AA,2));

% Dump on file
ofile = fopen('TestMats_Categories','w');
fprintf(ofile,'%d    ! cat_size\n',cat_size);
fprintf(ofile,'%d    ! ncat\n',ncat);
fprintf(ofile,'pt_cat:\n');
fprintf(ofile,' %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d\n',pt_cat);
if mod(numel(pt_cat),10) > 0
   fprintf(ofile,'\n');
end
fprintf(ofile,'pt_blk:\n');
fprintf(ofile,' %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d\n',pt_BB);
if mod(numel(pt_BB),10) > 0
   fprintf(ofile,'\n');
end
fclose(ofile);

print_SpMat('BLK_matrix.csr',BB);
