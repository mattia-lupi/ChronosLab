function print_SpMat(filename,A)

fileID = fopen(filename,'w');

[i,j,aa] = find(A);

[i,perm] = sort(i);
j = j(perm);
aa = aa(perm);

P = [i,j,aa]';
fprintf(fileID,'%12d %12d\n',size(A,1),nnz(A));
fprintf(fileID,'%10d %10d %24.15e\n',P);

fclose(fileID);
