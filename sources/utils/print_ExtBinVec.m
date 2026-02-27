function print_ExtBinVec(filename,A);

    % Get number of rows and entries
    nrows = size(A,1);
    ncols = size(A,2);
    val = A';
    val = val(:);

    urval = reshape(typecast(double(val),'uint8'),[8,nrows*ncols]);

    fileID = fopen(filename,'w');
    fwrite(fileID,typecast(int64(nrows),'uint8'),'uint8','ieee-be');
    fwrite(fileID,typecast(int64(ncols),'uint8'),'uint8','ieee-be');
    fwrite(fileID,urval,'uint8','ieee-be');
    fclose(fileID);

return
