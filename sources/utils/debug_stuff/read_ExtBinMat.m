function A = read_ExtBinMat(filename, varargin)

    rect = false;
    if (nargin > 1)
        rect = varargin{1};
    end

    fileID = fopen(filename,'r');
    nrows = fread(fileID,1,'int64','ieee-le');
    if (rect)
        ncols = fread(fileID,1,'int64','ieee-le');
    else
        ncols = nrows;
    end
    nterm = fread(fileID,1,'int64','ieee-le');

    u = reshape(fread(fileID,24*nterm,'uint8','ieee-le'),[24,nterm])';
    fclose(fileID);

    hirow = char(zeros(nterm,0));
    hjcol = char(zeros(nterm,0));
    hval = char(zeros(nterm,0));
    for i = 1 : 4
        hirow = [hirow,dec2hex(u(:,9-i),2)];
        hjcol = [hjcol,dec2hex(u(:,17-i),2)];
        hval = [hval,dec2hex(u(:,25-i),2)];
    end
    for i = 5 : 8
        hirow = [hirow,dec2hex(u(:,9-i),2)];
        hjcol = [hjcol,dec2hex(u(:,17-i),2)];
        hval = [hval,dec2hex(u(:,25-i),2)];
    end

    A = sparse(hex2dec(hirow),hex2dec(hjcol),hex2num(hval),nrows,ncols);

return
