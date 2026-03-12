function vec = read_ExtBinVec(filename)

   fid = fopen(filename, 'rb');
   if fid < 0
     warning('Cannot open file: %s', fname);
   end
   
   % Read header
   N = fread(fid, 1, 'int64');
   M = fread(fid, 1, 'int64');
   
   % Sanity checks
   if M ~= 1
     warning('File %s expected M=1 but got M=%d. Continuing...', fname, M);
   end
   
   % Read vector
   vec = fread(fid, N*M, 'double');
   
   fclose(fid);
end
