function print_ExtBinMat(filename,A);

   if (size(A,1) == size(A,2))
      % Fix pattern to be symmetric
      At = A';
      mask = spones(At) > spones(A);
      mismatches = nnz(mask);

      if mismatches
         fprintf('Found %d pattern entries mismatch, padding...\n', mismatches*2);
         A = A + (A' .* mask) * eps;
         mask = spones(A') - spones(A);
         if nnz(mask)
            error();
         else
            fprintf('Fixed pattern\n');
         end
      end
   end

   % Get irow, jcol and coef of transposed matrix (Matlab uses CSC, we prefer CSR)
   [irow,jcol,val] = find(A');

   % Get number of rows and entries
   nrows = size(A,1);
   nterm = nnz(A);

   if (nrows == size(A,2))
      % Check if there are null diagonal terms and count them
      full_set = [1:nrows]';
      null_diag = setdiff(full_set,irow(irow==jcol));
      nnull = length(null_diag);

      if (nnull > 0)
          % Update the number of total entries
          nterm = nterm + nnull;

          % Update irow, jcol and coef
          irow = [irow;null_diag];
          jcol = [jcol;null_diag];
          val = [val;zeros(nnull,1)];
      end

      if (nnull > 0)
          % Save all arrays in a matrix
          P = [jcol,irow,val];

          % Sort new rows and columns
          P = sortrows(P,[1,2]);

          jcol = P(:,1);
          irow = P(:,2);
          val = P(:,3);
      end
   end

   urjcol = reshape(typecast(int64(jcol),'uint8'),[8,nterm]);
   urirow = reshape(typecast(int64(irow),'uint8'),[8,nterm]);
   urval = reshape(typecast(double(val),'uint8'),[8,nterm]);
   u = [urjcol;urirow;urval];

   fileID = fopen(filename,'w');
   fwrite(fileID,typecast(int64(nrows),'uint8'),'uint8','ieee-be');
   if (nrows ~= size(A,2))
       fwrite(fileID,typecast(int64(size(A,2)),'uint8'),'uint8','ieee-be');
   end
   fwrite(fileID,typecast(int64(nterm),'uint8'),'uint8','ieee-be');
   fwrite(fileID,u,'uint8','ieee-be');
   fclose(fileID);

return
