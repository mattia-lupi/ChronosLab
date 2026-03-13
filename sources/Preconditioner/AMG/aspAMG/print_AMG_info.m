%-----------------------------------------------------------------------------------------
%
% Function to print information about AMG hierarchy
%
%-----------------------------------------------------------------------------------------
function print_AMG_info(info,ofile)

n_lev = numel(info);

if nargin == 1

   % Print Operator information
   fprintf('\n')
   fprintf('-------------------------------------------------------------------\n')
   fprintf('                      OPERATOR INFORMATION\n')
   fprintf('-------------------------------------------------------------------\n')
   fprintf(' %5s %10s %10s %10s %13s %13s\n','level','nn','nc/nn %','nnz','avg nnz/row','max Lambda')
   fprintf(' %5d %10d %10s %10d %13.2f %13.2e\n',1,info(1).nn,'----',info(1).nt,info(1).nt/info(1).nn,info(1).lmax_A);
   for i = 2:n_lev
      fprintf(' %5d %10d %10.2f %10d %13.2f %13.2e\n',i,info(i).nn,100*info(i).nn/info(i-1).nn,info(i).nt,info(i).nt/info(i).nn,info(i).lmax_A);
   end
   fprintf('-------------------------------------------------------------------\n')

   % Print Smother information
   fprintf('\n')
   fprintf('-----------------------------------------------------\n')
   fprintf('               SMOOTHER INFORMATION\n')
   fprintf('-----------------------------------------------------\n')
   fprintf(' %5s %10s %10s %13s %10s\n','level','nn','nnz','avg nnz/row','omega')
   for i = 1:n_lev-1
      fprintf(' %5d %10d %10d %13.2f %10.4f\n',i,info(i).nn,info(i).nt_S,info(i).nt_S/info(i).nn,info(i).omega);
   end
   fprintf(' %5d %10d %10s %13s %10s\n',n_lev,info(n_lev).nn,'----','----','----');
   fprintf('-----------------------------------------------------\n')

   % Print Prolongation information
   fprintf('\n')
   fprintf('------------------------------------------\n')
   fprintf('       PROLONGATION INFORMATION\n')
   fprintf('------------------------------------------\n')
   fprintf(' %5s %10s %10s %13s\n','level','nn','nnz','avg nnz/row')
   for i = 1:n_lev-1
      fprintf(' %5d %10d %10d %13.2f\n',i,info(i).nn,info(i).nt_P,info(i).nt_P/info(i).nn);
   end
   fprintf(' %5d %10d %10s %13s\n',n_lev,info(n_lev).nn,'----','----');
   fprintf('------------------------------------------\n')

else

   % Print Operator information
   fprintf(ofile,'\n');
   fprintf(ofile,'-------------------------------------------------------------------\n');
   fprintf(ofile,'                      OPERATOR INFORMATION\n');
   fprintf(ofile,'-------------------------------------------------------------------\n');
   fprintf(ofile,' %5s %10s %10s %10s %13s %13s\n','level','nn','nc/nn %','nnz','avg nnz/row','max Lambda');
   fprintf(ofile,' %5d %10d %10s %10d %13.2f %13.2e\n',1,info(1).nn,'----',info(1).nt,info(1).nt/info(1).nn,info(1).lmax_A);
   for i = 2:n_lev
      fprintf(ofile,' %5d %10d %10.2f %10d %13.2f %13.2e\n',i,info(i).nn,100*info(i).nn/info(i-1).nn,info(i).nt,info(i).nt/info(i).nn,info(i).lmax_A);
   end
   fprintf(ofile,'-------------------------------------------------------------------\n');


   % Print Smother information
   fprintf(ofile,'\n');
   fprintf(ofile,'-----------------------------------------------------\n');
   fprintf(ofile,'               SMOOTHER INFORMATION\n');
   fprintf(ofile,'-----------------------------------------------------\n');
   fprintf(ofile,' %5s %10s %10s %13s %10s\n','level','nn','nnz','avg nnz/row','omega');
   for i = 1:n_lev-1
      fprintf(ofile,' %5d %10d %10d %13.2f %10.4f\n',i,info(i).nn,info(i).nt_S,info(i).nt_S/info(i).nn,info(i).omega);;
   end
   fprintf(ofile,' %5d %10d %10s %13s %10s\n',n_lev,info(n_lev).nn,'----','----','----');;
   fprintf(ofile,'-----------------------------------------------------\n');

   % Print Prolongation information
   fprintf(ofile,'\n');
   fprintf(ofile,'------------------------------------------\n');
   fprintf(ofile,'       PROLONGATION INFORMATION\n');
   fprintf(ofile,'------------------------------------------\n');
   fprintf(ofile,' %5s %10s %10s %13s\n','level','nn','nnz','avg nnz/row');
   for i = 1:n_lev-1
      fprintf(ofile,' %5d %10d %10d %13.2f\n',i,info(i).nn,info(i).nt_P,info(i).nt_P/info(i).nn);
   end
   fprintf(ofile,' %5d %10d %10s %13s\n',n_lev,info(n_lev).nn,'----','----');
   fprintf(ofile,'------------------------------------------\n');

end

end
