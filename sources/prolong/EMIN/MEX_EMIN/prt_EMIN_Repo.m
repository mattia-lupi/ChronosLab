function prt_EMIN_Repo(AMG_hrc)

% Open report file
rf = fopen('EMIN_Report.txt','w');

% Print header
fprintf(rf,'%5s %5s %10s %10s %8s','Level','niter','relres','avg_lfil','nnzr');
fprintf(rf,'  %10s  %10s','Init_Ener','Fin_Ener');
fprintf(rf,'  %10s  %10s  %10s','nnz_K','nnz_PK','nnz_ZQ');
fprintf(rf,' %8s %8s %8s %8s %8s','T_gathK','T_precK','T_gathB','T_PCG','T_ovh');
fprintf(rf,' %8s %8s\n','T_emin','T_prol');

level = 1;
while true
   if isstruct(AMG_hrc.next)

      % General level
      info = AMG_hrc.emin_info;
      fprintf(rf,'%5d',level);
      fprintf(rf,' %5d',info.iter);
      fprintf(rf,' %10.2e',info.relres);
      fprintf(rf,' %10.2e',info.avg_lfil);
      fprintf(rf,' %8.2f',info.nnzr);
      fprintf(rf,' %11.2e',info.Energy_init);
      fprintf(rf,' %11.2e',info.Energy_final);
      fprintf(rf,' %11d',info.nnz_K);
      fprintf(rf,' %11d',info.nnz_PK);
      fprintf(rf,' %11d',info.nnz_ZQ);
      fprintf(rf,' %8.2f',info.time_gath_K);
      fprintf(rf,' %8.2f',info.time_prec_K);
      fprintf(rf,' %8.2f',info.time_gath_B);
      fprintf(rf,' %8.2f',info.time_PCG);
      fprintf(rf,' %8.2f',info.time_overhead);
      fprintf(rf,' %8.2f',info.time_glob);
      fprintf(rf,' %8.2f',info.time_prolTot);
      fprintf(rf,'\n');
   else

      % Last level
      break

   end

   % Prepare for the next level
   level = level + 1;
   AMG_hrc = AMG_hrc.next;

end

   % Close report file
   fclose(rf);

end
