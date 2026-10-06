function [param] = readDefaultParams()
   fileIN = fopen('aspAMG.fnames','r');
   C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_AMG     = D{1};
   C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_SMOOTH  = D{1};
   C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_TSPACE  = D{1};
   C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_COARSEN = D{1};
   C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_PROLONG = D{1};
   C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_FILTER  = D{1};
   fclose(fileIN);
   
   % Read parameters for the AMG hierarchy
   param.amg = read_amg(file_AMG);
   
   % Read parameters for the smoother
   param.smoother = read_smoother(file_SMOOTH);
   param.smoother.nthread = maxNumCompThreads;
   
   % Read parameters for the testspace
   param.tspace = read_tspace(file_TSPACE);
   
   % Read parameters for the smoother
   param.coarsen = read_coarsen(file_COARSEN);
   
   % Read parameters for the prolongation
   param.prolong = read_prolong(file_PROLONG);
   param.smoother.np = maxNumCompThreads;
   
   % Read parameters for the filtering
   param.filter = read_filter(file_FILTER);
   param.smoother.np = maxNumCompThreads;
end
