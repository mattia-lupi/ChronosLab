function param = read_filter(filename)
  
ifile = fopen(filename,'r');
C  = textscan(fgetl(ifile),'%f'); np         = C{1};
C  = textscan(fgetl(ifile),'%f'); filt_wgt   = C{1};
C  = textscan(fgetl(ifile),'%f'); filt_tol   = C{1};
C  = textscan(fgetl(ifile),'%f'); filt_tau   = C{1};
C  = textscan(fgetl(ifile),'%f'); min_patt   = C{1};
fclose(ifile);

param.np        = np;
param.filt_wgt  = filt_wgt;
param.filt_tol  = filt_tol;
param.filt_tau  = filt_tau;
param.min_patt  = min_patt;

return
