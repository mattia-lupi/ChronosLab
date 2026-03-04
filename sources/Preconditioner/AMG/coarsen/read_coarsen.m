function param = read_coarsen(filename)

ifile = fopen(filename,'r');
C  = textscan(fgetl(ifile),'%s'); D  = C{1}; SoCtype  = D{1};
C  = textscan(fgetl(ifile),'%f'); tau        = C{1};
C  = textscan(fgetl(ifile),'%f'); nl_agg     = C{1};
fclose(ifile);

param.SoC_type = SoCtype;
param.tau      = tau;
param.nl_agg   = nl_agg;

return
