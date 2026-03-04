function param = read_amg(filename)

ifile = fopen(filename,'r');
C  = textscan(fgetl(ifile),'%f'); nLevMax     = C{1};
C  = textscan(fgetl(ifile),'%f'); maxCoarseSZ = C{1};
fclose(ifile);

param.nLevMax = nLevMax;
param.maxCoarseSZ = maxCoarseSZ;

return
