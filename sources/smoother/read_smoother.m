function param = read_smoother(filename)

ifile = fopen(filename,'r');
C = textscan(fgetl(ifile),'%f'); nthread         = C{1};
C = textscan(fgetl(ifile),'%f'); nupre           = C{1};
C = textscan(fgetl(ifile),'%f'); nupost          = C{1};
C = textscan(fgetl(ifile),'%f'); nstep           = C{1};
C = textscan(fgetl(ifile),'%f'); step_size       = C{1};
C = textscan(fgetl(ifile),'%f'); epsilon         = C{1};
C = textscan(fgetl(ifile),'%s'); method          = C{1}{1};
fclose(ifile);

param.nthread         = nthread;
param.nupre           = nupre;
param.nupost          = nupost;
param.nstep           = nstep;
param.step_size       = step_size;
param.epsilon         = epsilon;
param.method          = method;

return
