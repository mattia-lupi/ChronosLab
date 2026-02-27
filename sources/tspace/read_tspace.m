function param = read_tspace(filename)

ifile = fopen(filename,'r');
C = textscan(fgetl(ifile),'%f'); ntv         = C{1};
C = textscan(fgetl(ifile),'%f'); itmax       = C{1};
C = textscan(fgetl(ifile),'%f'); tol         = C{1};
C = textscan(fgetl(ifile),'%s'); method      = C{1}{1};
C = textscan(fgetl(ifile),'%f'); init_approx = C{1};
fclose(ifile);

% Check method and init_approx
switch lower(method)
   case 'none'
   case 'smoothing'
   case 'srqcg'
   case 'lanczos'
   case 'ng-srqcg'
   case 'ng-lanczos'
   case 'arnoldi'
   otherwise
     error('Not existing method');
end

param.ntv = ntv;
param.itmax = itmax;
param.tol = tol;
param.method = method;
param.init_approx = init_approx;

return
