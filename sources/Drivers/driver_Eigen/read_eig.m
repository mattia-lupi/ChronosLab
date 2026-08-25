function [ascii_input,precond,method,neig,reslambda_check,itmax,tol,largest_flag] = read_eig(filename)

ifile = fopen(filename,'r');
C = textscan(fgetl(ifile),'%s'); ascii_input     = C{1}{1};
C = textscan(fgetl(ifile),'%s'); precond         = C{1}{1};
C = textscan(fgetl(ifile),'%s'); method          = C{1}{1};
C = textscan(fgetl(ifile),'%f'); neig            = C{1};
C = textscan(fgetl(ifile), '%s'); reslambda_check = strcmpi(C{1}{1}, 'true');
C = textscan(fgetl(ifile),'%f'); itmax           = C{1};
C = textscan(fgetl(ifile),'%f'); tol             = C{1};
C = textscan(fgetl(ifile),'%s'); largest_flag    = C{1}{1};
fclose(ifile);

switch lower(ascii_input)
    case 'true'
       ascii_input = true;
    case 'false'
       ascii_input = false;
    otherwise
       err_msg = ['Wrong value for ascii_input in ' filename];
       error(err_msg);
end

switch lower(precond)
    case 'fsai'

    case 'amg'

    otherwise
       err_msg = ['Wrong value for precond in ' filename];
       error(err_msg);
end

switch lower(method)
    case 'srqcg'

    case 'defl_srqcg'

    case 'lanczos'

    case 'lobpcg'

   case 'block_arnoldi' 

    otherwise
       err_msg = ['Wrong value for method in ' filename];
       error(err_msg);
end

return
