function [ascii_input,rhs_build,sym_flag,solv_method,itmax_ruiz,tol_ruiz,itmax,tol,restart]...
         = read_general_RACP(filename)

ifile = fopen(filename,'r');
C = textscan(fgetl(ifile),'%s'); ascii_input = C{1}{1};
C = textscan(fgetl(ifile),'%s'); rhs_build   = C{1}{1};
C = textscan(fgetl(ifile),'%s'); symmetric   = C{1}{1};
C = textscan(fgetl(ifile),'%s'); solv_method = C{1}{1};
C = textscan(fgetl(ifile),'%f'); itmax_ruiz  = C{1};
C = textscan(fgetl(ifile),'%f'); tol_ruiz    = C{1};
C = textscan(fgetl(ifile),'%f'); itmax       = C{1};
C = textscan(fgetl(ifile),'%f'); tol         = C{1};
C = textscan(fgetl(ifile),'%f'); restart     = C{1};
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

switch lower(rhs_build)
    case 'unit_sol'

    case 'unit_rhs'

    case 'rand_sol'

    case 'rand_rhs'

    case 'rhs_in'

    otherwise
       err_msg = ['Wrong value for rhs_build in ' filename];
       error(err_msg);
end

switch lower(symmetric)
    case 'true'
       sym_flag = true;
    case 'false'
       sym_flag = false;
    otherwise
       err_msg = ['Wrong value for symmetric in ' filename];
       error(err_msg);
end

switch lower(solv_method)
    case 'pcg'

    case 'stat_amg'

    case 'bicgstab'

    case 'gmres'

    case 'gmres_cj'

    case 'sqmr'

    otherwise
       err_msg = ['Wrong value for solv_method in ' filename];
       error(err_msg);
end

return
