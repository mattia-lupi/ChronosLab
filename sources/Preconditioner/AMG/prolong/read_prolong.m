function param = read_prolong(filename)

ifile = fopen(filename,'r');
C  = textscan(fgetl(ifile),'%s'); D  = C{1}; proltype = D{1};
C  = textscan(fgetl(ifile),'%s'); D  = C{1}; prol_emin = D{1};
C  = textscan(fgetl(ifile),'%f'); np            = C{1};
C  = textscan(fgetl(ifile),'%f'); itmax_Vol     = C{1};
C  = textscan(fgetl(ifile),'%f'); toll_Vol      = C{1};
C  = textscan(fgetl(ifile),'%f'); dist_min      = C{1};
C  = textscan(fgetl(ifile),'%f'); dist_max      = C{1};
C  = textscan(fgetl(ifile),'%f'); maxcond       = C{1};
C  = textscan(fgetl(ifile),'%f'); maxrownrm     = C{1};
C  = textscan(fgetl(ifile),'%f'); eps_prol      = C{1};
C  = textscan(fgetl(ifile),'%f'); updateCF      = C{1};
C  = textscan(fgetl(ifile),'%f'); patt_pow      = C{1};
C  = textscan(fgetl(ifile),'%f'); patt_tau      = C{1};
C  = textscan(fgetl(ifile),'%f'); nnzr_max      = C{1};
C  = textscan(fgetl(ifile),'%f'); itmax_emin    = C{1};
C  = textscan(fgetl(ifile),'%f'); energ_tol     = C{1};
C  = textscan(fgetl(ifile),'%f'); condmax_emin  = C{1};
C  = textscan(fgetl(ifile),'%f'); maxwgt_emin   = C{1};
C  = textscan(fgetl(ifile),'%f'); prec_emin     = C{1};
C  = textscan(fgetl(ifile),'%f'); solv_emin     = C{1};
C  = textscan(fgetl(ifile),'%f'); min_lfil      = C{1};
C  = textscan(fgetl(ifile),'%f'); max_lfil      = C{1};
C  = textscan(fgetl(ifile),'%f'); D_lfil        = C{1};

fclose(ifile);

% Check that flags are ok
if ( ~( strcmpi(proltype,'BAMG') || strcmpi(proltype,'EXTI')||...
        strcmpi(proltype,'CLAS') || strcmpi(proltype,'HYBC') ) )
   error('Wrong proltype flag: %s\n',proltype);
end
if ( ~( strcmpi(prol_emin,'EMIN') || strcmpi(prol_emin,'SMOOTH') || ...
        strcmpi(prol_emin,'NONE') ) )
   error('Wrong proltype flag: %s\n',prol_emin);
end

% Check EMIN parameters
if strcmpi(prol_emin,'EMIN')

   % Check preconditioner choice is correct
   if (prec_emin < 1 || prec_emin > 3)
      error('Wrong prec_emin: %d\n',prol_emin);
   end

   % Check EMIN parameters are ok
   if (solv_emin >= 1 && solv_emin <= 3)
      % SGS is not supported
      if (prec_emin == 3)
         error('SGS (prol_emin: %d) is not supported for solv_emin: %d\n',...
                prec_emin,solv_emin);
      end
   elseif (solv_emin < 1 || solv_emin > 4)
      error('Wrong solv_emin: %d\n',solv_emin);
   end
   if (solv_emin ~= 4)
      fprintf('WARNING: recommended EMIN solver is solv_emin: %d\n',solv_emin);
      fprintf('         the other solvers are experimental');
   end
end

param.proltype      = proltype;
param.prol_emin     = prol_emin;
param.np            = np;
param.itmax_Vol     = itmax_Vol;
param.tol_Vol       = toll_Vol;
param.dist_min      = dist_min;
param.dist_max      = dist_max;
param.maxcond       = maxcond;
param.maxrownrm     = maxrownrm;
param.eps_prol      = eps_prol;
param.updateCF      = updateCF;
param.patt_pow      = patt_pow;
param.patt_tau      = patt_tau;
param.nnzr_max      = nnzr_max;
param.itmax_emin    = itmax_emin;
param.energ_tol     = energ_tol;
param.condmax_emin  = condmax_emin;
param.maxwgt_emin   = maxwgt_emin;
param.prec_emin     = prec_emin;
param.solv_emin     = solv_emin;
param.min_lfil      = min_lfil;
param.max_lfil      = max_lfil;
param.D_lfil        = D_lfil;

return
