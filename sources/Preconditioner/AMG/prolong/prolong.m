function [P,clist,fcnode,emin_info] = prolong(two_stg_prol,level,symm_flag,param,...
                                              clist,fclist,...
                                              fflist,fcnode,S,S_patt,smootherOp,A,TV,verb)
%-----------------------------------------------------------------------------------------
%
% Main driver for the prolongation set-up
%
% Input
%
% two_stg_prol: compute a two-stage prolongation (because of aggressive coarsening)
% param:        structure containing all the parameters for AMG construction
% clist:        list of coarse nodes
% fclist:       list of fine/coarse nodes (coarse in the first stage)
% fflist:       list of fine nodes (fine from the first stage)
% fcnode:       C/F indicator
% S:            filtered strength of connection matrix
% S_patt:       pattern of strength of connection matrix indicating S/W connections
% smootherOP:   smoothing operator
% A:            system matrix
% TV:           initial testspace which is eventually padded to the desired size
%
% Output
%
% P:            prolongation matrix
% emin_info:    some detailed info on EMIN computation
%
%-----------------------------------------------------------------------------------------
if nargin < 14
   verb = 1;
end

% Compute prolongation
nc = numel(clist);
if ~two_stg_prol
   % Standard prolongation
   switch upper(param.prolong.proltype)
      case 'BAMG'
         [P,c_mark] = cpt_ProlMEX_BAMG(level,param.prolong,nc,fcnode,S,TV,verb);
         if verb 
            fprintf('Number of nodes that are not perfectly interpolated: %d\n',sum(c_mark));
         end
         %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if sum(c_mark) > 0
            if param.prolong.updateCF == 1
               if verb
                  fprintf('Update the coarse nodes only\n');
               end
               % Update coarsening
               ind_new = find(c_mark==1);
               n_add = numel(ind_new);
               [ii,jj,pp] = find(P);
               ii = [ii; ind_new];
               jj = [jj; (nc+1:nc+n_add)'];
               pp = [pp; ones(n_add,1)];
               % Add coatse nodes to prolongation
               P = sparse(ii,jj,pp,numel(fcnode),nc+n_add);
               fcnode(ind_new) = (nc+1:nc+n_add)';
               nc = nc + n_add;
               clist = [clist; ind_new];
            elseif param.prolong.updateCF == 2
               if verb
                  fprintf('Update the coarse nodes and recompute prolongation\n');
               end
               % Update coarsening
               fcnode(fcnode>0) = 1;
               fcnode(fcnode<0) = 0;
               fcnode = fcnode + c_mark;
               nc = sum(fcnode);
               nf = numel(fcnode) - nc;
               fcnode(fcnode>0) = 1:nc;
               fcnode(fcnode==0) = -(1:nf);
               clist = find(fcnode>0);
               % Recompute prolongation
               [P,c_mark] = cpt_ProlMEX_BAMG(level,param.prolong,nc,fcnode,S,TV);
               if verb
                  fprintf('Number of nodes that are not perfectly interpolated: %d\n',sum(c_mark));
               end
            elseif (param.prolong.updateCF ~= 0)
               error('Wrong value of updateCF: %d\n',updateCF);
            end
         end
         %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         v = vecnorm(P');
         [amax,ind] = max(v);
         if verb
            fprintf('Max row norm: %15.6e\n',full(amax));
            fprintf('On row:       %d\n',ind);
         end
      case 'CLAS'
         [P] = cpt_ProlMEX_Classical(param.prolong,nc,fcnode,S_patt,A);
      case 'HYBC'
         [P] = cpt_ProlMEX_HybridClas(param.prolong,nc,fcnode,S_patt,A);
      case 'EXTI'
         fcnode_exti = fcnode;
         fcnode_exti(fcnode_exti<0)= -1;
         % Create a compressed A with only strong connections
         A_comp = A;
         A_comp(S_patt <0) = 0;
         tic;
         [P] = cpt_ProlMEX_EXTI(param.prolong,nc,fcnode_exti,S_patt,A,A_comp);
      otherwise
         err_msg = [char(param.prolong.proltype) ' is not a valid key for prolongation'];
         error(err_msg);
   end
else
   % Two stage prolongation
   switch upper(param.prolong.proltype)
      case 'BAMG'
         [P,c_mark] = cpt_ProlMEX_BAMG(level,param.prolong,nc,fcnode,S,TV);
         if verb
            fprintf('Number of nodes that are not perfectly interpolated: %d\n',sum(c_mark));
         end
         %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if sum(c_mark) > 0
            if param.prolong.updateCF == 1
               if verb
                  fprintf('Update the coarse nodes only\n');
               end
               % Update coarsening
               ind_new = find(c_mark==1);
               n_add = numel(ind_new);
               [ii,jj,pp] = find(P);
               ii = [ii; ind_new];
               jj = [jj; (nc+1:nc+n_add)'];
               pp = [pp; ones(n_add,1)];
               % Add coatse nodes to prolongation
               P = sparse(ii,jj,pp,numel(fcnode),nc+n_add);
               fcnode(ind_new) = (nc+1:nc+n_add)';
               nc = nc + n_add;
               clist = [clist; ind_new];
            elseif param.prolong.updateCF == 2
               if verb
                  fprintf('Update the coarse nodes and recompute prolongation\n');
               end
               % Update coarsening
               fcnode(fcnode>0) = 1;
               fcnode(fcnode<0) = 0;
               fcnode = fcnode + c_mark;
               nc = sum(fcnode);
               nf = numel(fcnode) - nc;
               fcnode(fcnode>0) = 1:nc;
               fcnode(fcnode==0) = -(1:nf);
               clist = find(fcnode>0);
               % Recompute prolongation
               [P,c_mark] = cpt_ProlMEX_BAMG(level,param.prolong,nc,fcnode,S,TV);
               if verb
                  fprintf('Number of nodes that are not perfectly interpolated: %d\n',sum(c_mark));
               end
            elseif (param.prolong.updateCF ~= 0)
               error('Wrong value of updateCF: %d\n',updateCF);
            end
         end
         %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         v = vecnorm(P');
         [amax,ind] = max(v);
         if verb
            fprintf('Max row norm: %15.6e\n',full(amax));
            fprintf('On row:       %d\n',ind);
         end
      case 'EXTI'
         % Create a compressed A with only strong connections
         A_comp = A;
         A_comp(S_patt <0) = 0;
         [P] = cpt_TwoStageProlMEX_ExtI(param.prolong,clist,fclist,fflist,A,S_patt,A_comp);
      otherwise
         err_msg = [char(param.prolong.proltype)...
                   ' is not a valid key for prolongation with aggressive coarsening'];
         error(err_msg);
   end
end
if verb
   fprintf('Number of nonzeroes per row in P: %f\n',(nnz(P)-size(P,2))/(size(P,1)-size(P,2)));
end

% If the problem is not symmetric force smooth instead of emin
prol_emin = param.prolong.prol_emin;
%if ~symm_flag && strcmpi(prol_emin,'EMIN')
%   if verb
%      fprintf('Since the problem is not symmetric smoothing is forced instead of EMIN\n');
%   end
%   prol_emin = 'smooth';
%end

% Prolongation smoothing
if strcmpi(prol_emin,'SMOOTH')
   if verb
      fprintf('Smooting Prolongation\n');
   end
   [nn,~] = size(P);
   if verb
      fprintf('Unsmoothed prolongation non-zeroes per row: %10.2f\n',(nnz(P)-nc)/(nn-nc));
   end
   DA = diag(diag(A))\A;
   lambda_max = eigs(DA,1,'lm','Tolerance',1.e-3,'FailureTreatment','keep');
   omega = 1.9 / lambda_max;
   if verb
      fprintf('Max eigenvalue and omega: %10.2f %10.2f\n',lambda_max,omega);
   end
   P = P - omega*DA*P;
   if verb
      fprintf('Smoothed prolongation non-zeroes per row: %10.2f\n',(nnz(P)-nc)/(nn-nc));
   end
end

if verb
   fprintf('Computing initial energy\n');
end
energy = trace(P'*(A*P));
if verb
   fprintf('Energy before orthog: %e\n',energy);
end

% Minimize Energy in Prolongation
if strcmpi(prol_emin,'EMIN')

   if verb
      fprintf('Minimizing Energy (MEX)\n');
   end

   % Compute a sparser Strength of connection before EMIN if needed
   coarsen2 = param.coarsen;
   coarsen2.tau = param.prolong.patt_tau;
   coarsen2.tau = 0.01;
   coarsen2.SoC_type = 'DOM';
   [~,~,~,S_emin,~] = coarsen(coarsen2,A,smootherOp,TV,verb);
   if verb
      fprintf('\nDensity of S_emin over S: %10.2f\n',nnz(S_emin)/nnz(S));
   end
   [nn,~] = size(P);
   if verb
      fprintf('Initial prolongation non-zeroes per row: %10.2f\n',(nnz(P)-nc)/(nn-nc));
   end
   patt_pow = param.prolong.patt_pow;
   ntv = size(TV,2);
   Ppatt = mk_prolPatt(patt_pow,S_emin,P,fcnode,A,param.prolong.nnzr_max,ntv,verb);

   % Apply EMIN algorithm
   [P,info] = MEX_EMIN_enhance_prol(level,param.prolong,symm_flag,A,Ppatt,P,TV,fcnode,verb);
   if verb
      fprintf('Ener. Min. prolongation non-zeroes per row: %10.2f\n',(nnz(P)-nc)/(nn-nc));
   end
end

% Store in structure
if strcmpi(prol_emin,'EMIN')
   emin_info.time_gath_K   = info(1);
   emin_info.time_prec_K   = info(2);
   emin_info.time_gath_B   = info(3);
   emin_info.time_PCG      = info(4);
   emin_info.time_overhead = info(5);
   emin_info.time_glob     = info(6);
   emin_info.avg_lfil      = info(7);
   emin_info.iter          = info(8);
   emin_info.relres        = info(9);
   emin_info.nnzr          = (nnz(P) - size(P,2)) / ( size(P,1) - size(P,2) );
   emin_info.Energy_init   = 0.0;
   emin_info.Energy_final  = 0.0;
   emin_info.nnz_K         = info(10);
   emin_info.nnz_PK        = info(11);
   emin_info.nnz_ZQ        = info(12);
else
   emin_info = 0;
end

end
