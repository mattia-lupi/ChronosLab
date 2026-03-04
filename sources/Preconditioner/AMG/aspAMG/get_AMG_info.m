%-----------------------------------------------------------------------------------------
%
%  Function to get information of the AMG hierarchy
%
%-----------------------------------------------------------------------------------------
function [AMG_info] = get_AMG_info(AMG_hrc,A)

% Init structure
AMG_info = struct;

% Loop over levels
n_lev = 1;
grid_comp = 0;
oper_comp = 0;
curr_AMG_hrc = AMG_hrc;
while true

   if isstruct(curr_AMG_hrc.next)

      % General level
      if n_lev == 1
         if nargin == 2
            nn = size(A,1);
            nt = nnz(A);
         else
            nn = 0;
            nt = 0;
         end
      else
         nn = size(curr_AMG_hrc.A,1);
         nt = nnz(curr_AMG_hrc.A);
      end
      grid_comp = grid_comp  + nn;
      oper_comp = oper_comp  + nt;
      AMG_info(n_lev).nn     = nn;
      AMG_info(n_lev).nt     = nt;
      AMG_info(n_lev).lmax_A = curr_AMG_hrc.lmax_A;
      AMG_info(n_lev).nt_S   = curr_AMG_hrc.Snnz;
      AMG_info(n_lev).omega  = curr_AMG_hrc.omega;
      AMG_info(n_lev).nc     = size(curr_AMG_hrc.P,2);
      AMG_info(n_lev).nt_P   = nnz(curr_AMG_hrc.P);

   else

      % Last level
      if n_lev == 1
         if nargin == 2
            nn = size(A,1);
            nt = nnz(A);
         else
            nn = 0;
            nt = 0;
         end
      else
         nn = size(curr_AMG_hrc.A,1);
         nt = nnz(curr_AMG_hrc.A);
      end
      %AMG_info(n_lev).nn     = size(curr_AMG_hrc.L,1);
      AMG_info(n_lev).nn     = nn;
      AMG_info(n_lev).nt     = nt;
      AMG_info(n_lev).lmax_A = curr_AMG_hrc.lmax_A;
      AMG_info(n_lev).nt_L   = nnz(curr_AMG_hrc.L);
      grid_comp = grid_comp  + nn;
      oper_comp = oper_comp  + nt;
      break

   end

   % Prepare for the next level
   n_lev = n_lev + 1;
   curr_AMG_hrc = curr_AMG_hrc.next;

end
grid_comp = grid_comp / size(A,1);
oper_comp = oper_comp / nnz(A);
AMG_info(1).grid_comp = grid_comp;
AMG_info(1).oper_comp = oper_comp;

return
