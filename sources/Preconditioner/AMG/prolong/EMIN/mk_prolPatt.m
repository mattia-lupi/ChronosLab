function P_patt = mk_prolPatt(kpow,S,P,fcnode,A,nnzr,nnzr_min,verb)
%-----------------------------------------------------------------------------------------
%
% Creates the power pattern for prolongation set-up. The output pattern do not include
% coarse nodes.
%
%-----------------------------------------------------------------------------------------

% Flag wether to use or not a unitary prolongation
FLAG_unitary = false;

% Extract dimensions
[nn,mm] = size(P);
nf = nn - mm;

% Init prolongation pattern
P_patt = P;
if FLAG_unitary
   [ii,jj,pp] = find(P_patt);
   pp = 1;
   P_patt = sparse(ii,jj,pp,nn,mm);
end
if verb
   fprintf('Power 1 nnzr: %10.2f\n',(nnz(P_patt)-mm) / (nn-mm));
end

if nargin == 4

   % Compute standard prolongation pattern
   S(S ~= 0) = 1;
   for ipow = 2:kpow
      P_patt = S*P_patt;
      if verb
         fprintf('Power %d nnzr: %10.2f\n',ipow,nnz(P_patt)/nn);
      end
   end
   % Remove coarse nodes
   ind_c = find(fcnode > 0);
   tmp = P_patt';
   tmp(:,ind_c) = 0;
   P_patt = tmp';
   if verb
      fprintf('Final nnzr: %10.2f\n',(nnz(P_patt)-mm) / (nn-mm));
   end
else

   % Save initial pattern
   P_patt_0 = P_patt;

   % Compute prolongation pattern with maximum fill-in
   nnzr = max(nnzr,nnzr_min);
   for ipow = 2:kpow
      P_patt = S*P_patt;
      P_patt = drop_patt(nnzr,P_patt,P_patt_0);
      if verb
         fprintf('Power %d nnzr: %10.2f\n',ipow,nnz(P_patt)/nf);
      end
   end
   % Remove coarse nodes
   ind_c = find(fcnode > 0);
   tmp = P_patt';
   tmp(:,ind_c) = 0;
   P_patt = tmp';
   if verb
      fprintf('Final nnzr: %10.2f\n',(nnz(P_patt)-mm) / (nn-mm));
   end

   % Make sure that each prolongation row has at least nnzr_min entries
   % Use matrix A to increase connectivity

   % Find poorly connected nodes
   [ii,jj,pp] = find(P_patt);
   pp = 1;
   nn_P = size(P_patt,1);
   count = zeros(nn_P,1);
   for i = 1:numel(ii)
      count(ii(i)) = count(ii(i)) + 1;
   end
   count(fcnode>0) = nnzr_min + 100;
   ind_poor = ones(size(count));
   ind_poor(count ==0) = 0;
   ind_poor(count >= nnzr_min) = 0;
   ind_poor = find(ind_poor)';
   nn_poor = numel(ind_poor);

   % If there are poorly connected nodes
   if verb
      fprintf('There are %d poorly connected nodes\n',nn_poor);
   end
   %if (nn_poor > 0)
   if (nn_poor > 0 && false) %@@@@@@@@ STARE ATTENTI A QUESTA CORREZIONE
      of = fopen('BAD_LIST','w');

      % Increase connectivity of each node
      PT = P_patt';
      ind_list = 1;
      ii_add = zeros(3*nn_poor*nnzr_min,1);
      jj_add = zeros(3*nn_poor*nnzr_min,1);
      % Loop over poorly connected nodes
      for i = 1:nn_poor
         inod = ind_poor(i);
         % Get the list of currently connected nodes
         [list,~,~] = find(PT(:,inod));
         % Get the list of neighbors through A
         [A_list,~,~] = find(A(:,ind_poor(i)));
         k = 0;
         while numel(list) < nnzr_min
            % Get entries to add
            [add,~,~] = find(PT(:,A_list));
            % Update the list
            list = union(list,add);
            % Get a new list of neighbors
            [new_A_list,~,~] = find(A(:,A_list));
            A_list = setdiff(unique(new_A_list),A_list);
            % Check that the number of while iterations does not increase too much
            k = k + 1;
            if k > 10
               fprintf('ERROR for node %d, unable to create connectivity within %d loops\n',...
                       inod,k);
               fprintf(of,'%6d %6d\n',inod,nnz(PT(:,inod)));
               break
            end
         end
         nn_list = numel(list);
         ii_add(ind_list:ind_list+nn_list-1) = inod;
         jj_add(ind_list:ind_list+nn_list-1) = list;
         ind_list = ind_list + nn_list;
      end
      nn_add = ind_list-1;
      % Merge new indices with old pattern
      ii_add = ii_add(1:nn_add);
      jj_add = jj_add(1:nn_add);
      aa_add = ones(nn_add,1);
      [ii,jj,aa] = find(P_patt);
      ii_new = [ii; ii_add];
      jj_new = [jj; jj_add];
      aa_new = [aa; aa_add];
      % Create the new pattern matrix
      [nn,mm] = size(P);
      P_patt = sparse(ii_new,jj_new,aa_new,nn,mm);

   end

end

end
