function [fcnode, clist, flist] = cpt_PMIS(S, fcnode,verb)

%  0 --> not assigned
% -1 --> FINE
%  1 --> COARSE

% Get the number of vertices
nn = size(S,1);

% Remove diagonal from S
S = S - spdiags(diag(S), 0, nn, nn);

% Create the undirected graph
S_sym = abs(S) + abs(S');

% Compute initial weights (number of nodes influenced by)
nod_weight = sum(S ~= 0, 1)';

% Set Dirichlet nodes as FINE
fcnode(fcnode == -2) = -1;

% Add a random number to avoid ties
nod_weight = nod_weight + rand(nn,1);

% Mask free weights (-1 indicates assigned)
w_free = nod_weight;
w_free(fcnode == -1) = -1;

% Initialize the list of nodes to assign
list = find(fcnode == 0);
len_list = numel(list);

tic;
% Loop until the list is not empty
ind = 1;
while len_list > 0

   % Get a node from the list
   inod = list(ind);

   % Check that the node is free
   if fcnode(inod) ~= 0

      % Remove it from the search list
      list(ind) = list(len_list);
      len_list = len_list - 1;

   else

      % Check if it can become coarse
      is_coarse = true;
      inod_wgt = nod_weight(inod);

      % Get the strong neighbours of inod
      j_neigh = find(S_sym(:,inod));

      if any(w_free(j_neigh) > inod_wgt)

         % One of its strong-free neighbours has a larger weight ==> inod cannot be coarse
         is_coarse = false;
         ind = ind + 1;

      end

      if is_coarse

         % This node will become a coarse node
         % Mark the strongly influenced neighbours of inod as Fine
         j_neigh = find(S(:,inod));
         fcnode(j_neigh) = -1;
         w_free(j_neigh) = -1;

         % Mark inod as coarse
         fcnode(inod) = 1;
         w_free(inod) = -1;

         % Remove it from the list
         list(ind) = list(len_list);
         len_list = len_list - 1;

      end

   end

   % If ind is out of the list length, start from the beginning
   if ind > len_list
      ind = 1;
   end

end
time = toc;
if verb
   fprintf('Time for MIS: %f\n',time);
end

% Count and number coarse nodes (vectorized)
clist = find(fcnode > 0);
flist = find(fcnode < 0);
fcnode(clist) = 1:numel(clist);
fcnode(flist) = -(1:numel(flist));

return
