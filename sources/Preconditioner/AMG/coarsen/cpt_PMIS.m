function [fcnode, clist, flist] = cpt_PMIS(S, fcnode,verb)

%  0 --> not assigned
% -1 --> FINE
%  1 --> COARSE

% Get the number of vertices
nn = size(S,1);

% Remove diagonal from S
S = S - diag(diag(S));

% Create the undirected graph
S_sym = abs(S) + abs(S');

% Initialize the list of nodes to assign
list = (1:nn)';

% Compute initial weigths (number of nodes influenced by)
nod_weight = sum(S ~= 0)';

% Set Dirichlet nodes and non influencing nodes as FINE
list(fcnode == -2) = nn + 100;
fcnode(fcnode == -2) = -1;
%list(nod_weight == 0) = nn + 100;
%fcnode(nod_weight == 0) = -1;
%%%%%%%%%%%%%%%%%%%
%save('FFF','fcnode','nod_weight','S','list','nn');
%%%%%%%%%%%%%%%%%%%

% Add a random number to avoid ties
nod_weight = nod_weight + rand(nn,1);

% Sort the list and get its effective length
list = sort(list);
len_list = find(list == nn + 100,1,'first');
if isempty(len_list)
   len_list = nn;
else
   len_list = len_list - 1;
end

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

      if sum(fcnode(j_neigh) == 0 & nod_weight(j_neigh) > inod_wgt) > 0

         % One of its strong-free neighbours has a larger weight ==> inod cannot be coarse
         is_coarse = false;
         ind = ind + 1;

      end

      if is_coarse

         % This node will become a coarse node
         % Mark the strongly influenced neighbours of inod as Fine
         j_neigh = find(S(:,inod));
         fcnode(j_neigh) = -1;

         % Mark inod as coarse
         fcnode(inod) = 1;

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

% Count and number coarse nodes
nc = 0;
nf = 0;
clist = zeros(nn,1);
flist = zeros(nn,1);
for i = 1:nn
   if fcnode(i) > 0
      nc = nc + 1;
      fcnode(i) = nc;
      clist(nc) = i;
   elseif fcnode(i) < 0
      nf = nf + 1;
      fcnode(i) = -nf;
      flist(nf) = i;
   end
end
clist = clist(1:nc);
flist = flist(1:nf);

return
