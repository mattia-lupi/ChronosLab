function [S_Fc] = mk_SoC_FcC(S, Fclist,fcnode)

   % Remove diagonal from S
   S = S - spdiags(diag(S), 0, size(S,1), size(S,2));

   % Create the undirected graph
   S_sym = abs(S) + abs(S');

   % Restrict connection graph to Fclist directly without full S2 allocation
   S_sub = S_sym(:, Fclist);
   S_Fc = S_sub' * S_sub;

end
