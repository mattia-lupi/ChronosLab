function [S_Fc] = mk_SoC_FcC(S, Fclist,fcnode)

   % Get the number of coarse nodes obtained after the first stage 
   tic;
   nFc = numel(Fclist);

   % Remove diagonal from S
   S = S - diag(diag(S));

   % Create the undirected graph
   S_sym = abs(S) + abs(S');

   % Get the square of S
   S2 = S_sym*S_sym;

   % Restric connection graph to Fclist only
   opt1 = false;
   if opt1
      S2_fc = S2(:,Fclist);
      S2_fc = S2_fc';
      S_Fc = S2_fc(:,Fclist);
   else
      S_Fc = S2(Fclist,Fclist);
   end

end
