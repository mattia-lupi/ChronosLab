function y = apply_MCP(MCP_prec,A,B,C,x)

   n1 = size(A,1);
   n2 = size(C,1);

   % Apply preconditioner for 11 block
   y1 = AMG_Vcycle(MCP_prec.AMG_prec,A,x(1:n1));

   % Transfer on the second block
   tmp2 = B'*y1 - x(n1+1:n1+n2);

   % Apply Schur complement preconditioner
   if true
      G22 = MCP_prec.G22;
      y2 = G22'*(G22*tmp2);
   else
      y2 = -C\tmp2;
   end

   % Transfer on the first block
   y1 = B*y2;
   tmp1 = x(1:n1) - y1;

   % Apply preconditioner for 11 block
   y1 = AMG_Vcycle(MCP_prec.AMG_prec,A,tmp1);

   % Compose solution
   y = [y1;y2];

end
