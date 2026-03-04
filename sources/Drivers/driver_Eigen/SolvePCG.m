function x = SolvePCG(A,PREC,rhs,itmax,tol)

   global count;
   global tot_iter;
   prodA = @(x) A*x;
   [x,~,~,iter] = pcg(prodA,rhs,tol,itmax,PREC);
   count = count + 1;
   tot_iter = tot_iter + iter;
   %fprintf('%d %d\n',count,iter);

end
