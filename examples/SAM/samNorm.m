function normm = samNorm(Ak,A0)
   AkT = Ak';
   A0T = A0';
   normCol = NaN(size(A0,1),1);

   for i = 1:size(A0,1)
      normCol(i) = 2*norm(AkT(:,i) - A0T(:,i))/(norm(AkT(:,i)) + norm(A0T(:,i)));
   end
   
   normm = mean(normCol);
end