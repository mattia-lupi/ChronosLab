function [normm] = normSAM(Ak,A0)
normCol = zeros(size(A0,1),1);
for k = 1:size(A0,1)
   normCol(k) = 2*norm(Ak(:,k) - A0(:,k)) / (norm(Ak(:,k)) + norm(A0(:,k)));
end

normm = mean(normCol);