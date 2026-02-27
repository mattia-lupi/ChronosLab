function V = normc(V)

for i = 1:size(V,2)
   V(:,i) = V(:,i) / norm(V(:,i));
end

return
