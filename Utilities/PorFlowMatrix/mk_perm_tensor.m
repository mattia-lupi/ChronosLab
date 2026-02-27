function K_list = mk_perm_tensor(k_perm)

K_list = {};
for i = 1:size(k_perm,1);
   K = zeros(3);
   K(1,2) = k_perm(i,4);
   K(1,3) = k_perm(i,5);
   K(2,3) = k_perm(i,6);
   K = K + K';
   K = K + diag(k_perm(i,1:3));
   K_list{i} = K;
end

end
