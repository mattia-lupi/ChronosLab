function A = cpt_loc_ANISOMAT(tetra,coord,K)

% Form the basis matrix
updt_idx = @(p) mod(p,4)+1;
B = [1 coord(tetra(1),:); 1 coord(tetra(2),:); 1 coord(tetra(3),:); 1 coord(tetra(4),:)];

Delta = det(B);
a = zeros(4,1);
b = zeros(4,1);
c = zeros(4,1);
d = zeros(4,1);
pr = [2,3,4];
pc = [2,3,4];
fac = 1.0;
for i = 1:4
   a(i) =  det(B(pr,pc))*fac; pc = updt_idx(pc);
   b(i) = -det(B(pr,pc))*fac; pc = updt_idx(pc);
   c(i) =  det(B(pr,pc))*fac; pc = updt_idx(pc);
   d(i) = -det(B(pr,pc))*fac; pc = updt_idx(pc);
   pr = updt_idx(pr);
   fac = -fac;
end

A = zeros(4);
for i = 1:4
   for j = 1:4
      A(i,j) = [b(i) c(i) d(i)]*K*[b(j); c(j); d(j)];
   end
end
A = A / abs(6*Delta);

end
