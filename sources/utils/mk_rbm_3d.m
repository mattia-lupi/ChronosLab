function RBM = mk_rbm_3d(X)

nn = size(X,1);

RBM = zeros(3*nn,6);

RBM(1:3:end,1) = 1;
RBM(2:3:end,2) = 1;
RBM(3:3:end,3) = 1;

XB = sum(X) ./ nn;

ind = 0;
for i = 1:nn
   RBM(ind+2,4) = -X(i,3) + XB(3);
   RBM(ind+3,4) =  X(i,2) - XB(2);

   RBM(ind+1,5) =  X(i,3) - XB(3);
   RBM(ind+3,5) = -X(i,1) + XB(1);

   RBM(ind+1,6) = -X(i,2) + XB(2);
   RBM(ind+2,6) =  X(i,1) - XB(1);

   ind = ind + 3;
end

end
