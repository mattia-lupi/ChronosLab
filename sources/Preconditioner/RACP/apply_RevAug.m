function y = apply_RevAug(AMG_prec,AA,B,inv_C,x)

n1 = size(AA,1);
n2 = size(inv_C,1);

x1 = x(1:n1,:);
x2 = x(n1+1:end,:);

%  AA y1 + B y2 = x1 ---> 2 ---> [AA + B*inv(C)*B'] y1 = x1 +B*inv(C)*x2
%  B' y1 - C y2 = x2 ---> 1 ---> y2 = inv(C)*(B' y1 - x2) 

b1 = x1 + B*(inv_C*x2);
if true
   y1 = AMG_Vcycle(AMG_prec,AA,b1);
else
   tol = 1.e-8;
   itmax = 100;
   PREC = @(x) AMG_Vcycle(AMG_prec,AA,x);
   [y1,flag,relres,iter,resvec] = pcg(AA,b1,tol,itmax,PREC);
end
y2 = inv_C*(B'*y1 - x2);

y = [y1; y2];

end
