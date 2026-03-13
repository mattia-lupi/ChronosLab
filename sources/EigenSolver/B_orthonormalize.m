function [X_out,BX_out,U] = B_orthonormalize(prodB,X_in)

if isempty(prodB)
   BX = X_in;
else
   BX = prodB(X_in);
end
XBX = X_in'*BX; XBX = 0.5*(XBX+XBX)';
cond_XBX = cond(XBX);
if cond_XBX > 1/(max(size(XBX))*eps)
   fprintf('B_orthonormalize: XBX highly ill-conditioned with cond %15.6e\n',cond_XBX);
   X_out = [];
   BX_out = [];
   U = [];
   return;
end
U = chol(XBX);
X_out = X_in/U;
if ~isempty(prodB)
   BX_out = BX/U;
else
   BX_out = [];
end

end
