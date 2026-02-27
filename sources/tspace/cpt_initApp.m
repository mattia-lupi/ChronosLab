function TV_out = cpt_initApp(k_type,FT,TV_in)
%-----------------------------------------------------------------------------------------
%
% Function to set-up the initial solution in case of LANCZOS and SRQCG. Both methods try
% to solve the generalized eigenproblem A*v = lambda*M*v by splitting the inverse of M
% around A. Thus, if a near kernel is known for A (e.g. A*v0 simeq 0), it must be
% pre-processed to become a tentative near kernel for F*A*FT (with FT*F = inv(M)).
% The new tentative solution is computed as w0 simeq FT\v0.
% 
% ktype <= -2   ==>   random w0
%       == -1   ==>   w0 = v0
%       ==  0   ==>   w0 = FT\v0
%       ==  1   ==>   w0 = diag(FT)\v0
%       >=  2   ==>   perform ktype steps of Jacobi
%
%-----------------------------------------------------------------------------------------

if k_type == 0
   TV_out = FT\TV_in;
   return
elseif k_type == -1
   TV_out = TV_in;
   return
elseif k_type < 0
   TV_out = randn(size(TV_in));
   return
end

% Compute diagonal preconditioner
D_inv = 1 ./ diag(FT);
D_inv = diag(sparse(D_inv));
lmax = eigs(D_inv*FT,1,'lm','Display',1,'Tolerance',1.e-3);
omega = min(1.0,1.9 / abs(lmax));
TV_out = omega*D_inv*TV_in;
for i = 2:k_type
   TV_out = TV_out + omega*(D_inv*(TV_in - FT*TV_out));
end


end
