function TV_out = cpt_initApp(k_type,FT,TV_in,verb)
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
if nargin < 4
   verb = 1;
end
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
d_inv = 1 ./ diag(FT);

% If FT is triangular or diagonal (e.g. FSAI factor transpose or Jacobi),
% the eigenvalues of D_inv * FT are identically 1.0 on the diagonal, so lmax = 1.0.
% Otherwise compute lmax using an operator function without forming D_inv * FT.
if istriu(FT) || istril(FT) || isdiag(FT)
   omega = 1.0;
else
   D_inv_op = @(x) d_inv .* (FT * x);
   opts.disp = verb;
   opts.tol = 1.e-3;
   [~,lmax,~] = eigs(D_inv_op, size(FT, 1), 1, 'lm', opts);
   omega = min(1.0, 1.9 / abs(lmax));
end

scale = omega .* d_inv;
TV_out = scale .* TV_in;
for i = 2:k_type
   TV_out = TV_out + scale .* (TV_in - FT*TV_out);
end


end
