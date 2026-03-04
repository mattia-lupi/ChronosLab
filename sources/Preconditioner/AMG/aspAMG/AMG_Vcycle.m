function x = AMG_Vcycle(AMG_hrc,A,b,x0)

% Get the next level of the hierarchy
next_AMG_hrc = AMG_hrc.next;

if ~isstruct(next_AMG_hrc)

   % Solve coarsest system
   if AMG_hrc.symm
      x = AMG_hrc.S*(AMG_hrc.L'\(AMG_hrc.L\(AMG_hrc.S'*b)));
   else
      x = AMG_hrc.Q*(AMG_hrc.U\(AMG_hrc.L\(AMG_hrc.P*b)));
   end

else

   % Presmooth
   if nargin < 4
      x = AMG_hrc.Minv1(b);
   else
      r = b - A*x0;
      x = x0 + AMG_hrc.Minv1(r);
   end
   for k = 2:AMG_hrc.nupre
      r = b - A*x;
      x = x + AMG_hrc.Minv1(r);
   end

   % Compute residual
   r = b - A*x;

   % Restrict residual
   r = AMG_hrc.P'*r;

   % Recursive call to the next level
   if isstruct(next_AMG_hrc.next)
      next_A = next_AMG_hrc.A;
   else
      next_A = 0;
   end
   hc = AMG_Vcycle(next_AMG_hrc,next_A,r);

   % Prolong and correct solution
   x = x + AMG_hrc.P*hc;

   % Post-smooth
   for k = 1:AMG_hrc.nupost
      r = b - A*x;
      x = x + AMG_hrc.Minv2(r);
   end

end

return
