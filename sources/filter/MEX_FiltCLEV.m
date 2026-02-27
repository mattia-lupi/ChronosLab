%-----------------------------------------------------------------------------------------
% patt_min_flag:  Use or not the minimal pattern
% tau:            dropping tolerance
%-----------------------------------------------------------------------------------------
function AC = MEX_FiltCLEV(np,patt_min_flag,tau,A0,A,fcnode,P,TV)

% Get dimensions
[nf,nc] = size(P);

if patt_min_flag
   % Create simple injection operator
   PI = sparse(find(fcnode>0),fcnode(fcnode>0),ones(nc,1));
   if size(PI,1) < nf
      PI(nf,nc) = 0;
   end
   % Create minimal pattern
   Patt_min = P'*A0*PI;
   Patt_min=0.5*(Patt_min+Patt_min');
   Patt_min(Patt_min~=0) = 1;
   Patt_min = Patt_min - speye(nc);
else
   Patt_min = sparse(nc,nc);
end

% Filter and compensate rows
AC = MEX_FiltComp(A,Patt_min,np,tau,TV);

end
