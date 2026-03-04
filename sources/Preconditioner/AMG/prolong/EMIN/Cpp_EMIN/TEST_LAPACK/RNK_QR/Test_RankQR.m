clc
clear
%rng('default');

nn = 8;
mm = 6;
rnk = 5;
% Create a rank deficient matrix
B = rand(nn,mm);
[U,S,V] = svd(B);
for i = rnk+1:min(nn,mm)
   S(i,i) = 0;
end
B = U*S*V';

% Choose g in range of B'
tmp = rand(nn,1) .* 10.^(6*rand(nn,1)-3);
g = B'*tmp;

p0 = rand(nn,1);

z = rand(nn,1);

% Dump B
of = fopen('MATRIX','w');
fprintf(of,'%d %d\n',nn,mm);
for i = 1:nn
   fprintf(of,' %20.11e',B(i,:));
   fprintf(of,'\n');
end
fprintf(of,'%20.11e\n',g);
fprintf(of,'%20.11e\n',p0);
fprintf(of,'%20.11e\n',z);
fclose(of);

% Perform a rank revealing QR fact of B' ==> B*P = Q*R
[Q,R,P] = qr(B');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diag(R) / R(1,1)
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q_rid = Q(:,1:rnk);
R_rid = R(1:rnk,:);
res = g - B'*p0;
Delta = P*(R_rid'*((R_rid*R_rid')\(Q_rid'*res)));
p2 = p0+Delta;
fprintf('|B''*p2-g|: %e - |p2|: %e\n',norm(B'*p2-g), norm(p2));

%R_rid
%X = R_rid*R_rid'
%Q_rid'*res
%(R_rid*R_rid')\(Q_rid'*res)
%R_rid'*((R_rid*R_rid')\(Q_rid'*res))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U = chol(R_rid*R_rid'); % inv(U'*U) = inv(U)*inv(U');
%(U')\R_rid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orthogonal projector
RR = R_rid*P';
Orth = eye(size(B,1)) - RR'*(inv(RR*RR'))*RR;
%[tmp,~] = qr(R_rid',0); tmp = tmp';
%O2 = eye(size(B,1)) - tmp'*(inv(tmp*tmp'))*tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = chol(RR*RR'); % inv(U'*U) = inv(U)*inv(U');
(U')\RR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BT = B';

% Get something that may give a non-zero B'*z
norm(B'*z)
z_new = Orth*z;
norm(B'*z_new)


