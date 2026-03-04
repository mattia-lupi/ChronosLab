%-----------------------------------------------------------------------------------------
%
% Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG).
%
% References
% ----------  
% [1] A. V. Knyazev (2001),
%     Toward the Optimal Preconditioned Eigensolver: Locally Optimal
%     Block Preconditioned Conjugate Gradient Method.
%     SIAM Journal on Scientific Computing 23, no. 2,
%     pp. 517-541. :doi:`10.1137/S1064827500366124`
%                                 
% [2] A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov (2007),
%     Block Locally Optimal Preconditioned Eigenvalue Xolvers (BLOPEX)
%     in hypre and PETSc.  :arxiv:`0705.2626`
%    
% [3] A. V. Knyazev's C and MATLAB implementations:
%     https://github.com/lobpcg/blopex
%-----------------------------------------------------------------------------------------

function [niter,lambda,blockX,resnorm_vec,lambda_vec,ierr] = ...
         lobpcg(prodA,prodB,prodM,blockY,blockX0,largest_flag,...
                reslambda_check,itmax,tol,restartControl,verb)

if nargin < 11
   verb = 1;
end
global DEBUG

% Init error code
ierr = 0;

% Extract blockX0 size
nn = size(blockX0,1);

% Extract number of desired eigenpairs
neig = size(blockX0,2);

% Initialize blockX
blockX = blockX0;
%@@@@@@@@@@@@@@@@@@@@@@@@
if DEBUG
   prt_mat('Initial blockX',blockX);
end
%@@@@@@@@@@@@@@@@@@@@@@@@

% Initialize lambda old
lambda_old = zeros(neig,1);

% Check inputs
B_flag = ~isempty(prodB);
M_flag = ~isempty(prodM);
Y_flag = ~isempty(blockY);

% Init the best results
smallest_resnorm = realmax;
best_blockX = blockX;

% Init the vectors of residual norms and lambda history
resnorm_vec = zeros(itmax+1,neig);
lambda_vec = zeros(itmax+1,neig);

% Apply constraints to blockX
if Y_flag
   if B_flag
      blockBY = prodB(blockY);
   else
      blockBY = blockY;
   end
   YBY = blockY'*blockBY;
   U_YBY = chol(YBY);
   % Define the B-orthogonal projector of Y
   orthog_Y = @(x) x - blockY*(U_YBY\(U_YBY'\(blockBY'*x)));
   blockX = orthog_Y(blockX);
end 
%@@@@@@@@@@@@@@@@@@@@@@@@
if DEBUG
   prt_mat('Constrained blockX',blockX);
end
%@@@@@@@@@@@@@@@@@@@@@@@@

% B-orthonormalize X
[blockX,blockBX,~] = B_orthonormalize(prodB,blockX);

% Check failure in B-orthonormalization
if isempty(blockX)
   fprintf('Initial approximation is not full-rank\n');
   ierr = 1;
end

% Compute initial Ritz vectors
blockAX = prodA(blockX);
XAX = blockX'*blockAX; XAX = 0.5*(XAX+XAX');
if B_flag
   XBX = blockX'*blockBX; XBX = 0.5*(XBX+XBX');
else
   XBX = blockX'*blockX; XBX = 0.5*(XBX+XBX');
end
[blockEig,lambda] = eig(XAX,XBX);
lambda_vec(1,:) = diag(lambda);

% Update blockX, blockAX and blockBX
blockX = blockX*blockEig;
blockAX = blockAX*blockEig;
if B_flag
   blockBX = blockBX*blockEig;
end

% Set active set mask (all active)
activeMask = ones(neig,1);

%@@@@@@@@@@@@@@@@@@@@@@@@@@@
%save('XXX');
%error('AAA');
%@@@@@@@@@@@@@@@@@@@@@@@@@@@

% Allocate blockP, blockAP, blockBP, blockR
blockP = zeros(nn,neig);
blockAP = zeros(nn,neig);
if B_flag
   blockBP = zeros(nn,neig);
end
blockR = zeros(nn,neig);
ActBlockR = zeros(nn,neig);
ActBlockP = zeros(nn,neig);
ActBlockAP = zeros(nn,neig);
if B_flag
   ActBlockBP = zeros(nn,neig);
end
ActBlockBR = zeros(nn,neig);
blockAR = zeros(nn,neig);
ActBlockAP = zeros(nn,neig);

% Init the main loop
if verb
   fprintf('%4s | %5s | %18s | %18s\n','iter','# act','avg resnorm','max delta lambda');
end
iter = 0;
restart = true;
forcedRestart = false;
while iter < itmax-1 % MAIN LOOP START
   iter = iter + 1;

   % Compute residual
   if B_flag
      blockR = blockAX - blockBX*lambda;
   else
      blockR = blockAX - blockX*lambda;
   end
   
   % Compute residual norm
   resnorms = sqrt(diag(blockR'*blockR))./diag(lambda);
   resnorm_vec(iter,:) = resnorms;
   avg_resnorm = mean(resnorms);

   % Compute lambda relative delta
   lambda_relDelta = abs(lambda_old-diag(lambda))./diag(lambda);
   lambda_old = diag(lambda);

   if avg_resnorm < smallest_resnorm
      % Update best residual norm
      smallest_resnorm = avg_resnorm;
      best_blockX = blockX;
   else
      % Check if restart is required
      if avg_resnorm > (2^restartControl)*smallest_resnorm
         forcedRestart = true;
         blockAX = prodA(blockX);
         if B_flag
            blockBX = prodB(blockX);
         end
      end
   end

   % Update the active mask
   if reslambda_check
      ind = find(lambda_relDelta < tol);
   else
      ind = find(resnorms < tol);
   end
   activeMask(:) = 1;
   activeMask(ind) = 0;
   n_active = neig-numel(ind);

   % Dump the number of active vectors
   if mod(iter-1,1) == 0
      if verb
         fprintf('%4d | %5d | %18.10e | %18.10e\n',iter-1,n_active,avg_resnorm,...
               max(lambda_relDelta));
      end
   end

   % If no active vectors, break
   if n_active == 0
      break
   end

   % Create subblocks
   ind = find(activeMask);
   ActBlockR = blockR(:,ind);
   if iter > 1
      ActBlockP = blockP(:,ind);
      ActBlockAP = blockAP(:,ind);
      if B_flag
         ActBlockBP = blockBP(:,ind);
      end
   end
   %@@@@@@@@@@@@@@@@@@@@@@@@@@@
   if DEBUG
      fprintf('PRIMA PREC\n');
      prt_mat('ActBlockR',ActBlockR);
   end
   %@@@@@@@@@@@@@@@@@@@@@@@@@@@

   % Apply preconditioner to active residuals
   if M_flag
      ActBlockR = prodM(ActBlockR);
   end
   %@@@@@@@@@@@@@@@@@@@@@@@@@@@
   if DEBUG
      fprintf('DOPO PREC\n');
      prt_mat('ActBlockR',ActBlockR);
   end
   %@@@@@@@@@@@@@@@@@@@@@@@@@@@

   % Apply constraint to preconditioned residuals
   if Y_flag
      ActBlockR = orthog_Y(ActBlockR);
   end

   % B-orthogonalize the preconditioned residuals to X
   if B_flag
      ActBlockR = ActBlockR - blockX*(blockBX'*ActBlockR);
   else
      ActBlockR = ActBlockR - blockX*(blockX'*ActBlockR);
   end

   % B-orthonormalize the preconditioned residuals
   [ActBlockR,ActBlockBR,~] = B_orthonormalize(prodB,ActBlockR);
   
   % Check failure in B-orthonormalization
   if isempty(ActBlockR)
      fprintf('Failure in B-orthonormalization of blockR at iteration %d\n',iter);
      ierr = 2;
      break
   end
   %@@@@@@@@@@@@@@@@@@@@@@@@@@@
   if DEBUG
      prt_mat('ActBlockR dopo B-ortho',ActBlockR)
   end
   %@@@@@@@@@@@@@@@@@@@@@@@@@@@

   ActBlockAR = prodA(ActBlockR);

   % Update blockP
   if iter > 1
      if B_flag
         [ActBlockP,ActBlockBP,U_PBP] = B_orthonormalize(prodB,ActBlockP);
      else
         [ActBlockP,~,U_PBP] = B_orthonormalize(prodB,ActBlockP);
      end
      % Check failure in B_orthonormalization
      if isempty(ActBlockP)
         % Restart iteration
         restart = true;
      else
         ActBlockAP = ActBlockAP/U_PBP;
         restart = forcedRestart;
      end
   end
   %@@@@@@@@@@@@@@@@@@@@@@@@
   if restart
      if verb
         fprintf('restart %d\n',restart);
      end
   end
   %@@@@@@@@@@@@@@@@@@@@@@@@

   %----------------------
   % Perform Rayleigh Ritz
   %----------------------
   if ~B_flag
      blockBX = blockX;
      ActBlockBR = ActBlockR;
      if ~restart
         ActBlockBP = ActBlockP;
      end
   end

   % Compute Gram matrices
   gramXAR = blockX'*ActBlockAR;
   gramRAR = ActBlockR'*ActBlockAR; gramRAR = 0.5*(gramRAR+gramRAR');
   gramXAX = blockX'*blockAX; gramXAX = 0.5*(gramXAX+gramXAX');
   gramXBX = blockX'*blockBX; gramXBX = 0.5*(gramXBX+gramXBX');
   gramRBR = ActBlockR'*ActBlockBR; gramRBR = 0.5*(gramRBR+gramRBR');
   gramXBR = blockX'*ActBlockBR;
   if ~restart
      gramXAP = blockX'*ActBlockAP;
      gramRAP = ActBlockR'*ActBlockAP;
      gramPAP = ActBlockP'*ActBlockAP; gramPAP = 0.5*(gramPAP+gramPAP');
      gramXBP = blockX'*ActBlockBP;
      gramRBP = ActBlockR'*ActBlockBP;
      gramPBP = ActBlockP'*ActBlockBP; gramPBP = 0.5*(gramPBP+gramPBP');
      gramA = [ gramXAX , gramXAR , gramXAP;
                gramXAR', gramRAR , gramRAP;
                gramXAP', gramRAP', gramPAP ];
      gramB = [ gramXBX , gramXBR , gramXBP;
                gramXBR', gramRBR , gramRBP;
                gramXBP', gramRBP', gramPBP ];

   else
      gramA = [ gramXAX , gramXAR;
                gramXAR', gramRAR ];
      gramB = [ gramXBX , gramXBR;
                gramXBR', gramRBR ];

   end
   %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   if DEBUG
      prt_mat('gramA',gramA);
      prt_mat('gramB',gramB);
   end
   %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   % Compute Ritz eigenpairs
   [blockEig,lambda] = eig(gramA,gramB);
   %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   if DEBUG
      fprintf('Dopo Ritz:');
      prt_mat('lambda',diag(lambda));
      prt_mat('blockEig',blockEig);
      if (iter==3)
           %error('XXX')
      end
   end
   %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   % Select the neig largest (smallest)
   lambda = diag(lambda);
   if largest_flag
      ind = numel(lambda):-1:numel(lambda)-neig+1;
   else
      ind = 1:neig;
   end
   lambda = lambda(ind);
   blockEig = blockEig(:,ind);
   % Store lambda history
   lambda_vec(iter+1,:) = lambda;
   % Make lambda diagonal again
   lambda = diag(lambda);

   % Compute Ritz vectors
   if ~restart
      blockEigX = blockEig(1:neig,:);
      blockEigR = blockEig(neig+1:neig+n_active,:);
      blockEigP = blockEig(neig+n_active+1:end,:);
      pp = ActBlockR*blockEigR + ActBlockP*blockEigP;
      app = ActBlockAR*blockEigR + ActBlockAP*blockEigP;
      if B_flag
         bpp = ActBlockBR*blockEigR + ActBlockBP*blockEigP;
      end
   else
      blockEigX = blockEig(1:neig,:);
      blockEigR = blockEig(neig+1:end,:);
      pp = ActBlockR*blockEigR;
      app = ActBlockAR*blockEigR;
      if B_flag
         bpp = ActBlockBR*blockEigR;
      end
   end
   blockX = pp + blockX*blockEigX;
   blockAX = app + blockAX*blockEigX;
   blockP = pp;
   blockAP = app;
   if B_flag
      blockBX = bpp + blockBX*blockEigX;
      blockBP = bpp;
   end

end % MAIN LOOP END

if n_active > 0
   % There were still some active lambda
   iter = iter + 1;

   % Compute residuals
   if B_flag
      blockR = blockAX - blockBX*lambda;
   else
      blockR = blockAX - blockX*lambda;
   end

   % Compute residual norms
   resnorms = sqrt(diag(blockR'*blockR))./diag(lambda);

   % Store lambda history and residuals
   resnorm_vec(iter,:) = resnorms;
   lambda_vec(iter,:) = diag(lambda);

end

% Check if convergence has been achieved
if reslambda_check
   if max(lambda_relDelta) > tol
      ierr = -1;
      %retrieve best_blockX
      blockX = best_blockX;
   end
else
   if max(resnorms) > tol
      ierr = -1;
      %retrieve best_blockX
      blockX = best_blockX;
   end
end

%-----------------------------------------
% Last iteration (without preconditioning)
%-----------------------------------------

iter = iter + 1;

% Enforce constraint again
if Y_flag
   blockX = orthog_Y(blockX);
end

% Apply Rayleigh-Ritz again
blockAX = prodA(blockX);
XAX = blockX'*blockAX; XAX = 0.5*(XAX+XAX');
if B_flag
   blockBX = prodB(blockX);
   XBX = blockX'*blockBX; XBX = 0.5*(XBX+XBX');
else
   XBX = blockX'*blockX; XBX = 0.5*(XBX+XBX');
end

[blockEig,lambda] = eig(XAX,XBX);

% Update blockX, blockAX and blockBX
blockX = blockX*blockEig;
blockAX = blockAX*blockEig;
if B_flag
   blockBX = blockBX*blockEig;
end

% Compute final residuals
if B_flag
   blockR = blockAX - blockBX*lambda;
else
   blockR = blockAX - blockX*lambda;
end

% Compute residual norms
resnorms = sqrt(diag(blockR'*blockR))./diag(lambda);

% Store lambda history and residual norms
resnorm_vec(iter,:) = resnorms;
lambda_vec(iter,:) = diag(lambda);
lambda = diag(lambda);

% Resize resnorm and lambda history
resnorm_vec = resnorm_vec(1:iter,:);
lambda_vec = lambda_vec(1:iter,:);

niter = iter - 1;

end
