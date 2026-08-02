function [N,avg_resnorm] = sam_adaptive_left(A, A0, nstep, step_size, eps)

% Allocate 
N = speye(size(A,1));

mHatkOld = N(:,1);
avg_resnorm = 0;

% Precompute outside the main loop
col_norms_sq = full(sum(A.^2, 1));

for k = 1:size(A,1)
   % Diagonal initialization
   J = k;

   % Reset n2
   n2 = 0;

   A0kk = A0(:,k);
   normA0k = norm(A0kk);

   for t = 1:nstep
      % Get nonzero indices in all columns J of A
      [I, ~] = find(any(A(:,J),2));
      
      % Get sizes of the LS problem
      n2old = n2;
      n2 = size(J,2);

      % First iteration, solve full qr
      if t == 1
         % Get contracted initial matrix for column k
         a0k = full(A0(I,k));

         % Get contracted matrix A
         AHat = full(A(I,J));
      
         % Compute qr factorization of AHat
         [Q,R] = qr(AHat);
      
         % Make R to be the upper triangular it should be and discard the zeros
         % under
         R = R(1:n2,:);

         % Save I for the next iter
         Iold = I;
      else
         % Get the new columns
         Itilde = setdiff(I,Iold);
         n2tilde = length(J_add);

         % Get contracted initial matrix for column k
         a0k = [a0k; A0(Itilde,k)];

         % Find matrix 1,2 in equation 14 in the paper
         A_I_Jadd = full(A(Iold,J_add));

         % Save I for the next iter
         Iold = [Iold; Itilde];

         % Find matrix 2,2 in equation 14 in the paper
         A_Itilde_Jadd = full(A(Itilde,J_add));

         % Apply existing orthogonal transformations
         QAhat12 = Q' * A_I_Jadd;
         
         % Exact mathematical partitioning matching Equation (15)
         B1 = QAhat12(1:n2old, :);
         B2_top = QAhat12(n2old+1:end, :);
         B2 = [B2_top; A_Itilde_Jadd];

         % Compute QR decomposition of the small sub-block B2
         [Qtilde, Rtilde] = qr(B2);
         Rtilde = Rtilde(1:n2tilde,:);

         % Update Q and R
         R = [R B1; zeros(n2tilde,n2old) Rtilde];
         Q_big = blkdiag(Q, eye(length(Itilde)));
         Update_block = blkdiag(eye(n2old), Qtilde);
         Q = Q_big * Update_block;
      end

      % Compute cHat
      cHat = Q'*a0k;
      
      % Solve the least square problem and compute the values in mHatk
      mHatk = R\(cHat(1:n2));

      % Compute residual
      AjMh = A(:,J)*mHatk;
      res = AjMh - A0kk;

      % Save mHatk
      mHatkOld = mHatk;

      % Check exit residual
      res_norm = 2*norm(res)/(norm(AjMh)+normA0k);
      if res_norm < eps
         break;
      end
   
      if t < nstep
         % Get nonzero indices for res
         L = abs(res) > 0; 
      
         % Find potential new column indices
         [Jtilde, ~] = find(A(:,L));
         Jtilde = unique(Jtilde);
      
         % Take away the indices that already are in J
         Jtilde = setdiff(Jtilde,J);
   
         if isempty(Jtilde)
            break;
         end
      
         % Compute the residual improvement given by index j
         rTA = res' * A(:, Jtilde);
         sum_A2 = col_norms_sq(Jtilde);
         
         rho_j2 = norm(res, 2)^2 - (rTA.^2 ./ sum_A2); 
      
         % Select the minimum
         if step_size == 1
            [~,Jnew] = min(rho_j2);
      
            % Get the Jtilda index to add
            J_add = Jtilde(Jnew);
            J = [J J_add];
         else
            [~, idx] = sort(rho_j2);
   
            % Get the Jtilda indexes to add
            J_add = Jtilde(idx(1:step_size))';
            J = [J J_add]; 
         end
      end
   end

   % fprintf("avg %e t %d\n",res_norm,t);
   avg_resnorm = avg_resnorm + res_norm;

   % Update SPAI matrix
   N(J,k) = mHatkOld;

   clear J;
   
end
avg_resnorm = avg_resnorm / size(A,1);

end
