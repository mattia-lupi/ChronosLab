function [niter, lambda, blockX, resnorm_vec, lambda_vec, flag] = Lanczos(...
    A, blockX0, itmax, tol, largest_flag, verb)
% LANCZOS Partial Reorthogonalized Block Lanczos Eigensolver
%
% Inputs:
%   A            - Matrix (n x n) or function handle @(x) A*x
%   blockX0      - Starting block of vectors (n x p)
%   itmax        - Maximum number of block Lanczos iterations
%   tol          - Convergence tolerance on residual bounds (default: 1e-6)
%   largest_flag - true for largest algebraic, false for smallest (default: false)
%   verb         - true to print convergence progress (default: true)
%
% Outputs:
%   niter        - Number of iterations completed
%   lambda       - Converged Ritz values (p x 1)
%   blockX       - Converged orthonormal Ritz vectors (n x p)
%   resnorm_vec  - Final residual norms for all p eigenvalues (p x 1)
%   lambda_vec   - History of residual norms per iteration (niter x p)
%   flag         - 0 if all p eigenvalues converged within tolerance, 1 otherwise

    % 1. Argument Handling and Operator Setup
    if nargin < 4 || isempty(tol), tol = 1e-6; end
    if nargin < 5 || isempty(largest_flag), largest_flag = false; end
    if nargin < 6 || isempty(verb), verb = true; end

    if isa(A, 'function_handle')
        applyA = A;
    else
        applyA = @(x) A * x;
    end

    [n, p] = size(blockX0);
    itmax  = min(itmax, max(1, floor(n / p)));

    eps_mach = eps;
    sqrteps  = sqrt(eps_mach);
    eta      = eps_mach^(3/4);

    % 2. Initialization and Basis Preallocation
    V_cell = cell(itmax + 1, 1);
    A_diag = cell(itmax, 1);
    B_sub  = cell(itmax, 1);

    % Initial orthonormal basis
    [V_cell{1}, ~] = qr(blockX0, 0);

    % Orthogonality tracking matrix: omega(i, j) ~ ||V_i' * V_j||_2
    omega = eye(itmax + 1);

    lambda_vec  = zeros(itmax, p);
    resnorm_vec = zeros(p, 1);
    flag        = 1;
    niter       = 0;

    if verb
        fprintf('\n%5s | %12s | %18s\n', 'Iter', 'Not Conveged', 'Max Res Norm');
        fprintf('%s\n', repmat('-', 1, 55));
    end

    % 3. Main Lanczos Iterations
    for j = 1:itmax
        niter = j;

        % Apply matrix operator: W = A * V_j
        W = applyA(V_cell{j});

        % Subtract subdiagonal projection: W = W - V_{j-1} * B_{j-1}'
        if j > 1
            W = W - V_cell{j-1} * (B_sub{j-1}');
        end

        % Project diagonal block: A_j = V_j' * W
        A_j = V_cell{j}' * W;
        A_diag{j} = (A_j + A_j') / 2; % Enforce symmetry
        W = W - V_cell{j} * A_diag{j};

        % Generate next basis block via QR: W = V_{j+1} * B_j
        [V_next, B_j] = qr(W, 0);
        B_sub{j} = B_j;
        norm_Bj  = norm(B_j, 2);

        % 4. Partial Reorthogonalization (PRO) Recurrence
        omega(j, j+1) = eps_mach;
        omega(j+1, j) = eps_mach;

        if j > 1 && norm_Bj > eps_mach
            for i = 1:(j - 1)
                t_prev = 0;
                if i > 1
                    t_prev = norm(B_sub{i-1}, 2) * omega(i-1, j);
                end

                if i == j - 1
                    % Theoretical cancellation: omega(j,j) - omega(j-1,j-1) ~ eps_mach
                    t_diff = norm(B_sub{j-1}, 2) * eps_mach;
                    t_next = 0;
                else
                    t_diff = norm(B_sub{j-1}, 2) * omega(i, j-1);
                    t_next = norm(B_sub{i}, 2) * omega(i+1, j);
                end

                num = t_diff + norm(A_diag{j} - A_diag{i}, 2) * omega(i, j) + ...
                      t_prev + t_next + 2 * eps_mach * (norm(A_diag{j}, 2) + norm(A_diag{i}, 2));

                omega(i, j+1) = num / norm_Bj;
                omega(j+1, i) = omega(i, j+1);
            end

            % Selective Reorthogonalization (Vectorized DGKS)
            if max(omega(1:j-1, j+1)) > sqrteps
                reorth_idx = reshape(find(omega(1:j, j+1) > eta), 1, []);
                if isempty(reorth_idx)
                    reorth_idx = 1:j;
                end

                V_sub = [V_cell{reorth_idx}];
                for pass = 1:2
                    V_next = V_next - V_sub * (V_sub' * V_next);
                end

                [V_next, R_corr] = qr(V_next, 0);
                B_sub{j} = R_corr * B_sub{j};

                omega(reorth_idx, j+1) = eps_mach;
                omega(j+1, reorth_idx) = eps_mach;
            end
        end

        if j < itmax
            V_cell{j+1} = V_next;
        end

        % 5. Rayleigh-Ritz Projection
        Tk = zeros(j * p, j * p);
        for k = 1:j
            idx_k = (k-1)*p + 1 : k*p;
            Tk(idx_k, idx_k) = A_diag{k};
            if k < j
                idx_next = k*p + 1 : (k+1)*p;
                Tk(idx_next, idx_k) = B_sub{k};
                Tk(idx_k, idx_next) = B_sub{k}';
            end
        end

        [S, D] = eig((Tk + Tk') / 2);
        theta  = diag(D);

        if largest_flag
            [theta_sorted, sort_perm] = sort(theta, 'descend');
        else
            [theta_sorted, sort_perm] = sort(theta, 'ascend');
        end

        S_sorted = S(:, sort_perm);
        lambda   = theta_sorted(1:p);

        % Compute residual bounds: ||A*x_i - lambda_i*x_i||_2 <= ||B_j * S_bottom(:, i)||_2
        S_bottom = S_sorted(end-p+1:end, 1:p);
        for ip = 1:p
            resnorm_vec(ip) = norm(B_sub{j} * S_bottom(:, ip));
        end
        lambda_vec(j, :) = resnorm_vec';

        % Compute convergence counts
        num_converged = sum(resnorm_vec <= tol);
        num_remaining = p - num_converged;

        if verb && (mod(j,10) == 0 || j == 1)
            fprintf('%5d | %5d / %-2d | %18.10e\n', ...
                j, num_remaining, p, max(resnorm_vec));
        end

        % Exit when all p eigenvalues meet tolerance or invariant subspace found
        if num_converged == p || norm_Bj <= eps_mach
            flag = 0;
            break;
        end
    end

    % 6. Assemble Ritz Vectors
    V_total = [V_cell{1:niter}];
    blockX  = V_total * S_sorted(:, 1:p);

    % Normalize eigenvectors
    for ip = 1:p
        blockX(:, ip) = blockX(:, ip) / norm(blockX(:, ip));
    end

    lambda_vec = lambda_vec(1:niter, :);
end