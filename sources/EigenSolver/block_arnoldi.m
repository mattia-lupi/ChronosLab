function [iter, D, X, res_norm_X, res_norm_D, flag] = block_arnoldi(A, nev, P, V0, largest_flag, max_iter, tol, verbose)
    if nargin < 3 || isempty(P), P = {[], []}; end
    if nargin < 4 || isempty(V0), V0 = randn(size(A,1), nev); end
    if nargin < 5 || isempty(largest_flag), largest_flag = true; end
    if nargin < 6 || isempty(max_iter), max_iter = 100; end
    if nargin < 7 || isempty(tol), tol = 1e-6; end
    if nargin < 8 || isempty(verbose), verbose = false; end
    
    if largest_flag
       target = 'largest';
    else
       target = 'smallest';
    end
    
    n = size(V0, 1);
    r = size(V0, 2);
    
    if r < nev
        error('Starting block size must be >= nev.');
    end
    
    % Setup matrix-vector operator
    if isa(A, 'function_handle')
        opA = A;
    else
        opA = @(x) A * x;
    end
    
    % Setup matrix-free preconditioner application
    if isa(P, 'function_handle')
        applyP = P;
    elseif iscell(P)
        if numel(P) >= 2 && ~isempty(P{1}) && ~isempty(P{2})
            P1 = P{1};
            P2 = P{2};
            if isa(P1, 'function_handle') && isa(P2, 'function_handle')
                applyP = @(x) P2(P1(x));
            elseif isa(P1, 'function_handle')
                applyP = @(x) P2 * P1(x);
            elseif isa(P2, 'function_handle')
                applyP = @(x) P2(P1 * x);
            else
                applyP = @(x) P2 * (P1 * x);
            end
        elseif numel(P) >= 1 && ~isempty(P{1})
            P1 = P{1};
            if isa(P1, 'function_handle')
                applyP = P1;
            else
                applyP = @(x) P1 * x;
            end
        else
            applyP = @(x) x;
        end
    elseif ~isempty(P)
        applyP = @(x) P * x;
    else
        applyP = @(x) x;
    end
    
    V = zeros(n, (max_iter + 1) * r);
    H = zeros((max_iter + 1) * r, max_iter * r);
    [V(:, 1:r), ~] = qr(V0, 0);
    flag = 1;
    
    is_smallest = strcmpi(target, 'lowest') || strcmpi(target, 'smallest');
    
    for iter = 1:max_iter
        idx = (iter-1)*r + 1 : iter*r;
        V_curr = V(:, idx);
        
        % Matrix-free block Krylov step on preconditioned operator M^-1 * A
        Av = opA(V_curr);
        W = applyP(Av);
        
        % Gram-Schmidt orthogonalization against previous blocks
        for i = 1:iter
            idx_i = (i-1)*r + 1 : i*r;
            Vi = V(:, idx_i);
            Hij = Vi' * W;
            H(idx_i, idx) = Hij;
            W = W - Vi * Hij;
        end
        
        [V_next, H_next] = qr(W, 0);
        H(iter*r + 1 : (iter+1)*r, idx) = H_next;
        V(:, iter*r + 1 : (iter+1)*r) = V_next;
        
        curr_size = iter * r;
        H_curr = H(1:curr_size, 1:curr_size);
        
        [Y, D_mat] = eig(H_curr);
        diagD = diag(D_mat);
        if is_smallest
            [~, sort_idx] = sort(abs(diagD), 'ascend');
        else
            [~, sort_idx] = sort(abs(diagD), 'descend');
        end
        
        k = min(nev, curr_size);
        sel_idx = sort_idx(1:k);
        
        D = diagD(sel_idx);
        Y_k = Y(:, sel_idx);
        X = V(:, 1:curr_size) * Y_k;
        
        res_norm_X = zeros(k, 1);
        res_norm_D = zeros(k, 1);
        
        for idx_k = 1:k
            x = X(:, idx_k);
            lam = D(idx_k);
            res_norm_X(idx_k) = norm(applyP(opA(x)) - lam * x);
            y_last = Y_k(end-r+1:end, idx_k);
            res_norm_D(idx_k) = norm(H_next * y_last);
        end
        
        max_res = max(res_norm_X);
        if verbose
            fprintf('Iter: %d, Max Res X: %e, Max Res D: %e\n', iter, max_res, max(res_norm_D));
        end
        
        if max_res < tol && curr_size >= nev
            flag = 0;
            break;
        end
    end
end