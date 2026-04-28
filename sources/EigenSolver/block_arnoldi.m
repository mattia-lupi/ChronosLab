function [iter, D, X, res_norm_X, res_norm_D, flag] = block_arnoldi(A, nev, P, V0, target, max_iter, tol, verbose)
    if nargin < 3 || isempty(P), P = []; end
    if nargin < 4 || isempty(V0), V0 = randn(size(A,1), nev); end
    if nargin < 5 || isempty(target), target = 'largest'; end
    if nargin < 6 || isempty(max_iter), max_iter = 100; end
    if nargin < 7 || isempty(tol), tol = 1e-6; end
    if nargin < 8 || isempty(verbose), verbose = false; end

    n = size(V0, 1);
    r = size(V0, 2);

    if r < nev
        error('Starting block size must be >= nev.');
    end

    if isa(A, 'function_handle')
        opA = A;
        is_matrix_A = false;
    else
        opA = @(x) A * x;
        is_matrix_A = true;
    end

    if isempty(P)
        opP = @(x) x;
        has_P = false;
    elseif isa(P, 'function_handle')
        opP = P;
        has_P = true;
    else
        opP = @(x) P \ x;
        has_P = true;
    end

    V = zeros(n, (max_iter + 1) * r);
    H = zeros((max_iter + 1) * r, max_iter * r);

    [V(:, 1:r), ~] = qr(V0, 0);
    flag = 1; 

    is_smallest = strcmpi(target, 'lowest') || strcmpi(target, 'smallest');

    for iter = 1:max_iter
        idx = (iter-1)*r + 1 : iter*r;
        V_curr = V(:, idx);
        
        if is_smallest
            W = zeros(n, r);
            if is_matrix_A
                W = A \ V_curr; 
            else
                for col = 1:r
                    [W(:, col), ~] = gmres(A, V_curr(:, col), [], tol*1e-2, min(n, 100), P);
                end
            end
        else
            W = opA(V_curr);
            W = opP(W);
        end

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

        [~, sort_idx] = sort(abs(diagD), 'descend');
        
        k = min(nev, curr_size);
        sel_idx = sort_idx(1:k);
        
        if is_smallest
            D = 1 ./ diagD(sel_idx); 
        else
            D = diagD(sel_idx);
        end
        
        Y_k = Y(:, sel_idx);
        X = V(:, 1:curr_size) * Y_k;

        res_norm_X = zeros(k, 1);
        res_norm_D = zeros(k, 1);

        for idx_k = 1:k
            x = X(:, idx_k);
            lam = D(idx_k);
            
            if is_smallest || ~has_P
                res_norm_X(idx_k) = norm(opA(x) - lam * x);
            else
                res_norm_X(idx_k) = norm(opP(opA(x)) - lam * x);
            end
            
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