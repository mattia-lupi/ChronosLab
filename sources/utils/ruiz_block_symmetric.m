function [Abal, D] = balance_block(A, maxit, tol, verb)
%BALANCE_BLOCK  Iterative diagonal balancing for an NxN block sparse matrix.
%   A    : NxN cell array of sparse blocks (some entries may be empty)
%   D    : Nx1 cell array of diagonal scaling matrices
%   Abal : NxN cell array of balanced blocks

    if nargin<2, maxit = 10;   end
    if nargin<3, tol   = 1e-2; end
    if nargin<4, verb  = 1;    end

    nb = size(A, 1);

    % --- Infer block sizes from first non-empty block in each row ---
    blksize = zeros(nb, 1);
    for i = 1:nb
        for j = 1:nb
            if ~isempty(A{i,j})
                blksize(i) = size(A{i,j}, 1);
                break;
            end
        end
        if blksize(i) == 0
            error('balance_block: cannot determine size of block-row %d (entire row is empty)', i);
        end
    end

    % --- Initialize ---
    D    = arrayfun(@(n) speye(n), blksize, 'UniformOutput', false);  % Nx1 cell
    Abal = A;   % COW: no copy until a block is modified

    if maxit < 1
        if verb, fprintf('Made 0 iterations\n'); end
        return;
    end

    for k = 1:maxit

        % Row norms: sum over all non-empty blocks in each block-row
        r = cell(nb, 1);
        for i = 1:nb
            ri = zeros(blksize(i), 1);
            for j = 1:nb
                if ~isempty(Abal{i,j})
                    ri = ri + sum(abs(Abal{i,j}).^2, 2);
                end
            end
            r{i} = max(sqrt(ri), 1e-15);
        end

        s    = cellfun(@(ri) 1./sqrt(ri), r, 'UniformOutput', false);
        conv = max(abs(vertcat(s{:}) - 1));

        if verb
            fprintf('ITER: %d | Max(abs(s-1)): %e\n', k, conv);
        end
        if conv < tol, break; end

        % Scaling matrices
        S = cellfun(@(si, n) spdiags(si, 0, n, n), s, num2cell(blksize), ...
                    'UniformOutput', false);

        % Accumulate D
        for i = 1:nb
            D{i} = D{i} * S{i};
        end

        % Scale non-empty blocks: S{i} * Abal{i,j} * S{j}
        for i = 1:nb
            for j = 1:nb
                if ~isempty(Abal{i,j})
                    Abal{i,j} = S{i} * Abal{i,j} * S{j};
                end
            end
        end

    end

    if verb
        fprintf('Made %d iterations\n', k);
        fprintf('Error %e\n', conv);
    end
end
