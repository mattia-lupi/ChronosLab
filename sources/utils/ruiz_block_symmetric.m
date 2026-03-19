function [Abal,D] = ruiz_block_symmetric(A,maxit,tol,verb)
%ruiz_block_symmetric  Iterative diagonal balancing for a 2x2 block sparse matrix.
%   A    : 2x2 cell array { A00, A01; A10, A11 }
%   D    : 2x1 cell array of diagonal scaling matrices {D0; D1}
%   Abal : 2x2 cell array of balanced blocks

    if nargin<2, maxit = 10;   end
    if nargin<3, tol   = 1e-2; end
    if nargin<4, verb  = 1;    end

    % Unpack
    A00 = A{1,1};  A01 = A{1,2};
    A10 = A{2,1};  A11 = A{2,2};

    n0 = size(A00, 1);
    n1 = size(A11, 1);

    D0 = speye(n0);
    D1 = speye(n1);

    for k = 1:maxit

        % Row norms spanning both blocks in each block-row
        r0 = sqrt( sum(abs(A00).^2, 2) + sum(abs(A01).^2, 2) );
        r1 = sqrt( sum(abs(A10).^2, 2) + sum(abs(A11).^2, 2) );

        r0 = max(r0, 1e-15);
        r1 = max(r1, 1e-15);

        s0 = 1./sqrt(r0);
        s1 = 1./sqrt(r1);

        conv = max(abs([s0; s1] - 1));
        if verb
            fprintf('ITER: %d | Max(abs(s-1)): %e\n', k, conv);
        end
        if conv < tol, break; end

        S0 = spdiags(s0, 0, n0, n0);
        S1 = spdiags(s1, 0, n1, n1);

        D0 = D0 * S0;
        D1 = D1 * S1;

        A00 = S0 * A00 * S0;
        A01 = S0 * A01 * S1;
        A10 = S1 * A10 * S0;
        A11 = S1 * A11 * S1;

    end

    if verb
        fprintf('Made %d iterations\n', k);
        fprintf('Error %e\n', conv);
    end

    % Repack
    D    = {D0; D1};
    Abal = {A00, A01; A10, A11};

end
