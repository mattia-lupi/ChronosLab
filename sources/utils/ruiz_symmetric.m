function [Abal,D] = ruiz_symmetric(A,maxit,tol)
    if nargin<2, maxit=10; end
    if nargin<3, tol=1e-2; end
    n = size(A,1);
    D = speye(n);
    Abal = A;
    if maxit < 1
       return;
    end
    for k = 1:maxit
        r = sqrt(sum(abs(Abal).^2,2));    % 2-norm; usare max(abs(...),[],2) per ∞-norm
        r = max(r,1e-15);                 % evita divisione per zero
        s = 1./sqrt(r);
        fprintf('ITER: %d | ',k);
        fprintf('Max(abs(s-1)): %e\n',max(abs(s-1)));
        if max(abs(s-1)) < tol, break; end
        S = spdiags(s,0,n,n);
        D = D * S;
        Abal = S * Abal * S;
        %row = full(vecnorm(Abal'));
        %col = full(vecnorm(Abal));
        %fprintf('ROW: mean %e devStd %e max %e min %e\n',mean(row),std(row),max(row),min(row));;
        %fprintf('COL: mean %e devStd %e max %e min %e\n',mean(col),std(col),max(col),min(col));
    end
    fprintf('Fatte %d iterazioni\n',k);
    fprintf('Errore %e\n',max(abs(s-1)));
end
