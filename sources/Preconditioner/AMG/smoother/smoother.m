function smootherOp = smoother(A, symm_flag, param, verb)

if nargin < 4
   verb = 1;
end
amg_times;

% Mandatory parameters
nthread         = param.nthread;
nstep           = param.nstep;
step_size       = param.step_size;
epsilon         = param.epsilon;
method          = param.method;

% Init the outer smoother to void
smootherOp.left_out = [];
smootherOp.right_out = [];

% Force use of non-symmetric FSAI if the problem is not symmetric
if ~symm_flag && ~strcmpi(method,'jacobi')
   method = 'afsai_nsy';
end

switch lower(method)

    case 'afsai_sym'
        % Set-up AFSAI (afsai with mex-cpp code)
        F = afsai_cpp(A,nthread,nstep,step_size,epsilon);
        % Correct NaN for very ill-conditioned problems
        irow = find(isnan(diag(F)));
        if numel(irow) > 0
           if verb
              fprintf('WARNING: Correcting %d diagonals\n',numel(irow));
           end
        end
        for i = 1:numel(irow)
           ii = irow(i);
           F(ii,:) = 0;
           F(ii,ii) = 1 / sqrt(A(ii,ii));
           A(ii,ii)
        end
        % Compute damping parameter
        FAFT = @(x) F*(A*(F'*x));
        opts.issym = 1;
        opts.disp = verb;
        opts.tol = 5.e-4;
        lambda = eigs(FAFT,size(A,1),1,'la',opts);
        if verb
           fprintf('Max Lambda: %10.4f\n',lambda);
        end
        omega = min(1,1.9 / lambda);
        % Append the smoother
        smootherOp.left = F;
        smootherOp.right = F';
        smootherOp.omega = omega;
        smootherOp.lambda = lambda;

    case 'afsai_nsy'
        if verb
           fprintf('Non-Symmetric AFSAI is used\n');
        end
        % Set-up AFSAI_NSY (afsai for nsy systems with mex-cpp code)
        [FL,FU] = NSY_rfsai_cpp(nstep,step_size,epsilon,A);
        % Compute damping parameter
        FAFT = @(x) FL*(A*(FU*x));
        opts.issym = 0;
        opts.disp = verb;
        opts.tol = 5.e-4;
        lambda = eigs(FAFT,size(A,1),1,'lm',opts);
        if verb
           fprintf('Max Lambda: %10.4f\n',lambda);
        end
        omega = min(1,1.9 / lambda);
        % Append the smoother
        smootherOp.left = FL;
        smootherOp.right = FU;
        smootherOp.omega = omega;
        smootherOp.lambda = lambda;

    case 'jacobi'
        % Compute Diagonal
        F = 1 ./ sqrt(full(diag(A)));
	     F = diag(sparse(F));
        % Compute damping parameter
        FAFT = @(x) F*(A*(F'*x));
        % opts.issym = 1;
        if true
           lambda = eigs(FAFT,size(A,1),1,'lm','IsFunctionSymmetric',1,...
                         'Tolerance',1.e-2,'Display',verb,'FailureTreatment','keep');
           if verb
              fprintf('Max Lambda: %10.4f\n',lambda);
           end
           omega = min(1,1.9 / lambda);
        else
           lambda = 2.0;
           omega = 1.0;
        end
        % Append the smoother
        smootherOp.left = F;
        smootherOp.right = F';
        smootherOp.omega = omega;
        smootherOp.lambda = lambda;

    otherwise
        error('Not existing method');

end

return
