clear
clc

tWall = tic();      % start wall-clock timer
tCPU  = cputime;    % start CPU timer

testFiles = {
   % Test DOM + BAMG + prolongation smoothing + symFSAI
   fullfile('CubeBAMG_NONE.m')
   fullfile('CubeBAMG_SMOOTH.m')
   fullfile('CubeBAMG_EMIN.m')

   % Test CLA + EXTI + Aggressive coarsening (and not)  + Jacobi +
   % % SRQCG/LOBPCG/Lanczos
   % fullfile('FlowCLAS_0.m')
   % fullfile('FlowCLAS_AGG.m')
   % fullfile('FlowCLAS_AGG_SRQCG.m')
   % fullfile('FlowCLAS_AGG_LOBPCG.m')
   % fullfile('FlowCLAS_AGG_LANCZOS.m')

   % Test DOM(-) + BAMG + nsyFSAI + arnoldi + filterProl/Oper
   fullfile('dpNsyFSAI_NONE.m')
   fullfile('dpNsyFSAI_ARNOLDI.m')
   fullfile('dpNsyFSAI_FILTERpro.m')
   fullfile('dpNsyFSAI_FILTERoper.m')

   % Test nsyRACP + nsyFsai + DOM + BAMG + SMOOTH
   fullfile('stickSlipBAMG_SMOOTH.m')

};

results = runtests(testFiles);

elapsedWall = toc(tWall);
elapsedCPU  = cputime - tCPU;

fprintf("Elapsed wall-clock time: %1.2f s\n", elapsedWall);
fprintf("Elapsed CPU time:        %1.2f s\n", elapsedCPU);

if any([results.Failed])
   error("Some tests did not pass");
else
   disp("All tests passed");
end
