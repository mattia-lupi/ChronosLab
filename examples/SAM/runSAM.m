clc;
close all;

chronos_root = getappdata(0, 'chronos_root');

Name = "stickSlipFast";
files = dir(fullfile(chronos_root, 'examples', 'matrices', Name, strcat(Name, '_*.mat')));
fileNames = {files.name};

lagrange = true;

% Extract trailing numbers to enforce numeric sort (prevents 1, 10, 2 ordering)
tokens = regexp(fileNames, Name + '_(\d+)\.mat', 'tokens', 'once');
nums = cellfun(@(x) str2double(x{1}), tokens);
[sortedNums, sortIdx] = sort(nums);
fileNames = fileNames(sortIdx);

% Filter by positional index (e.g., 10th file in sequence onwards)
startVal = 1; 
fileNames = fileNames(startVal:end);
sortedNums = sortedNums(startVal:end);

sizeSeq = length(fileNames);

A   = cell(sizeSeq, 1);
TV0 = cell(sizeSeq, 1);
b   = cell(sizeSeq, 1);

for i = 1:sizeSeq
    fname = fullfile("..","matrices", Name, fileNames{i});
    data = load(fname);
    
    A{i} = data.Amat;
    b{i} = data.b;
    TV0{i} = data.TV0;
end

%% =========================================================================
%% Precompute Preconditioners Outside of the Cycle Once and for All
%% =========================================================================
M_all = cell(sizeSeq, 1);
t_AMG = zeros(sizeSeq, 1);

% Read the parameters in the input files in this directory
param = readDefaultParams();

for i = 1:sizeSeq
    [M_all{i}, t_AMG(i)] = buildPrecond(A{i}, TV0{i}, param, lagrange);
    disp(i)
end

%% =========================================================================
%% Solver Parameters & State Initialization (Set AFTER AMG Computation)
%% =========================================================================
gmres_restart = 100;
gmres_tol     = 1e-6;
gmres_maxit   = 10;
nstep         = 5;
step_size     = 1;
epss          = 1e-4;

% Strategy Columns: 1 = Recompute, 2 = Reuse (Recycle), 3 = SAM Adaptive
num_strategies = 3;
iter    = zeros(sizeSeq, num_strategies);
t_setup = zeros(sizeSeq, num_strategies);
t_solve = zeros(sizeSeq, num_strategies);

% Pre-assign setup times measured during AMG precomputation for Strategy 1
t_setup(:, 1) = t_AMG(1:sizeSeq);

% Step schedule configuration:
q  = [2 3];     % Strategy 2: Steps to recompute preconditioner
tP = [3];     % Strategy 3: Steps to recompute base preconditioner
tS = [2];     % Strategy 3: Steps to recompute SAM matrix

% State variables for recycling across the sequence cycle
P_reuse       = []; % Active preconditioner for Strategy 2
P_sam         = []; % Active base preconditioner for Strategy 3
idx_prec_base = 1;  % Index of matrix onto which P_sam was computed
sm            = []; % Active SAM matrix for Strategy 3

%% =========================================================================
%% Main Solution Cycle
%% =========================================================================
for i = 1:sizeSeq
   
   % --- Strategy 1: Recompute Preconditioner Every Step ---
   % Retrieve precomputed preconditioner
   M_recomp = M_all{i};
   
   t0 = tic;
   [~, ~, ~, it_recomp] = gmres_RIGHT(@(v) A{i}*v, b{i}, gmres_restart, gmres_tol, gmres_maxit, M_recomp, []);
   t_solve(i, 1) = toc(t0);
   iter(i, 1)    = (it_recomp(1)-1)*gmres_restart + it_recomp(2);
   
   
   % --- Strategy 2: Reuse/Recycle Preconditioner ---
   if i == 1 || ismember(i, q)
      % Reuse precomputed preconditioner for the current step and record AMG setup time
      P_reuse = M_recomp;
      t_setup(i, 2) = t_setup(i, 1);
   else
      % Reuse existing active preconditioner without setup overhead
      t_setup(i, 2) = 0;
   end
   
   Pfun_reuse = @(v) sam_apply_left(speye(size(A{i}, 1)), P_reuse, v);
   
   t0 = tic;
   [~, ~, ~, it_reuse] = gmres_RIGHT(@(v) A{i}*v, b{i}, gmres_restart, gmres_tol, gmres_maxit, Pfun_reuse, []);
   t_solve(i, 2) = toc(t0);
   iter(i, 2)    = (it_reuse(1)-1)*gmres_restart + it_reuse(2);
   
   
   if ~isempty(tS)
      % --- Strategy 3: SAM Adaptive ---
      % 1. Update Base Preconditioner if at start or at designated step (tP)
      if i == 1 || ismember(i, tP)
         P_sam = M_recomp;
         t_prec_sam = t_setup(i, 1);
         t_setup(i, 3) = t_setup(i, 3) + t_prec_sam;
      
         idx_prec_base = i;  % Track matrix index used for base preconditioner
         sm = speye(size(A{i}, 1)); % Reset SAM to identity since matrix matches
      end
      
      % 2. Update SAM if at designated step (tS)
      if ismember(i, tS)
         if i == idx_prec_base
            sm = speye(size(A{i}, 1));
         else
            t0_sam = tic;
            % Compute SAM between current matrix A{i} and base matrix A{idx_prec_base}
            [sm, normm] = MEX_sam_adaptive_left(A{i}, A{idx_prec_base}, 8, nstep, step_size, epss);
            t_setup(i, 3) = t_setup(i, 3) + toc(t0_sam);
         end
      end
      
      % 3. Apply SAM + Base Preconditioner
      Pfun_sam = @(v) sam_apply_left(sm, P_sam, v);
      
      t0 = tic;
      [~, ~, ~, it_sam] = gmres_RIGHT(@(v) A{i}*v, b{i}, gmres_restart, gmres_tol, gmres_maxit, Pfun_sam, []);
      t_solve(i, 3) = toc(t0);
      iter(i, 3)    = (it_sam(1)-1)*gmres_restart + it_sam(2);
   end
end

%% Visualization of Performance Results
figure('Name', "SAM Benchmark - " + Name, 'Color', 'w', 'Position', [100, 100, 1200, 800]);
labels = {'Recompute Precond', 'Reuse Precond', 'SAM Adaptive'};
colors = [0.8500 0.3250 0.0980;   % Red/Orange
          0.9290 0.6940 0.1250;   % Yellow/Gold
          0.4660 0.6740 0.1880];  % Green

steps = 1:sizeSeq;
tlayout = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tlayout, "Performance Analysis: Dataset [" + Name + "]", 'FontSize', 14, 'FontWeight', 'bold');

% --- Panel 1: Number of GMRES Iterations ---
nexttile;
hold on;
plot(steps, iter(:, 1), '-o', 'LineWidth', 1.8, 'Color', colors(1, :), 'MarkerSize', 3, 'MarkerFaceColor', colors(1, :));
plot(steps, iter(:, 2), '-o', 'LineWidth', 1.8, 'Color', colors(2, :), 'MarkerSize', 3, 'MarkerFaceColor', colors(2, :));
plot(steps, iter(:, 3), '-o', 'LineWidth', 1.8, 'Color', colors(3, :), 'MarkerSize', 3, 'MarkerFaceColor', colors(3, :));

if exist('q', 'var') && ~isempty(q)
   xline(q, ':', 'Color', colors(2, :), 'LineWidth', 1.5);
   labels{end+1} = 'Recomp. Precond (Reuse)';
end
if exist('tP', 'var') && ~isempty(tP)
   xline(tP, ':', 'Color', colors(3, :), 'LineWidth', 1.5);
   labels{end+1} = 'Recomp. Precond (SAM)';
end
if exist('tS', 'var') && ~isempty(tS)
   xline(tS, '--', 'Color', colors(3, :), 'LineWidth', 1.5);
   labels{end+1} = 'Compute SAM';
end
hold off;
grid on;
title('GMRES Iterations');
xlabel('Sequence Step'); ylabel('Iterations');
legend(labels, 'Location', 'best');

% --- Panel 2: Setup Time (Preconditioner + SAM) ---
nexttile;
b_setup = bar(steps, t_setup, 'grouped');
for k = 1:3
   b_setup(k).FaceColor = colors(k, :);
end
grid on;
title('Setup Time (Preconditioner + SAM)');
xlabel('Sequence Step'); ylabel('Time (s)');
legend(labels(1:3), 'Location', 'best');

% --- Panel 3: Time for Solving System ---
nexttile;
hold on;
for k = 1:3
   plot(steps, t_solve(:, k), '-s', 'LineWidth', 1.8, 'Color', colors(k, :), ...
        'MarkerSize', 3, 'MarkerFaceColor', colors(k, :));
end
hold off;
grid on;
title('GMRES Solve Time');
xlabel('Sequence Step'); ylabel('Time (s)');
legend(labels, 'Location', 'best');

% --- Panel 4: Cumulative Total Time (Setup + Solve) ---
nexttile;
start = 1;
t_total = cumsum(t_setup + t_solve, 1);
hold on;
for k = start:3
   plot(steps, t_total(:, k), '-^', 'LineWidth', 2.0, 'Color', colors(k, :), ...
        'MarkerSize', 3, 'MarkerFaceColor', colors(k, :));
end
hold off;
grid on;
title('Cumulative Runtime (Setup + Solve)');
xlabel('Sequence Step'); ylabel('Total Accumulated Time (s)');
legend(labels(start:3), 'Location', 'northwest');


%% Helper Functions
function [M, t_prec] = buildPrecond(A_mat, TV0_vec, param, lagrange)
    symm = (norm(A_mat-A_mat','f')/norm(A_mat,'f') < 1e-14);
    param.symm = symm;
    t0 = tic;
    if ~lagrange
       [M, ~] = computeAMG(A_mat, TV0_vec, param, 1);
    else
       [M, ~] = computeRACP(A_mat, TV0_vec, param, 1);
    end
    t_prec = toc(t0);
end
