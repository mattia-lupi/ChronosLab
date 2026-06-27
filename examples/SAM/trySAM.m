
clc;
clear;

Name = "richardsBig";%richards richardsBig stickSlipFast StickSlipOpenBig nicolas
files = dir('mats/'+Name+'/'+Name+'_*.mat');
fileNames = {files.name};
if contains(Name,"richards") || contains(Name,"dp")
   lagrange = false;
else
   lagrange = true;
end


% Extract trailing numbers to enforce numeric sort (prevents 1, 10, 2 ordering)
tokens = regexp(fileNames, Name+'_(\d+)\.mat', 'tokens', 'once');
nums = cellfun(@(x) str2double(x{1}), tokens);
[~, sortIdx] = sort(nums);
fileNames = fileNames(sortIdx);

sizeSeq = length(fileNames);

if strcmp(Name,"nicolas")
   sizeSeq = 22;
end

A = cell(sizeSeq, 1);
TV0 = cell(sizeSeq, 1);
b = cell(sizeSeq, 1);

for i = 1:sizeSeq
    % Load into struct to prevent namespace collision
    fname = "mats/"+Name+"/"+fileNames{i};
    data = load(fname);
    
    % Replace 'mat_name' and 'rhs_name' with actual internal variable names
    if ~strcmp(Name,"nicolas")
       A{i} = data.Amat;
    else
       A{i} = data.A;
    end
    if ~strcmp(Name,"nicolas")
       b{i} = data.b;
    end

    if ~contains(Name,"richards")
       if strcmp(Name,"StickSlipOpenBig")
          TV0{i} = data.TV;
       else
          TV0{i} = data.TV0;
       end
    else
       TV0{i} = ones(size(A{i},1),1);
    end
end




%% See how much they differ in norm

normm = zeros(sizeSeq,sizeSeq);

for i = 1:sizeSeq
   for j = 1:i-1
      normm(i,j) = 2 * norm(A{j}-A{i},'f')/(norm(A{j},'f')+norm(A{i},'f'));
   end
end

%% See if pattern has changed

normPatt = zeros(sizeSeq,sizeSeq);
nnzrA = zeros(sizeSeq,1);

for i = 1:sizeSeq
   nnzrA(i) = nnz(A{i})/size(A{i},1);
   for j = 1:i-1
      normPatt(i,j) = norm(logical(A{j})-logical(A{i}),'f');
   end
end

%% Compute Sam

if ~contains(Name,"Slip") && ~contains(Name,"dp")
   normmSam1 = zeros(sizeSeq,1);

   for i = 2:sizeSeq
      tic
      pre_mex1 = MEX_sam_preprocess_left(A{i});
      time_prepr(i,1) = toc;
      tic
   
      tic
      sam1{i} = MEX_sam_compute_left(A{i},A{1},pre_mex1);
      time_compute1(i) = toc;
      nnzSam(i,1) = nnz(sam1{i});
      normmSam1(i) = 2 * norm(A{1}-sam1{i}*A{i},'f')/(norm(A{1},'f')+norm(sam1{i}*A{i},'f'));
   end
else

   normmSam1 = zeros(sizeSeq,sizeSeq);
   normmSam2 = zeros(sizeSeq,sizeSeq);
   nnzSam1 = zeros(sizeSeq,sizeSeq);
   nnzSam2 = zeros(sizeSeq,sizeSeq);
   time_compute1 = zeros(sizeSeq,sizeSeq);
   time_compute2 = zeros(sizeSeq,sizeSeq);
   sam1 = cell(sizeSeq,sizeSeq);
   sam2 = cell(sizeSeq,sizeSeq);
   
   for i = 1:sizeSeq
      tic
      pre_mex1 = MEX_sam_preprocess_left(A{i});
      time_prepr(i,1) = toc;
   
      % tic
      % nn = size(TV0{i},1);
      % pre_mex2 = MEX_sam_preprocess_left(A{i},[speye(nn), A{i}(1:nn,nn+1:end); A{i}(nn+1:end,1:nn) A{i}(nn+1:end,nn+1:end)]);
      % time_prepr(i,2) = toc;
   
      for j = 1:i-1   
         tic
         sam1{i,j} = MEX_sam_compute_left(A{i},A{j},pre_mex1);
         time_compute1(i,j) = toc;
         nnzSam1(i,j) = nnz(sam1{i,j});
         normmSam1(i,j) = 2 * norm(A{j}-sam1{i,j}*A{i},'f')/(norm(A{j},'f')+norm(sam1{i,j}*A{i},'f'));
   
         % tic
         % sam2{i,j} = MEX_sam_compute_left(A{i},A{j},pre_mex2);
         % time_compute2(i,j) = toc;
         % nnzSam2(i,j) = nnz(sam2{i,j});
         % normmSam2(i,j) = 2 * norm(A{j}-sam2{i,j}*A{i},'f')/(norm(A{j},'f')+norm(sam2{i,j}*A{i},'f'));
      end
   end
end

%% Compute the Preconditioner

gmres_restart = 100;
gmres_tol     = 1e-6;
gmres_maxit   = 20;

if ~lagrange
   [M,time(1)] = computeAMG(A{1},TV0{1},false);
else
   [M,time(1)] = computeRACP(A{1},TV0{1},false);
end

t_amg0    = time(1);
P0_amg    = M;
fprintf('done in %.3f s\n\n', t_amg0);

%% Solve the system

% for i = 1:sizeSeq
% 
%    uni = ones(size(A{i},1),1);
% 
%    if i > 1
%       if ~lagrange
%          [MM{i},time(i)] = computeAMG(A{i},TV0{i},false,0);
%       else
%          [MM{i},time(i)] = computeRACP(A{i},TV0{i},false,0);
%       end
%    else
%       MM{i} = P0_amg;
%    end
%    t0 = tic;
%    [~, ~, ~, it_recomp(i,:)] = gmres_RIGHT(@(v) A{i}*v, uni, gmres_restart, gmres_tol, gmres_maxit, MM{i}, []);
%    t_solve_recomp(i) = toc(t0);
%    it(i,1) = (it_recomp(i,1)-1)*gmres_restart + it_recomp(i,2);
% 
%    Pfun_new = @(v) sam_apply_left(speye(size(A{i},1)), P0_amg, v);
%    t0 = tic;
%    [~, ~, ~, it_reuse(i,:)] = gmres_RIGHT(@(v) A{i}*v, uni, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
%    t_solve_reuse(i) = toc(t0);
%    it(i,2) = (it_reuse(i,1)-1)*gmres_restart + it_reuse(i,2);
% 
%    for j = 1:i-1
%       Pfun_new = @(v) sam_apply_left(sam1{i,j}, MM{j}, v);
%       t0 = tic;
%       [~, ~, ~, it_sam1(i,:)] = gmres_RIGHT(@(v) A{i}*v, uni, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
%       t_solve_sam1(i,j) = toc(t0);
%       itsam1(i,j) = (it_sam1(i,1)-1)*gmres_restart + it_sam1(i,2);
% 
% 
%       % Pfun_new = @(v) sam_apply_left(sam2{i,j}, MM{j}, v);
%       % t0 = tic;
%       % [~, ~, ~, it_sam2(i,:)] = gmres_RIGHT(@(v) A{i}*v, uni, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
%       % t_solve_sam2(i,j) = toc(t0);
%       % itsam2(i,j) = (it_sam2(i,1)-1)*gmres_restart + it_sam2(i,2);
%    end
% end


%% Solve the system

for i = 1:sizeSeq

   uni = ones(size(A{i},1),1);

   if i > 1
      if ~lagrange
         [MM,time(i)] = computeAMG(A{i},TV0{i},false,0);
      else
         [MM,time(i)] = computeRACP(A{i},TV0{i},false,0);
      end
   else
      MM = P0_amg;
   end
   t0 = tic;
   [~, ~, ~, it_recomp(i,:)] = gmres_RIGHT(@(v) A{i}*v, uni, gmres_restart, gmres_tol, gmres_maxit, MM, []);
   t_solve_recomp(i) = toc(t0);
   it(i,1) = (it_recomp(i,1)-1)*gmres_restart + it_recomp(i,2);

   if contains(Name,"richards") 
      q = [5 8 9 10 11 12 13];
      if ismember(i,q)
         [P0_amg1,timm] = computeAMG(A{i},TV0{i},false,0);
      elseif i < q(1)
         P0_amg1 = P0_amg;
         timm = 0;
      end
   else
      P0_amg1 = P0_amg; 
      timm = 0;
   end

   Pfun_new = @(v) sam_apply_left(speye(size(A{i},1)), P0_amg1, v);
   t0 = tic;
   [~, ~, ~, it_reuse(i,:)] = gmres_RIGHT(@(v) A{i}*v, uni, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
   t_solve_reuse(i) = toc(t0) + timm;
   it(i,2) = (it_reuse(i,1)-1)*gmres_restart + it_reuse(i,2);

   if i > 1
      if contains(Name,"richards") 
         if i < 5
            sm = speye(size(A{i},1));
         elseif i >= 5 && i < 8
            sm = sam1{5};
         elseif i >= 8 
            sm = sam1{i};
         end
      else
         sm = sam1{i};
      end

      Pfun_new = @(v) sam_apply_left(sm, P0_amg, v);
      t0 = tic;
      [~, ~, ~, it_sam1(i,:)] = gmres_RIGHT(@(v) A{i}*v, uni, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
      t_solve_sam1(i) = toc(t0);
      it(i,3) = (it_sam1(i,1)-1)*gmres_restart + it_sam1(i,2);


      if strcmp("stickSlipFast",Name)
         Pfun_new = @(v) sam_apply_left(sam2{i}, P0_amg, v);
         t0 = tic;
         [~, ~, ~, it_sam2(i,:)] = gmres_RIGHT(@(v) A{i}*v, uni, gmres_restart, gmres_tol, gmres_maxit, Pfun_new, []);
         t_solve_sam2(i) = toc(t0);
         it(i,4) = (it_sam2(i,1)-1)*gmres_restart + it_sam2(i,2);
      end
   end
end


%%

hold on
for i = 1:4
   if contains("stickSlipFast",Name,'IgnoreCase',true)
      plot(1:sizeSeq, it(:,1), 'b-o', 1:sizeSeq, it(:,2), 'g-s', 1:sizeSeq, it(:,3), 'm-d', 1:sizeSeq, it(:,4), 'c-^');
   else
      plot(1:sizeSeq, it(:,1), 'b-o', 1:sizeSeq, it(:,2), 'g-s', 1:sizeSeq, it(:,3), 'm-d');
   end
   xlabel('Index');
   ylabel('Iterations');
   grid on;
   % ylim([0 150])
end

if strcmp(Name,"stickSlipFast")
   legend('RecompPrec','ReusePrec1','SAM1+P1','SAM2+P1','Location','best')
else
   legend('RecompPrec','ReusePrec1','SAM1+P1','Location','best')
end


%% ========================================================================
%% ---- Richards Post-Processing: Summary Table & Visualizations --------
%% ========================================================================
idx = [5 8 9 10 11 12 13];%1:13;%[9 11 12 13];
% 1. Map existing computation variables to explicit reporting terms
Nseq            = sizeSeq;
t_amg0          = time(1);          % Initial preconditioner construction
t_pre_L1        = zeros(1,sizeSeq);
t_pre_L1(idx) = time_prepr(idx,1); % Preprocessing overhead for SAM L1

% Execution setups
t_prec_new      = time;             % Recomputed AMG times per step
t_prec_reuse    = zeros(1, Nseq);   % Pure reuse has zero setup overhead after k=1
t_prec_L1       = zeros(1, Nseq);
t_prec_L1(idx) = time_compute1(idx);    % SAM step calculation time

% GMRES solution times
t_solve_new     = t_solve_recomp;
t_solve_reuse   = t_solve_reuse;
t_solve_L1      = t_solve_sam1;

% Total execution times (Sum of solving + local setups)
t_total_new     = sum(t_prec_new) + sum(t_solve_new);
t_total_reuse   = t_amg0 + sum(t_solve_reuse);
t_total_L1      = t_amg0 + t_pre_L1 + sum(t_prec_L1) + sum(t_solve_L1);

% Iteration allocations mapping from your 'it' matrix
iter_new        = it(:,1).';
iter_reuse      = it(:,2).';
iter_L1         = it(:,3).';

% Fake/Placeholder values for P_0 step 1 handling to align visualization steps
t_solve_L1(1)   = t_solve_reuse(1); 
iter_L1(1)      = iter_reuse(1);


%% ---- High-Quality Plots -----------------------------------------------
ks    = 1:Nseq;
c_new = [0.1 0.1 0.1];   % Dark Charcoal
c_reu = [0.8 0.1 0.1];   % Crimson Red
c_L1  = [0.1 0.3 0.8];   % Sapphire Blue

figure('Name','Richards Sequence Performance Profiles','Color','w','Position',[100 60 900 800]);

% ---- Subplot 1: GMRES Iteration Performance
subplot(3,1,1);
plot(ks, iter_new,   '-o','Color',c_new,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Recompute AMG'); hold on;
plot(ks, iter_reuse, '-s','Color',c_reu,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Reuse P_0');
plot(ks, iter_L1,    '-^','Color',c_L1, 'LineWidth',1.5,'MarkerSize',6,'DisplayName','Left SAM L1');
ylabel('GMRES iterations');
title('GMRES Iterations Count per Nonlinear Iteration Step');
legend('Location','best'); grid on;
xlim([0.5, Nseq+0.5]);

% ---- Subplot 2: Grouped Operational Costs
subplot(3,1,2);
t_prec_mat  = [t_prec_new(:) t_prec_reuse(:) t_prec_L1(:)];
t_solve_mat = [t_solve_new(:) t_solve_reuse(:) t_solve_L1(:)];

b = bar(ks, t_prec_mat + t_solve_mat, 'grouped');
b(1).FaceColor = c_new; b(2).FaceColor = c_reu; b(3).FaceColor = c_L1;
hold on;

% Translucent overlay to display setup vs solution shares distinctly
b2 = bar(ks, t_solve_mat, 'grouped');
b2(1).FaceColor = min(c_new+0.45,1); b2(1).FaceAlpha = 0.55;
b2(2).FaceColor = min(c_reu+0.45,1); b2(2).FaceAlpha = 0.55;
b2(3).FaceColor = min(c_L1+0.45,1);  b2(3).FaceAlpha = 0.55;

ylabel('Time (s)');
title('Stepwise Computation Time Breakdown (Solid Base = Setup Time, Overlay = Solve Time)');
legend([b(1) b(2) b(3)], {'Recompute','Reuse','SAM-L1'}, 'Location','best');
grid on; xlim([0.5, Nseq+0.5]);

% ---- Subplot 3: Cumulative Operational Lifespan
subplot(3,1,3);
cum_new   = cumsum(t_prec_new   + t_solve_new);
cum_reuse = cumsum(t_prec_reuse + t_solve_reuse) + t_amg0;
cum_L1    = cumsum(t_prec_L1    + t_solve_L1 + t_pre_L1)    + t_amg0;

plot(ks, cum_new,   '-o','Color',c_new,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Recompute AMG'); hold on;
plot(ks, cum_reuse, '-s','Color',c_reu,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Reuse P_0');
plot(ks, cum_L1,    '-^','Color',c_L1, 'LineWidth',1.5,'MarkerSize',6,'DisplayName','Left SAM L1');
xlabel('System Matrix Index (Step k)');
ylabel('Cumulative Time (s)');
title('Aggregated Operational Time Trajectory');
legend('Location','best'); grid on;
xlim([0.5, Nseq+0.5]);

%% ---- Separate Map Quality Plot -----------------------------------------
figure('Name','Richards Structural Mapping Deviation','Color','w','Position',[150 120 700 350]);
semilogy(2:Nseq, normmSam1(2:Nseq), '-^','Color',c_L1,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Left SAM L1');
xlabel('System Matrix Index (Step k)');
ylabel('2 \cdot ||A_1 - M_{sam} A_k||_F / (||A_1||_F + ||M_{sam} A_k||_F)');
title('Relative SAM Operator Residual Trace Quality');
legend('Location','best'); grid on;
xlim([1.5, Nseq+0.5]);