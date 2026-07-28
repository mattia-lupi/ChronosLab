clc;
clear;

nstep = 5;
step_size = 1;
epss = 1e-3;

Name = "dpMedium";%richards richardsBig stickSlipFast StickSlipOpenBig nicolas
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
      normm(i,j) = samNorm(A{j},A{i});
   end
end