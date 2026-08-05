function [A, TV0, b, lagrange, sizeSeq] = loadMatrices(Name)
    % LOADMATRICES Loads dataset files and formats vectors/matrices
    files = dir(fullfile('mats', Name, [char(Name) '_*.mat']));
    fileNames = {files.name};
    if isempty(fileNames)
        error('No files found for dataset "%s" in mats/%s/', Name, Name);
    end

    lagrange = ~contains(Name, ["richards", "dp", "symMech"]);

    % Extract trailing numbers for numerical sorting
    tokens = regexp(fileNames, Name + '_(\d+)\.mat', 'tokens', 'once');
    nums = cellfun(@(x) str2double(x{1}), tokens);
    [~, sortIdx] = sort(nums);
    fileNames = fileNames(sortIdx);

    startVal = 1;
    fileNames = fileNames(startVal:end);

    sizeSeq = min(24, length(fileNames));
    if strcmp(Name, "nicolas")
       sizeSeq = min(30, sizeSeq);
    end

    A   = cell(sizeSeq, 1);
    TV0 = cell(sizeSeq, 1);
    b   = cell(sizeSeq, 1);

    for i = 1:sizeSeq
        fname = fullfile('mats', Name, fileNames{i});
        data = load(fname);

        if ~strcmp(Name, "nicolas")
           A{i} = data.Amat;
           b{i} = data.b;
        else
           A{i} = data.A;
           b{i} = ones(size(A{1}, 1), 1);
        end

        if ~contains(Name, "richards")
           if strcmp(Name, "StickSlipOpenBig")
              TV0{i} = data.TV;
           else
              TV0{i} = data.TV0;
           end
        else
           TV0{i} = ones(size(A{i}, 1), 1);
        end
    end
end