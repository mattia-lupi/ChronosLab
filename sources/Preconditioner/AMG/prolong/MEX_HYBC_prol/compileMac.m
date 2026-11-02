function compileMac(omp_inc, omp_lib)
if nargin < 2 || isempty(omp_inc) || isempty(omp_lib)
    if exist('findOpenMPMacOs', 'file') ~= 2
        this_dir = fileparts(mfilename('fullpath'));
        curr = this_dir;
        while ~isempty(curr) && ~isfile(fullfile(curr, 'compileAll.m'))
            parent = fileparts(curr);
            if strcmp(parent, curr), break; end
            curr = parent;
        end
        if isfile(fullfile(curr, 'compileAll.m'))
            addpath(fullfile(curr, 'sources', 'utils'));
        end
    end
    [omp_inc, omp_lib] = findOpenMPMacOs();
end

% Determine host macOS deployment target
[status, ver_str] = system('sw_vers -productVersion');
if status == 0
    mac_ver = strtrim(ver_str);
else
    mac_ver = '15.0';
end
ver_flag = ['-mmacosx-version-min=' mac_ver];

% Compile command for cpt_Prolongation_HybridClas
mex('-silent', ...
    '-O', ...
    '-R2018a',...
    omp_inc, ...
    ['CXXFLAGS="$CXXFLAGS -std=c++14 -O2 -Xpreprocessor -fopenmp ' ver_flag ' -fPIC -I./include/"'], ...
    ['LDFLAGS="$LDFLAGS ' ver_flag ' -O2"'], ...
    'cpt_Prolongation_HybridClas.cpp', ...
    'HybridClas_prolongation.cpp', ...
    'ProlStripe_HybridClas.cpp', ...
    'ir_heapsort.cpp', ...
    'find_bad_FS.cpp', ...
    'find_dist2CC.cpp', ...
    omp_lib);
end
