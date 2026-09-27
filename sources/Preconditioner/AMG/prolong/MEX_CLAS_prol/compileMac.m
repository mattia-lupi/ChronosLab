% Get Homebrew Paths
[status, cmdout] = system('/opt/homebrew/bin/brew --prefix libomp');
if status ~= 0
    error('Error: Could not find libomp via Homebrew. Make sure it is installed.');
end
basePath = strtrim(cmdout);
omp_inc = ['-I' fullfile(basePath, 'include')];

% Define the STATIC Library Path
% We link against libomp.a directly to hide it from MATLAB's runtime.
omp_static = fullfile(basePath, 'lib', 'libomp.a');
if ~isfile(omp_static)
    error('Static library libomp.a not found at %s.', omp_static);
end

% Determine host macOS deployment target
[status, ver_str] = system('sw_vers -productVersion');
if status == 0
    mac_ver = strtrim(ver_str);
else
    mac_ver = '15.0';
end
ver_flag = ['-mmacosx-version-min=' mac_ver];

% Compile command for cpt_Prolongation_Classical
mex('-silent', ...
    '-O', ...
    '-R2018a',...
    omp_inc, ...
    ['CXXFLAGS="$CXXFLAGS -std=c++14 -O2 -Xpreprocessor -fopenmp ' ver_flag ' -fPIC -I./include/"'], ...
    ['LDFLAGS="$LDFLAGS ' ver_flag ' -O2"'], ...
    'cpt_Prolongation_Classical.cpp', ...
    'Classical_prolongation.cpp', ...
    'ProlStripe_Classical.cpp', ...
    'ir_heapsort.cpp', ...
    omp_static);
