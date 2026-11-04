function compileWin(win_cfg)
if nargin < 1 || isempty(win_cfg)
    if exist('findCompilerWindows', 'file') ~= 2
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
    win_cfg = findCompilerWindows();
end

if strcmp(win_cfg.compiler, 'mingw')
    cflags = sprintf('CXXFLAGS="$CXXFLAGS %s %s -ffast-math"', win_cfg.cxx14_flags, win_cfg.omp_cflags);
    lflags = sprintf('LDFLAGS="$LDFLAGS %s"', win_cfg.omp_lflags);
else
    cflags = sprintf('COMPFLAGS="$COMPFLAGS %s %s /fp:fast"', win_cfg.cxx14_flags, win_cfg.omp_cflags);
    lflags = sprintf('LINKFLAGS="$LINKFLAGS %s"', win_cfg.omp_lflags);
end

% Compile
mex('-O', ...
    '-R2018a', ...
    '-lmwblas', ...
    '-lmwlapack', ...
    cflags, ...
    lflags, ...
    'sam_adaptive_left_mex.cpp', ...
    'cpt_sam_adaptive_left.cpp', ...
    'qr_functions.cpp', ...
    'find_stuff.cpp', ...
    'cpt_resRho.cpp', ...
    win_cfg.extra_libs{:});
end
