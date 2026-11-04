function win_cfg = findCompilerWindows(silent)
% FINDCOMPILERWINDOWS Resolves C/C++ compiler and OpenMP settings for Windows.
%
%   win_cfg = findCompilerWindows()
%   win_cfg = findCompilerWindows(silent)
%
% Returns a struct win_cfg with fields:
%   win_cfg.compiler: 'msvc' or 'mingw'
%   win_cfg.omp_cflags: OpenMP compiler flag
%   win_cfg.omp_lflags: OpenMP linker flag
%   win_cfg.cxx14_flags: C++14/standards flags
%   win_cfg.extra_libs: cell array of extra libraries to link (e.g. {omp_lib} for MSVC)

    if nargin < 1
        silent = false;
    end

    arch = computer('arch');
    if ~strcmp(arch, 'win64')
        error('findCompilerWindows is intended for win64 architectures. Detected: %s', arch);
    end

    % Query configured C++ mex compiler in MATLAB
    cc = mex.getCompilerConfigurations('C++', 'Selected');
    if isempty(cc)
        % Fallback to C compiler configuration if C++ is not explicitly selected
        cc = mex.getCompilerConfigurations('C', 'Selected');
    end

    compiler_name = '';
    compiler_id = '';
    if ~isempty(cc)
        compiler_name = cc(1).Name;
        compiler_id = cc(1).Manufacturer;
    end

    % Detect if MinGW / GCC or MSVC
    is_mingw = ~isempty(regexpi(compiler_name, 'mingw|gcc|gnu')) || ...
               ~isempty(regexpi(compiler_id, 'gnu|mingw'));

    if is_mingw
        compiler_type = 'mingw';
        % For MinGW-w64 (GCC), OpenMP uses -fopenmp.
        % Link statically against libgcc, libstdc++, and libgomp to avoid runtime
        % DLL version conflicts with MATLAB's bundled libraries.
        % Note: Do NOT use global '-static' because that prevents linking against
        % MATLAB's dynamic import libraries (mwblas, mwlapack, mex, mx).
        mingw_libdir = fullfile(matlabroot, 'extern', 'lib', arch, 'mingw64');
        omp_cflags = '-fopenmp';
        if isfolder(mingw_libdir)
            omp_lflags = sprintf('-L"%s" -fopenmp -static-libgcc -static-libstdc++ -Wl,-Bstatic -lgomp -Wl,-Bdynamic', mingw_libdir);
        else
            omp_lflags = '-fopenmp -static-libgcc -static-libstdc++ -Wl,-Bstatic -lgomp -Wl,-Bdynamic';
        end
        cxx14_flags = '-std=c++14';
        extra_libs = {};
    else
        % Default to MSVC
        compiler_type = 'msvc';
        omp_lib = fullfile(matlabroot, 'bin', arch, 'libiomp5md.lib');
        if ~isfile(omp_lib)
            error('MATLAB internal OpenMP library not found at: %s', omp_lib);
        end
        omp_cflags = '/openmp:llvm';
        omp_lflags = '/NODEFAULTLIB:libomp.lib /NODEFAULTLIB:vcomp.lib';
        cxx14_flags = '/std:c++14';
        extra_libs = {omp_lib};
    end

    win_cfg.compiler = compiler_type;
    win_cfg.omp_cflags = omp_cflags;
    win_cfg.omp_lflags = omp_lflags;
    win_cfg.cxx14_flags = cxx14_flags;
    win_cfg.extra_libs = extra_libs;

    if ~silent
        fprintf('[Windows Compiler] Detected %s (%s)\n', win_cfg.compiler, compiler_name);
        fprintf('                   omp_cflags: %s\n', win_cfg.omp_cflags);
        fprintf('                   omp_lflags: %s\n', win_cfg.omp_lflags);
        if ~isempty(win_cfg.extra_libs)
            fprintf('                   Extra Lib: %s\n', win_cfg.extra_libs{1});
        end
    end
end
