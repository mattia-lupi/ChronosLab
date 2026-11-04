function win_cfg = findCompilerWindows(silent)
% FINDCOMPILERWINDOWS Resolves C/C++ compiler and OpenMP/BLAS settings for Windows.
%
%   win_cfg = findCompilerWindows()
%   win_cfg = findCompilerWindows(silent)
%
% Returns a struct win_cfg with fields:
%   win_cfg.compiler: 'msvc' or 'mingw'
%   win_cfg.omp_cflags: OpenMP compiler flag
%   win_cfg.omp_lflags: OpenMP linker flag
%   win_cfg.cxx14_flags: C++14/standards flags
%   win_cfg.blas_lib: full path or flag for mwblas
%   win_cfg.lapack_lib: full path or flag for mwlapack
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
        % Link statically against libgcc, libstdc++, libgomp, and libwinpthread to avoid runtime
        % DLL dependencies (libgomp-1.dll, libstdc++-6.dll) when MinGW's bin directory is not on Windows PATH.
        % Note: Static libraries must be passed via LINKLIBS so they appear AFTER the
        % object files on the linker command line, ensuring undefined symbols from the object files
        % are resolved from the static archives rather than dynamic DLLs.
        omp_cflags = '-fopenmp';
        omp_lflags = '-static-libgcc -static-libstdc++';
        cxx14_flags = '-std=c++14 -include stdint.h';
        extra_libs = {'LINKLIBS=$LINKLIBS -Wl,-Bstatic -lstdc++ -lgomp -lpthread -Wl,-Bdynamic'};

        % MinGW requires the explicit path to libmwblas and libmwlapack libraries
        % to avoid "MEX cannot find library 'mwblas'" when using -l shorthand.
        mingw_libdir = fullfile(matlabroot, 'extern', 'lib', arch, 'mingw64');
        msvc_libdir  = fullfile(matlabroot, 'extern', 'lib', arch, 'microsoft');

        if isfile(fullfile(mingw_libdir, 'libmwblas.lib'))
            blas_lib = fullfile(mingw_libdir, 'libmwblas.lib');
            lapack_lib = fullfile(mingw_libdir, 'libmwlapack.lib');
        elseif isfile(fullfile(mingw_libdir, 'libmwblas.a'))
            blas_lib = fullfile(mingw_libdir, 'libmwblas.a');
            lapack_lib = fullfile(mingw_libdir, 'libmwlapack.a');
        elseif isfile(fullfile(msvc_libdir, 'libmwblas.lib'))
            blas_lib = fullfile(msvc_libdir, 'libmwblas.lib');
            lapack_lib = fullfile(msvc_libdir, 'libmwlapack.lib');
        else
            blas_lib = '-lmwblas';
            lapack_lib = '-lmwlapack';
        end
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
        blas_lib = '-lmwblas';
        lapack_lib = '-lmwlapack';
    end

    win_cfg.compiler = compiler_type;
    win_cfg.omp_cflags = omp_cflags;
    win_cfg.omp_lflags = omp_lflags;
    win_cfg.cxx14_flags = cxx14_flags;
    win_cfg.blas_lib = blas_lib;
    win_cfg.lapack_lib = lapack_lib;
    win_cfg.extra_libs = extra_libs;

    if ~silent
        fprintf('[Windows Compiler] Detected %s (%s)\n', win_cfg.compiler, compiler_name);
        fprintf('                   omp_cflags: %s\n', win_cfg.omp_cflags);
        fprintf('                   omp_lflags: %s\n', win_cfg.omp_lflags);
        fprintf('                   BLAS:       %s\n', win_cfg.blas_lib);
        fprintf('                   LAPACK:     %s\n', win_cfg.lapack_lib);
        if ~isempty(win_cfg.extra_libs)
            fprintf('                   Extra Lib:  %s\n', win_cfg.extra_libs{1});
        end
    end
end
