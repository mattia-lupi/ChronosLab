function [omp_inc, omp_lib] = findOpenMPMacOs(silent)
% FINDOPENMPMACOS Resolves OpenMP include flags and library path for macOS.
%
%   [omp_inc, omp_lib] = findOpenMPMacOs()
%   [omp_inc, omp_lib] = findOpenMPMacOs(silent)
%
% Priority:
%   1. MATLAB internal OpenMP shared library (libomp.dylib) + MATLAB/Homebrew omp.h
%   2. Homebrew OpenMP static library (libomp.a) fallback

    if nargin < 1
        silent = false;
    end

    omp_inc = '';
    omp_lib = '';
    used_type = '';
    arch = computer('arch');

    % -----------------------------------------------------------------
    % 1. Try MATLAB internal OpenMP (Shared Library: libomp.dylib)
    % -----------------------------------------------------------------
    matlab_lib_candidates = {
        fullfile(matlabroot, 'bin', arch, 'libomp.dylib');
        fullfile(matlabroot, 'toolbox', 'eml', 'externalDependency', 'omp', arch, 'lib', 'libomp.dylib')
    };

    matlab_lib = '';
    for i = 1:numel(matlab_lib_candidates)
        if isfile(matlab_lib_candidates{i})
            matlab_lib = matlab_lib_candidates{i};
            break;
        end
    end

    matlab_hdr_candidates = {
        fullfile(matlabroot, 'toolbox', 'eml', 'externalDependency', 'omp', arch, 'include');
        fullfile(matlabroot, 'toolbox', 'coder', 'clang_api', 'llvm-include', arch);
        fullfile(matlabroot, 'extern', 'include')
    };

    matlab_inc = '';
    for i = 1:numel(matlab_hdr_candidates)
        if isfile(fullfile(matlab_hdr_candidates{i}, 'omp.h'))
            matlab_inc = matlab_hdr_candidates{i};
            break;
        end
    end

    % If header not found in primary candidate locations, do a targeted toolbox search
    if isempty(matlab_inc) && ~isempty(matlab_lib)
        hits = dir(fullfile(matlabroot, 'toolbox', '**', 'omp.h'));
        if ~isempty(hits)
            matlab_inc = hits(1).folder;
        end
    end

    if ~isempty(matlab_lib)
        if ~isempty(matlab_inc)
            omp_lib = matlab_lib;
            omp_inc = ['-I' matlab_inc];
            used_type = 'MATLAB internal OpenMP (shared)';
        else
            % Fallback for omp.h only: check Homebrew
            brew_prefix = '';
            if isfolder('/opt/homebrew/opt/libomp')
                brew_prefix = '/opt/homebrew/opt/libomp';
            else
                [status, cmdout] = system('/opt/homebrew/bin/brew --prefix libomp');
                if status == 0
                    brew_prefix = strtrim(cmdout);
                end
            end
            if ~isempty(brew_prefix) && isfile(fullfile(brew_prefix, 'include', 'omp.h'))
                omp_lib = matlab_lib;
                omp_inc = ['-I' fullfile(brew_prefix, 'include')];
                used_type = 'MATLAB internal OpenMP (shared with Homebrew header)';
            end
        end
    end

    % -----------------------------------------------------------------
    % 2. Fallback: Homebrew OpenMP (Static Library: libomp.a)
    % -----------------------------------------------------------------
    if isempty(omp_lib) || isempty(omp_inc)
        brew_prefix = '';
        if isfolder('/opt/homebrew/opt/libomp')
            brew_prefix = '/opt/homebrew/opt/libomp';
        else
            [status, cmdout] = system('/opt/homebrew/bin/brew --prefix libomp');
            if status ~= 0
                error('OpenMP not found in MATLAB and Homebrew libomp is unavailable.');
            end
            brew_prefix = strtrim(cmdout);
        end

        brew_static = fullfile(brew_prefix, 'lib', 'libomp.a');
        if ~isfile(brew_static)
            error('Static library libomp.a not found at %s.', brew_static);
        end

        brew_header = fullfile(brew_prefix, 'include');
        if ~isfile(fullfile(brew_header, 'omp.h'))
            error('OpenMP header omp.h not found at %s.', brew_header);
        end

        omp_inc = ['-I' brew_header];
        omp_lib = brew_static;
        used_type = 'Homebrew OpenMP (static)';
    end

    if ~silent
        fprintf('[OpenMP] Using %s\n', used_type);
        fprintf('         Include: %s\n', omp_inc);
        fprintf('         Library: %s\n', omp_lib);
    end
end
