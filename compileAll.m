clc;
clear;
close all;

list = { ...
    fullfile('sources', 'Preconditioner', 'AMG', 'smoother', 'MEX_SYM_AFSAI'), ...
    fullfile('sources', 'Preconditioner', 'AMG', 'smoother', 'MEX_NSY_RFSAI'), ...
    fullfile('sources', 'Preconditioner', 'AMG', 'prolong', 'MEX_HYBC_prol'), ...
    fullfile('sources', 'Preconditioner', 'AMG', 'prolong', 'MEX_CLAS_prol'), ...
    fullfile('sources', 'Preconditioner', 'AMG', 'prolong', 'MEX_BAMG_Prol'), ...
    fullfile('sources', 'Preconditioner', 'AMG', 'prolong', 'MEX_EXTI_Prol'), ...
    fullfile('sources', 'Preconditioner', 'AMG', 'prolong', 'EMIN', 'MEX_EMIN'), ...
    fullfile('sources', 'Preconditioner', 'AMG', 'filter', 'MEX_CLev_Filter'), ...
    fullfile('sources', 'Preconditioner', 'SAM', 'MEX_Adaptive'), ...
    fullfile('sources', 'Preconditioner', 'AMG', 'filter', 'MEX_Prol_Filter') ...
};

home_dir = pwd;
[~, base_folder] = fileparts(home_dir);

% Ensures MATLAB returns to home_dir even on error or Ctrl+C
cleanUpObj = onCleanup(@() cd(home_dir)); 

sys_arch = computer('arch');

if strcmp(sys_arch, 'maca64')
    addpath(fullfile(home_dir, 'sources', 'utils'));
    [omp_inc, omp_lib] = findOpenMPMacOs();
end

for idx = 1:length(list)
    folder = fullfile(home_dir, list{idx});
    display_folder = fullfile(base_folder, list{idx});
    fprintf('Compiling MEX files in %s\n', display_folder);
    
    cd(folder);
    
    switch sys_arch
        case 'win64'
            compileWin;
        case 'maca64'
            compileMac(omp_inc, omp_lib);
        case 'glnxa64'
            compile;
        otherwise
            error('Architecture [%s] is unsupported. Allowed: win64, maca64, glnxa64.', sys_arch);
    end
end

cd(home_dir)

fprintf('\nAll MEX files compiled successfully.\n');