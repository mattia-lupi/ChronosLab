clc;
clear;
close all;

list = {'sources/Preconditioner/AMG/smoother/MEX_SYM_AFSAI/', ...
        'sources/Preconditioner/AMG/smoother/MEX_NSY_RFSAI/', ...
        'sources/Preconditioner/AMG/prolong/MEX_HYBC_prol/', ...
        'sources/Preconditioner/AMG/prolong/MEX_CLAS_prol/', ...
        'sources/Preconditioner/AMG/prolong/MEX_BAMG_Prol/', ...
        'sources/Preconditioner/AMG/prolong/MEX_EXTI_Prol/', ...
        'sources/Preconditioner/AMG/prolong/EMIN/MEX_EMIN/', ...
        'sources/Preconditioner/AMG/filter/MEX_CLev_Filter/', ...
        'sources/Preconditioner/SAM/MEX_Adaptive',...
        'sources/Preconditioner/AMG/filter/MEX_Prol_Filter/'};

home_dir = pwd;
sys_arch = computer('arch');

try
   for idx = 1:length(list)
       folder = list{idx};
       fprintf('Compiling MEX files in %s\n', folder);
       cd(folder);
       
       switch sys_arch
           case 'win64'
               compileWin
           case 'maca64'
               compileMac
           case 'glnxa64'
               compile
           otherwise
               error('Architecture [%s] is explicitly unsupported. Allowed configurations: win64, maca64, glnxa64.', sys_arch);
       end
       
       cd(home_dir);
   end
catch ME
   cd(home_dir);
   fprintf(2, 'Compilation failed: %s\n', ME.message);
end
