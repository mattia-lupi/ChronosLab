function initChronos(outputFlag)

   % List of folders to check
   foldList = ["sources", "Tests"];
   
   if nargin == 0
     outputFlag = true;
   end
   
   % Check if any of the folders is already on the MATLAB path
   needRestore = false;
   for k = 1:numel(foldList)
     if contains(path, fullfile('ChronosLab',foldList(k)))
       needRestore = true;
       break;
     end
   end
   
   % Restore path only if needed
   if needRestore
     if outputFlag
       fprintf('Restoring default paths \n');
     end
     restoredefaultpath;
   end
   
   % set the root to the gres directory
   chronos_root = fileparts(mfilename('fullpath'));
   setappdata(0,'chronos_root', chronos_root);
   
   % Add path
   for k = 1:numel(foldList)
     addpath(genpath(fullfile(chronos_root, foldList(k))));
   end
   
   if outputFlag
      fprintf('Chronos sources added to path\n');
   end

end

