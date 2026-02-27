function GB_used = GetWhosMemory()
   vars = evalin('caller', 'whos');
   totalBytes = sum([vars.bytes]);
   GB_used = totalBytes / (1024^3); % In GB
end
