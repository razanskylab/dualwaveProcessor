%% clear workspace, command window and open connections
clear, clc; 
close all;
rehash;
if ~isempty(instrfind)
  fclose(instrfind);
  delete(instrfind);
end


