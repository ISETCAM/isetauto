function rootPath=nngenRootPath()
% Return the path to the root NN_Camera_Generalization directory
%
% This function must reside in the directory at the base of the
% NN_Camera_Generalization directory structure.  It is used to determine the
% location of various sub-directories.
% 
% Example:
%   fullfile(isetRootPath,'data')

rootPath=which('nngenRootPath');

rootPath=fileparts(rootPath);

end
