function [ classMap, instanceMap ] = mergeMetadata( fileName, classesToMerge, varargin )

% [ simplifiedMap ] = mergeMetadata( fileName, classesToMerge, ... )
%
% This function accepts a path to a .mat file that represents the output
% from a 'mesh' or 'material' metadata PBRT renderer. It outputs a
% simplified map with fewer meshes or materials than in the original file.
% The classes are groupped accoridng the presence of common strings in
% their class descriptions in the accompanying .txt file. For example if
% classesToMerge = {'car','city'}, then all class names contaning the
% substring 'car' will be assigned one class and all class names containing
% the substring 'city' the other. 
%
% It is assumed tha strings in classesToMerge are mutually exclusive, i.e.
% there is no class that could contain two or more strings from
% classesToMerge cell array.
%
% Copyright, Henryk Blasinski 2017


p = inputParser;
p.addRequired('fileName');
p.addRequired('classesToMerge',@isstruct);
p.addOptional('mode','mesh');

p.parse(fileName,classesToMerge,varargin{:});
mode = p.Results.mode;

load(fileName);
img = multispectralImage(:,:,1);
map = uint32(img/radiometricScaleFactor);

[path, fName] = fileparts(fileName);
switch mode
    case 'mesh'
        textFile = fullfile(path,sprintf('%s_mesh.txt',fName));
    case 'texture'
        textFile = fullfile(path,sprintf('%s_texture.txt',fName));
end

fid = fopen(textFile);
mappings = textscan(fid,'%u %s');
fclose(fid);
nSimplified = length(classesToMerge);
   
classMap = zeros(size(map));
instanceMap = zeros(size(map));
instanceId=1;
for i=1:nSimplified

    loc = cellfun(@(x) isempty(x)==false, (strfind(lower(mappings{2}),lower(classesToMerge(i).name))));
    classMap(ismember(map,mappings{1}(loc))) = classesToMerge(i).id;
    
    % We assume that there is only one instance per class !!!!!
    instanceMap(ismember(map,mappings{1}(loc))) = instanceId;
    instanceId = instanceId + 1;
end

classMap = uint8(classMap);
instanceMap = uint8(instanceMap);

end

