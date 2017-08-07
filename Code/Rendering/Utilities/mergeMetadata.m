function [ classMap, instanceMap ] = mergeMetadata( fileName, labelMap, varargin )

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

p.parse(fileName,labelMap,varargin{:});
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
nSimplified = length(labelMap);
   
classMap = zeros(size(map));
instanceMap = zeros(size(map));
instanceId=1;
for i=1:nSimplified

    % Find all the meshes that contain the class in their name
    loc = cellfun(@(x) isempty(x)==false, (strfind(lower(mappings{2}),lower(labelMap(i).name))));
    classMap(ismember(map,mappings{1}(loc))) = labelMap(i).id;
    
    meshNames = mappings{2}(loc);
    meshIds = mappings{1}(loc);
    
    % Find all instances belonging to the category
    posA = strfind(lower(meshNames),lower('_inst_'));
    
    if any(cellfun(@isempty,posA))
        % We don't differentiate between instances
        instanceMap(ismember(map,meshIds)) = instanceId;
        instanceId = instanceId + 1;
    else
        % but if we find the keyword _inst_ then we have different
        % instances
        posB = cellfun(@(x,y) strfind(y(x(1)+length('_inst_'):end),'_'),posA,meshNames,'UniformOutput',false);
        posB = cellfun(@(x) x(1),posB,'UniformOutput',false);
        
        instances = cellfun(@(x,y,z) str2double(z(x+length('_inst_'):(x+length('_inst_')+y-2))),posA,posB,meshNames);
        instances = unique(instances);
        nInstances = length(instances);
        
        for j=1:nInstances
            loc = cellfun(@(x) isempty(strfind(lower(x),sprintf('_inst_%i',instances(j))))==false,meshNames);
            
            instanceMap(ismember(map,meshIds(loc))) = instanceId;
            instanceId = instanceId + 1;
        end
    end
    
end

classMap = uint8(classMap);
instanceMap = uint8(instanceMap);

end

