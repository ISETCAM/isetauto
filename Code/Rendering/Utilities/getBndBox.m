function [objects] = getBndBox(labelImage, instanceImage, labelMap)

objects = {};

h = size(labelImage,1);
w = size(labelImage,2);

% Get instance id's, we assume that 0 is the background class,
% which we ignore.
instanceIDs = unique(instanceImage);
instanceIDs = instanceIDs(instanceIDs > 0);

objectId = 1;
for i=1:length(instanceIDs)
   indicator = instanceImage == instanceIDs(i);
   
   xSpread = sum(indicator);
   xIndices = find(xSpread > 0);
   
   ySpread = sum(indicator,2);
   yIndices = find(ySpread > 0);
   
   classId = labelImage(indicator);
   classId = classId(1);
   
   objects{objectId}.name = labelMap([labelMap(:).id] == classId).name;
   objects{objectId}.bndbox.xmin = min(xIndices);
   objects{objectId}.bndbox.xmax = max(xIndices);
   objects{objectId}.bndbox.ymin = min(yIndices);
   objects{objectId}.bndbox.ymax = max(yIndices);
   
   
   % Occlusions   
   ccomp = bwconncomp(indicator);
   if ccomp.NumObjects > 1
       occluded = 1;
   else
       occluded = 0;
   end
   
   % Truncations
   if (objects{objectId}.bndbox.xmin == 0 || objects{objectId}.bndbox.ymin == 0 || ...
           objects{objectId}.bndbox.xmax == w || objects{objectId}.bndbox.ymax == h)
       truncated = 1;
   else
       truncated = 0;
   end
   
   objects{objectId}.difficult = 0;
   objects{objectId}.occluded = occluded;
   objects{objectId}.pose = 'Unspecified';
   objects{objectId}.truncated = truncated;
   
   objectId = objectId + 1; 
end




end