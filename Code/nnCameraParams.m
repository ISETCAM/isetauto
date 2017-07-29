function [cameraPosition, filmDistanceVec, cameraDistanceVec] = ...
    nnCameraParams(carPosition, cameraDefocus,cameraHeight,cameraDistance, cameraOrientation, lensFile)
% Set critical parameters for the neural network generalization
%
% We make all of these three parameters together because they depend on one
% another.  What comes back is a set of 
%  cameraPositions (N x 3), filmDistances (N x 1) and cameraDistances (N x 1)
% The N is for each condition.
% The number of conditions depends on the other parameters as well
%   defocus, height, carPosition, cameraOrientation

%%
sz = [length(cameraDefocus) length(cameraHeight) length(cameraDistance) length(cameraOrientation)];

cameraPosition    = zeros(prod(sz),3);  % Three camera position variables
filmDistanceVec   = zeros(prod(sz),1);
cameraDistanceVec = zeros(prod(sz),1);

for cdef=1:length(cameraDefocus)
    for ch=1:length(cameraHeight)
        for cdd=1:length(cameraDistance)
            for co=1:length(cameraOrientation)
                
                % This would be the counter variable
                loc = sub2ind(sz,cdef,ch,cdd,co);
                
                cx = cameraDistance(cdd)*sind(cameraOrientation(co));
                cy = cameraDistance(cdd)*cosd(cameraOrientation(co));
                
                cameraPosition(loc,1) = cx + carPosition(1);
                cameraPosition(loc,2) = cy + carPosition(2);
                cameraPosition(loc,3) = cameraHeight(ch);
                
                filmDistanceVec(loc) = focusLens(lensFile,(max(cameraDistance(cdd)+cameraDefocus(cdef),0.1))*1000);
                cameraDistanceVec(loc) = cameraDistance(cdd);
                
            end
        end
    end
end

end


% I think the code below can be an ndgrid
% I think what is being built here is
% cameraPosition, filmDistance and cameraDistance
%
%                 [X1,X2,X3,X4] = ndgrid(1:length(cameraDefocus),1:length(cameraHeight), 1:length(cameraDistance), 1:length(cameraOrientation));
%                 pList = [X1(:) X2(:) X3(:) X4(:)];
%                 % dim 1 is cameraDefocus
%                 % dim 2 is cameraHeight
%                 % dim 3 is cameraDistance
%                 % dim 4 is cameraOrientation
%
%                 for ii=1:size(pList,1)
%                     p = pList(ii,:);
%                     cx = cameraDistance(p(3))*sind(cameraOrientation(p(4)));
%                     cy = cameraDistance(p(3))*cosd(cameraOrientation(p(4)));
%                     cameraPosition(ii,1) =  cx + carPosition(ap,1);
%                    cameraPosition(ii,2)  =  cy + carPosition(ap,2);
%
%                    filmDistanceVec(ii)   = focusLens(lensFile,(max(cameraDistance(p(3)) + cameraDefocus(p(1)),0.1))*1000);
%
%                    cameraDistanceVec(ii) = cameraDistance(p(3));
%                 end
%
