function [ camPos, camLookAt, filmDist ] = nnCameraPos( objects, camHeight, camDist, camOr, camDef, lensFile, varargin )

p = inputParser;
p.addOptional('refId',1);
p.parse(varargin{:});
refId = p.Results.refId;


nHeight = length(camHeight);
nDist = length(camDist);
nOr = length(camOr);
nVp = length(refId);
nDef = length(camDef);

nPos = nHeight*nDist*nOr*nVp*nDef;

camPos = zeros(nPos,3);
filmDist = zeros(nPos,1);
camLookAt = zeros(nPos,1);

for ch=1:nHeight
    for cd=1:nDist
        for co=1:nOr
            
            cx = camDist(cd)*sind(camOr(co));
            cy = camDist(cd)*cosd(camOr(co));
            
            for la=1:nVp
                for de=1:nDef
                
                    id = sub2ind([nHeight, nDist, nOr, nVp, nDef],ch, cd, co, la, de);
                    
                    if isstruct(objects)
                        objPosition = eval(objects(refId(la)).position);
                    else
                        objPosition = objects(refId(la),:);
                    end
                    
                    
                    
                    camPos(id,1) = cx + objPosition(1);
                    camPos(id,2) = cy + objPosition(2);
                    camPos(id,3) = camHeight(ch);
                    
                    filmDist(id) = focusLens(lensFile,(camDist(cd) + camDef(de))*1000);
                    
                    camLookAt(id,1:2) = objPosition(1:2);
                    camLookAt(id,3) = camHeight(ch);
                end
            end
        end
    end
end


end

