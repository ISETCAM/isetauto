%% Compare lens renderings of different cars
%
% 
% Henryk Blasinski, SCIEN Stanford, 2017

%% Scene description

% Initialize ISET related variables
ieInit;

% Sets up related to the car renderings
nnConstants;

% Set-up RenderToolbox4

hints = nnHintsInit;

%% Simulation parameters
%
% Negative z is up.
% Scene is about 200x200m, units are mm.
% However we should specify meters, as they are automatically converted to
% mm in remodellers.

cameraType = {'pinhole','lens','lens'}; 
lensType = {'tessar.22deg.6.0mm','tessar.22deg.6.0mm','2el.XXdeg.6.0mm'};
microlens = {[0,0],[0,0],[0,0]};

mode = {'radiance'};

fNumber = 2.8;
filmDiag = (1/3.6)*25.4;

diffraction = {'false','true'};
chrAber = {'false','true'};

% diffraction = {'false','true'};
% chrAber = {'false','true'};

pixelSamples = 128;
shadowDirection = [-0.5 -1 1;];

cameraDistance = [10 20];
cameraOrientation = [0];
cameraPan = [0];
cameraTilt = [0];
cameraRoll = [0];

cameraHeight = -1.5;
cameraDefocus = 0;

nCarPositions = 1;
carOrientation = [30];

maxCars = 1;
maxCities = 1;

% These are the variable names used in the conditionsFile.  See
%  https://github.com/RenderToolbox/RenderToolbox4/wiki/Conditions-File-
% Some of these are standard.  Some are selected here.
names = {'imageName','cameraType','lensType','mode','pixelSamples','filmDist','filmDiag','cameraPosition',...
    'shadowDirection','microlensDim','cameraLookAt','fNumber','carPosition','carOrientation','fog',...
    'diffraction','chromaticAberration','cameraPan','cameraTilt','cameraRoll'};

%% Check
assert(length(cameraType)  == length(lensType));
assert(length(cameraType)  == length(microlens));
assert(length(diffraction) == length(chrAber));


%% Choose renderer options 
%
% Select the working folder and copy critical files there
% Not sure why the lens file is not handled automatically.  This caused a
% problem for me with the calibration rendering, too. (BW).

resourceFolder = rtbWorkingFolder('folderName','resources',...
    'rendererSpecific',false,...
    'hints',hints);

% Copy resources
lensFiles = fullfile(lensDir,strcat(lensType,'.dat'));
for i=1:length(lensFiles)
    copyfile(lensFiles{i},resourceFolder);
end

% Copy sky map
skyFile = fullfile(assetDir,'City','*.exr');
copyfile(skyFile,resourceFolder);

% Copy D65 spectrum
[wave, d65] = rtbReadSpectrum(fullfile(rtbRoot,'RenderData','D65.spd'));
d65 = 100*d65;
rtbWriteSpectrumFile(wave,d65,fullfile(resourceFolder,'D65.spd'));


%% Choose files to render
sceneID = 1;

for cityId=1:maxCities
    sceneFile = sprintf('City_%i.obj',cityId);
    parentSceneFile = fullfile(assetDir,'City',sceneFile);
    
    [cityScene, elements] = mexximpCleanImport(parentSceneFile,...
        'ignoreRootTransform',true,...
        'flipUVs',true,...
        'imagemagicImage','hblasins/imagemagic-docker',...
        'toReplace',{'jpg','png','tga'},...
        'options','-gamma 0.45',...
        'targetFormat','exr',...
        'makeLeftHanded',true,...
        'flipWindingOrder',true,...
        'workingFolder',resourceFolder);
    
    % Defines valid ranges for the car positions
    carPosition = zeros(nCarPositions,3);
    for i=2:nCarPositions
        carPosition(i,:) = drawCarPosition(cityId);
    end
    
    % There seem to be 1 car and 1 city
    for carId=1:maxCars
        carFile = sprintf('Car_%i.obj',carId);
        parentSceneFile = fullfile(assetDir,car2directory{carId},carFile);
        
        [carScene, elements] = mexximpCleanImport(parentSceneFile,...
            'ignoreRootTransform',true,...
            'flipUVs',true,...
            'imagemagicImage','hblasins/imagemagic-docker',...
            'toReplace',{'jpg','png','tga'},...
            'targetFormat','exr',...
            'makeLeftHanded',true,...
            'flipWindingOrder',true,...
            'workingFolder',resourceFolder);
        
        scene = mexximpCombineScenes(cityScene,carScene,...
            'insertTransform',mexximpTranslate([0 0 0]),...
            'cleanupTransform',mexximpTranslate([0 0 0]));
        

        for ap=1:nCarPositions;
            
            for lt=1:length(lensType)
                
                % It seems like we get a clean copy of the Conditions.txt file
                % here?  Later, we write it out.
                conditionsFile = fullfile(resourceFolder,'Conditions.txt');

                % This is one place where names is used.  
                values = cell(1,numel(names));
                cntr = 1;
                
                lensFile = fullfile(lensDir,sprintf('%s.dat',lensType{lt}));
                
                sz = [length(cameraDefocus) length(cameraHeight) length(cameraDistance) length(cameraOrientation)];
                
                cameraPosition = zeros(prod(sz),3);
                filmDistanceVec = zeros(prod(sz),1);
                cameraDistanceVec = zeros(prod(sz),1);
                
                % I think the code below can be an ndgrid
                % I think what is being built here is
                % cameraPosition, filmDistance and cameraDistance
                %
                [X1,X2,X3,X4] = ndgrid(1:length(cameraDefocus),1:length(cameraHeight), 1:length(cameraDistance), 1:length(cameraOrientation));
                pList = [X1(:) X2(:) X3(:) X4(:)];
                % dim 1 is cameraDefocus
                % dim 2 is cameraHeight
                % dim 3 is cameraDistance
                % dim 4 is cameraOrientation
                
                for ii=1:size(pList,1)
                    p = pList(ii,:);
                    cx = cameraDistance(p(3))*sind(cameraOrientation(p(4)));
                    cy = cameraDistance(p(3))*cosd(cameraOrientation(p(4)));
                    cameraPosition(ii,1) =  cx + carPosition(ap,1);
                   cameraPosition(ii,2)  =  cy + carPosition(ap,2);
                
                   filmDistanceVec(ii)   = focusLens(lensFile,(max(cameraDistance(p(3)) + cameraDefocus(p(1)),0.1))*1000);
                
                   cameraDistanceVec(ii) = cameraDistance(p(3));
                end

%                 for cdef=1:length(cameraDefocus)
%                     for ch=1:length(cameraHeight)
%                         for cdd=1:length(cameraDistance)
%                             for co=1:length(cameraOrientation)
%                                 
%                                 % This would be the counter variable
%                                 loc = sub2ind(sz,cdef,ch,cdd,co);
%                                 
%                                 cx = cameraDistance(cdd)*sind(cameraOrientation(co));
%                                 cy = cameraDistance(cdd)*cosd(cameraOrientation(co));
%                                 
%                                 cameraPosition(loc,1) = cx + carPosition(ap,1);
%                                 cameraPosition(loc,2) = cy + carPosition(ap,2);
%                                 cameraPosition(loc,3) = cameraHeight(ch);
%                                 
%                                 filmDistanceVec(loc) = focusLens(lensFile,(max(cameraDistance(cdd)+cameraDefocus(cdef),0.1))*1000);
%                                 cameraDistanceVec(loc) = cameraDistance(cdd);
%                                 
%                             end
%                         end
%                     end
%                 end
                
                % This is building up the values used for the Conditions file.
                % values = rtbConditionsCreate(...)
                %
                % We might, again, do this with a ndgrid() on the lengths and do
                % a single loop all the way through.  These multiple nested
                % loops are very hard to understand and read.
                for ao=1:length(carOrientation);
                    for p=1:size(cameraPosition,1)
                        for s=1:size(shadowDirection,1)
                            for fn=1:length(fNumber)
                                for cpan=1:length(cameraPan)
                                    for ctilt=1:length(cameraTilt)
                                        for croll=1:length(cameraRoll)
                                            for df=1:length(diffraction)
                                                
                                                for mo=1:length(mode)
                                                    
                                                    
                                                    if strcmp(cameraType{lt},'pinhole')
                                                        currentFilmDistance = effectiveFocalLength(lensFile);
                                                    else
                                                        currentFilmDistance = filmDistanceVec(p);
                                                    end
                                                    
                                                    % ap is the car position
                                                    % index
                                                    cameraLookAt = [carPosition(ap,1:2) cameraHeight];
                                                    
                                                    fName = sprintf('%05i_city_%02i_car_%02i_%s_%s_%s_fN_%.2f_diff_%s_chr_%s',...
                                                        sceneID,cityId,carId,cameraType{lt},lensType{lt},mode{mo},fNumber(fn),diffraction{df},chrAber{df});
                                                    
                                                    values(cntr,1) = {fName};
                                                    values(cntr,2) = cameraType(lt);
                                                    values(cntr,3) = lensType(lt);
                                                    values(cntr,4) = mode(mo);
                                                    values(cntr,5) = num2cell(pixelSamples,1);
                                                    values(cntr,6) = num2cell(currentFilmDistance,1);
                                                    values(cntr,7) = num2cell(filmDiag,1);
                                                    values(cntr,8) = {mat2str(cameraPosition(p,:))};
                                                    values(cntr,9) = {mat2str(shadowDirection(s,:))};
                                                    values(cntr,10) = {mat2str(microlens{lt})};
                                                    
                                                    values(cntr,11) = {mat2str(cameraLookAt)};
                                                    
                                                    values(cntr,12) = num2cell(fNumber(fn),1);
                                                    values(cntr,13) = {mat2str(carPosition(ap,:))};
                                                    values(cntr,14) = num2cell(carOrientation(ao));
                                                    values(cntr,15) = {0};
                                                    values(cntr,16) = diffraction(df);
                                                    values(cntr,17) = chrAber(df);
                                                    values(cntr,18) = num2cell(cameraPan(cpan),1);
                                                    values(cntr,19) = num2cell(cameraTilt(ctilt),1);
                                                    values(cntr,20) = num2cell(cameraRoll(croll),1);
                                                    
                                                    cntr = cntr + 1;
                                                    
                                                    
                                                end
                                                
                                                sceneID = sceneID+1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                   
                % Here is names and the values for the conditions file.
                rtbWriteConditionsFile(conditionsFile,names,values);
                
                % Generate files and render
                % We parallelize scene generation, not the rendering because PBRT
                % automatically scales the number of processes to equal the nubmer of
                % cores.
                %
                nativeSceneFiles = rtbMakeSceneFiles(scene, 'hints', hints,...
                    'conditionsFile',conditionsFile);
                               
                radianceDataFiles = rtbBatchRender(nativeSceneFiles, 'hints', hints);
               
                % We aren't saving the radianceDataFiles for all the conditions.
                %  This means we have to rerun too many times.
                % load('radianceDataFiles');
                for i=1:length(radianceDataFiles)
                    % chdir(fullfile(nnGenRootPath,'local'));
                    % save('radianceDataFiles','radianceDataFiles');

                    radianceData = load(radianceDataFiles{i});
                    
                    %% Create an oi
                    oiParams.lensType = lensType{lt};
                    oiParams.filmDistance = 10;
                    oiParams.filmDiag = 20;
                    
                    [~, label] = fileparts(radianceDataFiles{i});
                                                            
                    oi = buildOi(radianceData.multispectralImage, [], oiParams);
                    oi = oiSet(oi,'name',label);
                    
                    ieAddObject(oi);
                    oiWindow;
                    
                end
            end
        end
    end
end

%% Save out the oi if you like
if 0
    chdir(fullfile(nnGenRootPath,'local'));
    oiNames = vcGetObjectNames('oi');
    for ii=1:length(oiNames)
        thisOI = ieGetObject('oi',ii);
        save(oiNames{ii},'thisOI');
    end
end



