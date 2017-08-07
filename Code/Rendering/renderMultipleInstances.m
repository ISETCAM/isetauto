%%% RenderToolbox3 Copyright (c) 2012-2013 The RenderToolbox3 Team.
%%% About Us://github.com/DavidBrainard/RenderToolbox3/wiki/About-Us
%%% RenderToolbox3 is released under the MIT License.  See LICENSE.txt.

%% Scene description

% Henryk Blasinski
close all;
clear all;
clc;

ieInit;
constants;

%% Simulation parameters

cameraType = {'pinhole'}; %'dgauss.22deg.12.5mm'
lensType = {'dgauss.22deg.6.0mm'};
mode = {'radiance','mesh'};

% Negative z is up.
% Scene is about 200x200m, units are mm.
% However we should specify meters, as they are automatically converted to
% mm in remodellers.

pixelSamples = 128;

shadowDirection = [-0.5 -1 1];

cameraDistance = [50];
cameraOrientation = [0];
cameraHeight = [-1.5];
cameraPTR = [0, 0, 0];
cameraDefocus = [0, -40];

diffraction = {'false'};
chromaticAberration = {'false'};

fNumber = 2.8;
filmDiag = [(1/6.4)*25.4];
microlensDim = [0, 0];


%% Choose renderer options.
hints.imageWidth = 640;
hints.imageHeight = 480;
hints.recipeName = 'MultiInstance'; % Name of the render
hints.renderer = 'PBRT'; % We're only using PBRT right now
hints.copyResources = 1;
hints.batchRenderStrategy = RtbAssimpStrategy(hints);

% Change the docker container
hints.batchRenderStrategy.renderer.pbrt.dockerImage = 'vistalab/pbrt-v2-spectral';
hints.batchRenderStrategy.remodelPerConditionAfterFunction = @MexximpRemodellerMultipleObj;
hints.batchRenderStrategy.converter.remodelAfterMappingsFunction = @PBRTRemodeller;
hints.batchRenderStrategy.converter.rewriteMeshData = false;

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


%% Choose files to render
cityId = 1;
carId = 1;


sceneFile = sprintf('City_%i.obj',cityId);
parentSceneFile = fullfile(assetDir,'City',sceneFile);
cityScene = mexximpCleanImport(parentSceneFile,...
    'ignoreRootTransform',true,...
    'flipUVs',true,...
    'imagemagicImage','hblasins/imagemagic-docker',...
    'toReplace',{'jpg','png','tga'},...
    'targetFormat','exr',...
    'makeLeftHanded',true,...
    'flipWindingOrder',true,...
    'workingFolder',resourceFolder);


carFile = sprintf('Car_%i.obj',carId);
parentSceneFile = fullfile(assetDir,car2directory{carId},carFile);
carScene = mexximpCleanImport(parentSceneFile,...
    'ignoreRootTransform',true,...
    'flipUVs',true,...
    'imagemagicImage','hblasins/imagemagic-docker',...
    'toReplace',{'jpg','png','tga'},...
    'targetFormat','exr',...
    'makeLeftHanded',true,...
    'flipWindingOrder',true,...
    'workingFolder',resourceFolder);

carPosition = [-2.5 -2.5 0;
                2    25 0];
           
 carOrientation = [0 90];

scene = cityScene;
for i=1:size(carPosition,1)
    objects(i).prefix = sprintf('car_inst_%i_',i); % Note that spaces or : are not allowed
    objects(i).position = mat2str(carPosition(i,:));
    objects(i).orientation = sprintf('%i',carOrientation(i));
    objects(i).bndbox = mat2str(mexximpSceneBox(carScene));
    
    scene = mexximpCombineScenes(scene,carScene,...
        'insertTransform',mexximpTranslate([0 0 0]),...
        'cleanupTransform',mexximpTranslate([0 0 0]),...
        'insertPrefix',objects(i).prefix);
end


conditionsFile = fullfile(resourceFolder,sprintf('Conditions_city_%i_car_%i.txt',cityId,carId));
objectMotionFile = fullfile(resourceFolder,'Movements.json');
savejson('',objects,objectMotionFile);

%% Create a list of render conditions
names = {'imageName','cameraType','lensType','mode','pixelSamples','filmDist','filmDiag','cameraPosition',...
    'cameraLookAt','cameraPTR','microlensDim','fNumber','diffraction','chromaticAberration',...
    'shadowDirection','objPosFile'};


values = cell(1,numel(names));
cntr = 1;

sceneId=1;
for ct=1:length(cameraType)
    lensFile = fullfile(lensDir,sprintf('%s.dat',lensType{ct}));
    
    [camPos, camLookAt, filmDist] = nnCameraPos(objects,cameraHeight,...
        cameraDistance,...
        cameraOrientation,...
        cameraDefocus,...
        lensFile);
    
    for p=1:size(camPos,1)
        for s=1:size(shadowDirection,1)
            for fn=1:length(fNumber)
                for pan=1:size(cameraPTR,1)
                    for mo=1:length(mode)
                        for df=1:length(diffraction)
                            
                            
                            if strcmp(cameraType{ct},'pinhole')
                                currentFilmDistance = effectiveFocalLength(lensFile);
                            else
                                currentFilmDistance = filmDist(p);
                            end
                            
                            fName = sprintf('%03i_city_%i_car_%i_%s_%s_%s_fN_%.2f',sceneId,cityId,carId,cameraType{ct},lensType{ct},mode{mo},fNumber(fn));
                            
                            values(cntr,1) = {fName};
                            values(cntr,2) = cameraType(ct);
                            values(cntr,3) = lensType(ct);
                            values(cntr,4) = mode(mo);
                            values(cntr,5) = num2cell(pixelSamples(ct),1);
                            values(cntr,6) = num2cell(currentFilmDistance,1);
                            values(cntr,7) = num2cell(filmDiag,1);
                            values(cntr,8) = {mat2str(camPos(p,:))};
                            values(cntr,9) = {mat2str(camLookAt(p,:))};
                            values(cntr,10) = {mat2str(cameraPTR(pan,:))};
                            values(cntr,11) = {mat2str(microlensDim)};
                            values(cntr,12) = num2cell(fNumber(fn),1);
                            values(cntr,13) = diffraction(df);
                            values(cntr,14) = chromaticAberration(df);
                            values(cntr,15) = {mat2str(shadowDirection(s,:))};
                            values(cntr,16) = {objectMotionFile};
                            
                            cntr = cntr + 1;
                            
                        end
                        sceneId = sceneId+1;
                    end
                end
            end
        end
    end
end

rtbWriteConditionsFile(conditionsFile,names,values);

% Generate files and render
% We parallelize scene generation, not the rendering because PBRT
% automatically scales the number of processes to equal the nubmer of
% cores.
%%


nativeSceneFiles = rtbMakeSceneFiles(scene, 'hints', hints,...
    'conditionsFile',conditionsFile);


%%
radianceDataFiles = rtbBatchRender(nativeSceneFiles, 'hints', hints);

%%

labelMap(1).name = 'car';
labelMap(1).id = 7;

for i=1:2:length(radianceDataFiles)
    radianceData = load(radianceDataFiles{i});
    
    
    %% Create an oi
    oiParams.lensType = lensType{1};
    oiParams.filmDistance = 10;
    oiParams.filmDiag = 20;
    
    [path, condname] = fileparts(radianceDataFiles{i});
    
    label = condname;
    
    oi = buildOi(radianceData.multispectralImage, [], oiParams);
    oi = oiSet(oi,'name',label);
    
    
    ieAddObject(oi);
    oiWindow;
    
    [classMap, instanceMap] = mergeMetadata(radianceDataFiles{i+1},labelMap);

    objects = getBndBox(classMap,instanceMap,labelMap);
    
    figure;
    imshow(oiGet(oi,'rgb image'));
    for j=1:length(objects)
       pos = [objects{j}.bndbox.xmin objects{j}.bndbox.ymin ...
              objects{j}.bndbox.xmax-objects{j}.bndbox.xmin ...
              objects{j}.bndbox.ymax-objects{j}.bndbox.ymin];
       rectangle('Position',pos,'EdgeColor','red'); 
    end
    
end







