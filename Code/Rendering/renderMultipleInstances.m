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

cameras = nnGenCameras('type',{'pinhole'},...
                        'mode',{'radiance','mesh'},...
                        'distance',50);

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
hints.batchRenderStrategy.converter.remodelAfterMappingsFunction = @PBRTRemodellerV2;
hints.batchRenderStrategy.converter.rewriteMeshData = false;

resourceFolder = rtbWorkingFolder('folderName','resources',...
    'rendererSpecific',false,...
    'hints',hints);


% Copy resources
lensTypes = unique({cameras(:).lens});
lensFiles = fullfile(lensDir,strcat(lensTypes,'.dat'));
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

shadowDirection = [-0.5 -1 1];

carPosition = [-2.5 -2.5 0;
                1.5   25 0];
           
carOrientation = [0 90];

scene = cityScene;
for i=1:size(carPosition,1)
    objects(i).prefix = sprintf('car_inst_%i_',i); % Note that spaces or : are not allowed
    objects(i).position = carPosition(i,:);
    objects(i).orientation = carOrientation(i);
    objects(i).bndbox = mat2str(mexximpSceneBox(carScene));
    
    scene = mexximpCombineScenes(scene,carScene,...
        'insertTransform',mexximpTranslate([0 0 0]),...
        'cleanupTransform',mexximpTranslate([0 0 0]),...
        'insertPrefix',objects(i).prefix);
end

% Create a second arrangement, where the second car is moved by 10 meters.
objectArrangements = repmat({objects},[1 2]);
objectArrangements{2}(2).position = [1.5 15 0];


placedCameras = nnPlaceCameras(cameras,objectArrangements);



%% Create a list of render conditions
conditionsFile = fullfile(resourceFolder,'Conditions.txt');
names = cat(1,'imageName',fieldnames(placedCameras{1}),'objPosFile','shadowDirection');
values = cell(1,length(names));

cntr = 1;
sceneId=1;
for m=1:length(objectArrangements)
    objectArrangementFile = fullfile(resourceFolder,sprintf('Arrangement_%i.json',m));
    savejson('',objectArrangements{m},objectArrangementFile);
    
    currentCameras = placedCameras{m};
    
    for s=1:size(shadowDirection,1)
        for c=1:length(placedCameras{m});
            
            fName = sprintf('%03i_%s',sceneId,currentCameras(c).mode);
            
            values(cntr,1) = {fName};
            for i=2:(length(names)-2)
                values(cntr,i) = {currentCameras(c).(names{i})};
            end
            values(cntr,length(names)-1) = {objectArrangementFile};
            values(cntr,length(names)) = {shadowDirection(s,:)};
            
            if strcmp(currentCameras(c).mode,'radiance')
                sceneId = sceneId+1;
            end
            cntr = cntr + 1;
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
    oiParams.lensType = values{i,strcmp(names,'lens')};
    oiParams.filmDistance = values{i,strcmp(names,'filmDistance')};
    oiParams.filmDiag = values{i,strcmp(names,'filmDiagonal')};
    
    [~, label] = fileparts(radianceDataFiles{i});
        
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







