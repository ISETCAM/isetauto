% This script reads in the sensor irradiance data generated with RTB4+PBRT,
% simulates the images some camera would produce and arranges images and labels
% in a format that matches that of PASCAL VOC datasets.
%
% Copytight, Henryk Blasinski 2017.

close all;
clear all;
clc;

ieInit;

numTestImages = 629;
numImages = 2144;
recipe = 'Car-Complete-Pinhole';

lightLevels = [0.1, 1, 10, 100, 1000, 10000];
expTime = 0.015;
subMode = 'rawRGB';

%These class ids correspond to the ones from PASCAL VOC
labelMap(1).name = 'car';
labelMap(1).id = 7;

for ll=1:length(lightLevels)
    
    mode = sprintf('%s_luxLevel_%.1f',subMode,lightLevels(ll));
    
    
    dataDir = fullfile('/','share','wandell','data','NN_Camera_Generalization','Renderings',recipe);
    renderDir = fullfile('renderings','PBRTCloud');
    
    destDir = fullfile('/','scratch','Datasets',recipe);
    
    % Prepare the directory structure
    xVal = {'trainval','test'};
    for i=1:length(xVal)
        dirName = fullfile(destDir,xVal{i},mode,'Annotations');
        if exist(dirName,'dir') == false, mkdir(dirName); end;
        dirName = fullfile(destDir,xVal{i},mode,'JPEGImages');
        if exist(dirName,'dir') == false, mkdir(dirName); end;
        dirName = fullfile(destDir,xVal{i},mode,'ImageSets','Main');
        if exist(dirName,'dir') == false, mkdir(dirName); end;
    end
    
    % Prepare files listing which image contains which class
    fids = cell(length(xVal),length(labelMap));
    for x=1:length(xVal)
        for v=1:length(labelMap)
            
            fName = fullfile(destDir,xVal{x},mode,'ImageSets','Main',sprintf('%s_%s.txt',lower(labelMap(v).name),xVal{x}));
            fids{x,v} = fopen(fName,'w');
            
        end
    end
    
    % Create a label map file
    fName = fullfile(destDir,sprintf('%s_label_map.pbtxt',recipe));
    fid = fopen(fName,'w');
    fprintf(fid,'item {\n   id: 0\n   name: ''none_of_the_above''\n}\n\n');
    
    for i=1:length(labelMap)
        fprintf(fid,'item {\n   id: %i\n   name: ''%s''\n}\n\n',labelMap(i).id,lower(labelMap(i).name));
    end
    fclose(fid);
    
    %%
    
    fileNames = dir(fullfile(dataDir,renderDir,'*radiance*.mat'));
    nFiles = length(fileNames);
        
    rng(1);
    shuffling = randperm(nFiles);
    
    cntr = 1;
    for f=1:nFiles
        
        outputFileName = sprintf('%06i',cntr);
        outputXmlFileName = sprintf('%s.xml',outputFileName);
        outputJpegFileName = sprintf('%s.png',outputFileName);
        
        if cntr <= numTestImages;
            % Test set
            currentSet = 'test';
        else
            currentSet = 'trainval';
        end
        if cntr > numImages;
            break;
        end
        outputJpegFullFileName = fullfile(destDir,currentSet,mode,'JPEGImages',outputJpegFileName);
        outputXmlFullFileName = fullfile(destDir,currentSet,mode,'Annotations',outputXmlFileName);
        
        if exist(outputJpegFullFileName,'file') && exist(outputXmlFullFileName,'file')
            fprintf('File: %s exists, skipping\n',outputJpegFullFileName);
            cntr = cntr + 1;
            continue;
        end
        
        
        
        inputFileName = fileNames(shuffling(f)).name;
        
        [pth, name] = fileparts(inputFileName);
        
        %% Load image radiance data
        radianceDataFileName = fullfile(dataDir,renderDir,inputFileName);
        
        radianceData = load(radianceDataFileName);
        
        % Create an oi
        oiParams.lensType = 'pinhole';
        oiParams.filmDistance = 10;
        oiParams.filmDiag = 20;
        
        
        
        oi = BuildOI(radianceData.multispectralImage, [], oiParams);
        oi = oiSet(oi,'name',name);
        oi = oiAdjustIlluminance(oi,lightLevels(ll),'mean');
        
        ieAddObject(oi);
        % oiWindow();
        
        sensor = sensorCreate('bayer (rggb)');
        sensor = sensorSet(sensor,'name',name);
        sensor = sensorSet(sensor,'size',oiGet(oi,'size'));
        sensor = sensorSet(sensor,'pixel widthandheight',[oiGet(oi,'hres'), oiGet(oi,'wres')]);
        sensor = sensorSet(sensor,'analog gain',1);
        sensor = sensorSet(sensor,'exposure time',expTime);
        sensor = sensorSet(sensor,'quantizationmethod','8 bit');
        
        
        sensor = sensorCompute(sensor,oi);
        ieAddObject(sensor);
        % sensorWindow();
        
        
        
        ip = ipCreate();
        ip = ipSet(ip,'name',name);
        ip = ipCompute(ip,sensor);
        ieAddObject(ip);
        % ipWindow();
        
        switch subMode
            case 'sRGB'
                img = ipGet(ip,'data srgb');
            case 'fullResRGB'
                img = oiGet(oi,'rgb image');
            case 'linearRGB'
                img = uint8(ipGet(ip,'sensor channels'));
            case 'rawRGB'
                img = uint8(ipGet(ip,'sensor mosaic'));
                img = repmat(img,[1 1 3]);
        end
        
        img = uint8(255*double(img)/max(double(img(:))));
        
        
        %% Labels
        
        meshDataFileName = fullfile(dataDir,renderDir,sprintf('%s.mat',strrep(name,'radiance','mesh')));
        [labels, instances] = mergeMetadata(meshDataFileName,labelMap);
        
        objects = getBndBox(labels, instances, labelMap);
        
        annotation.folder = mode;
        annotation.filename = outputJpegFileName;
        annotation.source.annotation = sprintf('%s.mat',strrep(name,'radiance','mesh'));
        annotation.source.database = recipe;
        annotation.source.image = inputFileName;
        
        annotation.size.depth = size(img,3);
        annotation.size.height = size(img,1);
        annotation.size.width = size(img,2);
        
        annotation.object = objects(:);
        
        
        if isempty(objects)
            fprintf('No Objects in the image, skipping\n');
            continue;
        end
        
        
        
        
        %% Save data
        
        
        imwrite(img,outputJpegFullFileName);
        s.annotation = annotation;
        struct2xml(s,outputXmlFullFileName);
        
        sel = cellfun(@(x) strcmp(x,currentSet),xVal);
        
        for o=1:length(annotation.object)
            for c=1:length(labelMap)
                if strcmpi(labelMap(c).name,annotation.object{o}.name)
                    isPresent = strcmpi(annotation.object{o}.name,labelMap(c).name)*2-1;
                    fprintf(fids{sel,c},'%s %i\n',outputFileName,isPresent);
                end
            end
        end
        
        
        cntr = cntr + 1;
    end
    
    % Close files
    for i=1:numel(fids)
        fclose(fids{i});
    end
    
end



