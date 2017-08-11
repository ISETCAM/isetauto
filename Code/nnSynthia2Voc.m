% Convert Synthia dataset to VOC format
%
% HB 2017

close all;
clear all;
clc;

nTest = 629;
nTrain = 1515;

synthiaPath = fullfile('~','Downloads','SYNTHIA-RAND-CITYSCAPES');

labelMap(1).name = 'car';
labelMap(1).id = 8;

currentSet = 'RAND';
destDir = fullfile('/','Users','hblasins','Downloads','SYNTHIA-VOC');
if ~exist(destDir,'dir')
    mkdir(destDir);
end

xVal = {'trainval','test'};
for m=1:length(xVal)
    fName = fullfile(destDir,xVal{m},currentSet,'Annotations');
    mkdir(fName);
    fName = fullfile(destDir,xVal{m},currentSet,'ImageSets','Layout');
    mkdir(fName);
    fName = fullfile(destDir,xVal{m},currentSet,'ImageSets','Main');
    mkdir(fName);
    fName = fullfile(destDir,xVal{m},currentSet,'ImageSets','Segmentation');
    mkdir(fName);
    fName = fullfile(destDir,xVal{m},currentSet,'JPEGImages');
    mkdir(fName);
end

% Prepare files listing which image contains which class
fids = cell(length(xVal),length(labelMap)); 
for x=1:length(xVal)
for v=1:length(labelMap)

    fName = fullfile(destDir,xVal{x},currentSet,'ImageSets','Main',sprintf('%s_%s.txt',lower(labelMap(v).name),xVal{x}));
    fids{x,v} = fopen(fName,'w');
    
end
end



files = dir(fullfile(synthiaPath,'RGB','*.png'));
ids = randperm(length(files));

cntr = 1;
i=1;
while cntr<=(nTest+nTrain)
    if cntr<=nTest
        mode = 'test';
    else
        mode = 'trainval';
    end
    
    outputFileName = sprintf('%06i',cntr);
    outputJpegFileName = sprintf('%06i.png',cntr);
    outputXmlFileName = sprintf('%06i.xml',cntr);
    
    imageName = fullfile(synthiaPath,'RGB',files(ids(i)).name);
    labelsName = fullfile(synthiaPath,'GT','LABELS',files(ids(i)).name);
    i=i+1;
    
    image = imread(imageName);
    labels = imread(labelsName);
    
    objectLabels = labels(:,:,1);
    objectInstances = labels(:,:,2);
    
    objects = getBndBox(objectLabels, objectInstances, labelMap);
    
    if isempty(objects)
        continue;
    end
    
    %{
    figure;
    imshow(image);
    for j=1:length(objects)
       pos = [objects{j}.bndbox.xmin objects{j}.bndbox.ymin ...
              objects{j}.bndbox.xmax-objects{j}.bndbox.xmin ...
              objects{j}.bndbox.ymax-objects{j}.bndbox.ymin];
       rectangle('Position',pos,'EdgeColor','red'); 
    end
    %}
    
    annotation.folder = mode;
    annotation.filename = sprintf('%06i.png',cntr);
    annotation.source.annotation = fullfile('GT','LABELS',files(ids(i)).name);
    annotation.source.database = 'SYNTHIA-RAND';
    annotation.source.image = files(ids(i)).name;
    
    annotation.size.depth = size(image,3);
    annotation.size.height = size(image,1);
    annotation.size.width = size(image,2);
    
    annotation.object = objects(:);
    
    imwrite(image,fullfile(destDir,mode,currentSet,'JPEGImages',outputJpegFileName));
    s.annotation = annotation;
    struct2xml(s,fullfile(destDir,mode,currentSet,'Annotations',outputXmlFileName));
    
    sel = cellfun(@(x) strcmp(x,mode),xVal);
    for o=1:length(annotation.object)
        for c=1:length(labelMap)
            if strcmpi(labelMap(c).name,annotation.object{o}.name)
                isPresent = strcmpi(annotation.object{o}.name,labelMap(c).name)*2-1;
                fprintf(fids{sel,c},'%s %i\n',outputFileName,isPresent);
            end
        end
    end
    
    
    cntr = cntr+1;
end

% Close files
for i=1:numel(fids)
    fclose(fids{i});
end
