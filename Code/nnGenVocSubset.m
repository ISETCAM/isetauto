%% Create a sub-set from VOC dataset containing only selected classes

close all;
clear all;
clc;

mode = 'test';

classes = {'car'};

vocPath = fullfile('/','scratch','Datasets','PASCAL',mode,'VOC2007');
resPath = fullfile('/','scratch','Datasets','PASCAL',mode,sprintf('VOC2007%s',classes{1}));

datasetSubDirs = dir(vocPath);

for i=1:length(datasetSubDirs)
    if datasetSubDirs(i).isdir
        mkdir(fullfile(resPath,datasetSubDirs(i).name));
    end
    if strcmp(datasetSubDirs(i).name,'ImageSets')
        subDir = dir(fullfile(vocPath,'ImageSets'));
        for j=1:length(subDir)
            if subDir(j).isdir
                mkdir(fullfile(resPath,'ImageSets',subDir(j).name));
            end
        end
    end
end


sourceClassFile = fullfile(vocPath,'ImageSets','Main',sprintf('%s_%s.txt',classes{1},mode));
fid = fopen(sourceClassFile,'r');
data = textscan(fid,'%d %d');
fclose(fid);

destClassFile = fullfile(resPath,'ImageSets','Main',sprintf('%s_%s.txt',classes{1},mode));
fid = fopen(destClassFile,'w');

for i=1:length(data{1})
    if (data{2}(i) == -1)
        continue;
    end
    
    fileName = data{1}(i);
    
    srcFile = fullfile(vocPath,'Annotations',sprintf('%06i.xml',fileName));
    destFile = fullfile(resPath,'Annotations',sprintf('%06i.xml',fileName));
    copyfile(srcFile,destFile);
    
    srcFile = fullfile(vocPath,'JPEGImages',sprintf('%06i.jpg',fileName));
    destFile = fullfile(resPath,'JPEGImages',sprintf('%06i.jpg',fileName));
    copyfile(srcFile,destFile);
    
    fprintf(fid,'%06i %06i\n',fileName,data{2}(i));
    
end

fclose(fid);
fprintf('Done!\n');