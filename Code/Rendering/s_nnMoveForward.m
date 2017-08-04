%% Compare lens renderings of different cars
%
%  Set car and camera parameters
%  Set renderer options
%    There are some resources that are copied here.  Good if this could be
%    handled automatically
%  Make the combined scene
%    For each city and car (there is only one of each right now)
%     Import city scene
%     Import car scene
%     Combine the city and car scene
%  rtbWriteConditionsFile Combine parameters and write conditions file
%  rtbMakeSceneFiles - combined scene, conditions, hints
%  Run the batch renderer
%  Read the rendered radiance files into OI data structures
%
% TODO:
%    Illustrate impact of chromatic aberration and diffraction some how.
%    More colorful objects, scene complexity.
%    Save the geometry for efficiency, varying optics but not scene
%    Make a small video as the camera moves forward, say 120 frames comprising 4
%      secs?
%    Not sure we need to copoy both exr and jpg, but maybe.  Should check.
%
% Henryk Blasinski, SCIEN Stanford, 2017

%% Scene description

% Initialize ISET related variables
ieInit;

% Sets up related to the car renderings and local directory tree
% Maybe should be called nnDirectories
nnConstants;

% Small image size for debugging
% 640 x 480 when you care
hints = nnHintsInit('imageWidth',80,'imageHeight',60);

% Smaller for debugging
%% Simulation parameters
%
% Negative z is up. Scene is about 200x200m, units are mm. We  specify distances
% in meters; they are automatically converted to mm in remodellers.

% cameraType = {'pinhole','lens','lens'};
% lensType = {'tessar.22deg.6.0mm','tessar.22deg.6.0mm','2el.XXdeg.6.0mm'};
% microlens = {[0,0],[0,0],[0,0]};

cameraType = 'lens';
lensType = 'tessar.22deg.6.0mm';
microlens = [0,0];    % No microlens

mode = {'radiance'};

fNumber  = 2.8;
fov = 45;                  % Deg
meanIlluminance = 5;       % Lux

filmDiag = (1/3.6)*25.4;   % Millimeters

% Calculate the width from the diag assuming a 4x3 form factor
%
% filmDiag = sqrt(w^2 + h^2) = sqrt(w^2 + (0.75*w^2))
% w = filmDiag/sqrt(1 + (0.75)^2))
% w = filmDiag / sqrt( 1.5625) = filmDiag/1.25;
%
% We assume the angular (width) field of view is FOV deg.  The sensor size
% in mm can be calculated as
%
%    FOV/2 = atand( (width/2)/focalLength)
%    width = 2* focalLength*tand(FOV/2)
%
% So, if FOV = 45 and focalLength is 5mm, then the sensor size (width) is
%
%    2*5*tand(45/2), or 4.1421 mm
%
% The height is 0.75*width.  It is

% % Compute the horizontal field of view
% x = size(photons, 1);
% y = size(photons, 2);
% d = sqrt(x^2 + y^2);  % Number of samples along the diagonal
% fwidth= (oiParams.filmDiag / d) * x;    % Diagonal size by d gives us mm per step
%
% % multiplying by x gives us the horizontal mm
% % Calculate angle in degrees
% fov = 2 * atan2d(fwidth / 2, oiParams.filmDistance);
%
% % Store the horizontal field of view in degrees in the oi
% oi = oiSet(oi, 'fov', fov);

diffraction = {'true'};
chrAber     = {'true'};
pixelSamples    = 128;             % 256 is a lot.  192 is OK. Ray samples per pixel?
shadowDirection = [-0.5 -1 1;];

% Anything on the 'names' list can be a vector
cameraDistance    = [10 10.5 11];  % Meters
cameraOrientation = [0];
cameraPan     = 0;
cameraTilt    = 0;
cameraRoll    = 0;
cameraHeight  = -1.5;
cameraDefocus = 0;

carOrientation = [30];

maxCars = 1; maxCities = 1;

%% Choose renderer options
%
% Select the working folder and copy critical files there
% Not sure why the lens file is not handled automatically.  This caused a
% problem for me with the calibration rendering, too. (BW).

resourceFolder = rtbWorkingFolder('folderName','resources',...
    'rendererSpecific',false,...
    'hints',hints);

% Copy resources
lensFile = fullfile(lensDir,strcat(lensType,'.dat'));
copyfile(lensFile,resourceFolder);

% Copy sky map
skyFile = fullfile(assetDir,'City','*.exr');
copyfile(skyFile,resourceFolder);

% Use ISET, copy D65 spectrum
wave = 400:10:700;
d65  = ieReadSpectra('D65',wave);
d65 = 100*d65;

% [wave, d65] = rtbReadSpectrum(fullfile(rtbRoot,'RenderData','D65.spd'));
rtbWriteSpectrumFile(wave,d65,fullfile(resourceFolder,'D65.spd'));

%% Choose files to render
sceneID = 1;
cityId  = 1; 
carId   = 1;

% These are the variable names used in the conditionsFile.  See
%  https://github.com/RenderToolbox/RenderToolbox4/wiki/Conditions-File-
% Some of these are standard.  Some are selected here.
names = {'imageName','cameraType','lensType','mode','pixelSamples','filmDist','filmDiag','cameraPosition',...
    'shadowDirection','microlensDim','cameraLookAt','fNumber','carPosition','carOrientation','fog',...
    'diffraction','chromaticAberration','cameraPan','cameraTilt','cameraRoll'};
idxImageName = ismember('imageName',names);

sceneFile = sprintf('City_%i.obj',cityId);
parentSceneFile = fullfile(assetDir,'City',sceneFile);

fprintf('mexximpCleanImport cityScene\n');
[cityScene, ~] = mexximpCleanImport(parentSceneFile,...
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
% I don't understand this.
% nCarPositions  = 1;       % How do we want to control the car positions?
carPosition = zeros(1,3);
% for ii=2:nCarPositions
%     carPosition(ii,:) = drawCarPosition(cityId);
% end

% There seem to be 1 car and 1 city
carFile = sprintf('Car_%i.obj',carId);
parentSceneFile = fullfile(assetDir,car2directory{carId},carFile);

fprintf('mexximpCleanImport parentScene\n');
[carScene, elements] = mexximpCleanImport(parentSceneFile,...
    'ignoreRootTransform',true,...
    'flipUVs',true,...
    'imagemagicImage','hblasins/imagemagic-docker',...
    'toReplace',{'jpg','png','tga'},...
    'targetFormat','exr',...
    'makeLeftHanded',true,...
    'flipWindingOrder',true,...
    'workingFolder',resourceFolder);

disp('Combining city and car scenes');
scene = mexximpCombineScenes(cityScene,carScene,...
    'insertTransform',mexximpTranslate([0 0 0]),...
    'cleanupTransform',mexximpTranslate([0 0 0]));

lensFile = fullfile(lensDir,sprintf('%s.dat',lensType));

% Calculate the vectors of some of the parameters for the
% particular carPosition and lens.  Need better explanation.
[cameraPosition, filmDistanceVec, cameraDistanceVec] = ...
    nnCameraParams(carPosition(1,:),...
    cameraDefocus,cameraHeight,cameraDistance, cameraOrientation, ...
    lensFile);

% Current film distance
currentFilmDistance = filmDistanceVec(1);
   
cntr = 1;   % Counts the files

%% Make values used for the Conditions file.
%
% Parameters are placed in a struct that will be gridded for the
% conditions.

for ii=1:length(names), params.(names{ii}) = [];  end
fName = sprintf('move_%02i_%s_%s',carId,cameraType,lensType);

% Only list parameters that can be fully crossed by ndgrid
params.(names{1}) = fName;
params.(names{2}) = cameraType;
params.(names{3}) = lensType;
params.(names{4}) = mode;
params.(names{5}) = pixelSamples;
params.(names{6}) = currentFilmDistance;
params.(names{7}) = filmDiag;
params.(names{8}) = cameraPosition;
params.(names{9}) = shadowDirection;
params.(names{10}) = microlens;
params.(names{11}) = [carPosition(1,1:2) cameraHeight];  % cameraLookAt
params.(names{12}) = fNumber;
params.(names{13}) = carPosition;
params.(names{14}) = carOrientation;
params.(names{15}) = 0;
params.(names{16}) = diffraction;
params.(names{17}) = chrAber;
params.(names{18}) = cameraPan;
params.(names{19}) = cameraTilt;
params.(names{20}) = cameraRoll;

% This creates the first group of conditions
values = nnConditions(params);

% Now adjust the imageName, adding an integer index
for ii=1:size(values,1)
    values{ii,idxImageName} = sprintf('%02d-%s',cntr,values{ii,idxImageName});
    cntr = cntr+1;
end

% Write out the conditions file
conditionsFile = fullfile(resourceFolder,'Conditions.txt');
rtbWriteConditionsFile(conditionsFile,names,values);
% edit(conditionsFile);

%% Generate files and render
% We parallelize scene generation, not the rendering because
% PBRT automatically scales the number of processes to equal the
% number of cores.
%
nativeSceneFiles = rtbMakeSceneFiles(scene, 'hints', hints,...
    'conditionsFile',conditionsFile);

fprintf('Batch rendering %d files\n',length(nativeSceneFiles));
radianceDataFiles = rtbBatchRender(nativeSceneFiles, 'hints', hints);

% We aren't saving the radianceDataFiles for all the conditions.
% This means we have to rerun too many times.
%
% Also, we don't have the true irradiance level, just a
% noise-free irradiance.  So, we should aim to set the
% irradiance to a reasonable level here.
%
% load('radianceDataFiles');
fprintf('Creating OI\n');
for i=1:length(radianceDataFiles)
    % chdir(fullfile(nnGenRootPath,'local'));
    % save('radianceDataFiles','radianceDataFiles');
    
    radianceData = load(radianceDataFiles{i});
    
    % Create an oi and set the parameters
    clear oiParams;
    oiParams.optics_name = lensType;
    oiParams.optics_model = 'diffractionlimited';
    oiParams.fov = fov;
    switch ieParamFormat(lensType)
        case 'pinhole'
            oiParams.optics_fnumber = 999;
        otherwise
            oiParams.optics_fnumber = fNumber;
    end
    oiParams.optics_focalLength = filmDistanceVec(1)*1e-3; % In meters
    [~, label] = fileparts(radianceDataFiles{i});
    oiParams.name = label;
    
    oi = buildOi(radianceData.multispectralImage, [], oiParams);
    
    oi = oiAdjustIlluminance(oi,meanIlluminance);
    
    ieAddObject(oi);
    oiWindow;
    
end




%% Save out the oi if you like
if 0
    chdir(fullfile(nnGenRootPath,'local','tmp'));
    oiNames = vcGetObjectNames('oi');
    for ii=1:length(oiNames)
        thisOI = ieGetObject('oi',ii);
        save([oiNames{ii},'.mat'],'thisOI');
    end
end

%%
if 0
    %% Experiment with different camera renderings
    oi   = ieGetObject('oi');
    fov  = oiGet(oi,'fov');
    oi   = oiAdjustIlluminance(oi,10);   % The illuminance values are very small
    
    % Big sensor
    sensor = sensorCreate;
    sensor = sensorSet(sensor,'fov',fov);
    
    sensor = sensorCompute(sensor,oi);
    ieAddObject(sensor); sensorWindow;
    
    ip = ipCreate;
    ip = ipCompute(ip,sensor);
    ieAddObject(ip); ipWindow;
    
end

%%




