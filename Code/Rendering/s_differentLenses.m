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
% Henryk Blasinski, SCIEN Stanford, 2017

%% Scene description

% Initialize ISET related variables
ieInit;

% Sets up related to the car renderings and local directory tree
nnConstants;

% Small image size for debugging
hints = nnHintsInit('imageWidth',160,'imageHeight',120);

% Smaller for debugging
%% Simulation parameters
%
% Negative z is up. Scene is about 200x200m, units are mm. We  specify distances
% in meters; they are automatically converted to mm in remodellers.

% cameraType = {'pinhole','lens','lens'};
% lensType = {'tessar.22deg.6.0mm','tessar.22deg.6.0mm','2el.XXdeg.6.0mm'};
% microlens = {[0,0],[0,0],[0,0]};

cameraType = {'pinhole'};
lensType = {'tessar.22deg.6.0mm'};
microlens = {[0,0]};

mode = {'radiance'};

fNumber  = 2.8;
filmDiag = (1/3.6)*25.4;   % Millimeters
fov = 45;                  % Deg

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


diffraction = {'false','true'};
chrAber = {'false','true'};

% diffraction = {'false','true'};
% chrAber = {'false','true'};

pixelSamples = 128;        % Ray samples per pixel?
shadowDirection = [-0.5 -1 1;];

cameraDistance = [10 20];
cameraOrientation = [0];
cameraPan = [0];
cameraTilt = [0];
cameraRoll = [0];

cameraHeight = -1.5; cameraDefocus = 0;

nCarPositions = 1;  carOrientation = [30];

maxCars = 1; maxCities = 1;


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

% Use ISET, copy D65 spectrum
wave = 400:10:700;
d65  = ieReadSpectra('D65',wave);
d65 = 100*d65;

% [wave, d65] = rtbReadSpectrum(fullfile(rtbRoot,'RenderData','D65.spd'));
rtbWriteSpectrumFile(wave,d65,fullfile(resourceFolder,'D65.spd'));


%% Choose files to render
sceneID = 1;

% These are the variable names used in the conditionsFile.  See
%  https://github.com/RenderToolbox/RenderToolbox4/wiki/Conditions-File-
% Some of these are standard.  Some are selected here.
names = {'imageName','cameraType','lensType','mode','pixelSamples','filmDist','filmDiag','cameraPosition',...
    'shadowDirection','microlensDim','cameraLookAt','fNumber','carPosition','carOrientation','fog',...
    'diffraction','chromaticAberration','cameraPan','cameraTilt','cameraRoll'};

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
        
        disp('Combining city and car scenes');
        scene = mexximpCombineScenes(cityScene,carScene,...
            'insertTransform',mexximpTranslate([0 0 0]),...
            'cleanupTransform',mexximpTranslate([0 0 0]));
        
        
        for ap=1:nCarPositions;
            fprintf('Car position %d\n',ap);
            for lt=1:length(lensType)
                fprintf('Lens type %s\n',lensType{lt});
                
                % It seems like we get a clean copy of the Conditions.txt file
                % here?  Later, we write it out.
                conditionsFile = fullfile(resourceFolder,'Conditions.txt');
                
                % The conditions file has a structure of names by values.  We
                % set the size here.
                values = cell(1,numel(names));
                cntr = 1;
                
                % Here we make the vectors of some of the parameters for the
                % particular carPosition and lens.
                lensFile = fullfile(lensDir,sprintf('%s.dat',lensType{lt}));
                [cameraPosition, filmDistanceVec, cameraDistanceVec] = ...
                    nnCameraParams(carPosition(ap,:),...
                    cameraDefocus,cameraHeight,cameraDistance, cameraOrientation, ...
                    lensFile);
                
                % Here we make the values used for the Conditions file. This
                % could be a function
                % 
                %   values = rtbConditionsCreate(...)
                %
                % The function might be a ndgrid() call that uses the lengths() of the
                % variables. These multiple nested  
                % loops are very hard to understand and read.
                %
                % Maybe the rule is this.  There are multiple params, which is a
                % cellarray.  Each element of the cell aray is either a string
                % or a vector, or a matrix.
                %
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
                
                fprintf('Back rendering %d files\n',length(nativeSceneFiles));
                radianceDataFiles = rtbBatchRender(nativeSceneFiles, 'hints', hints);
                
                % We aren't saving the radianceDataFiles for all the conditions.
                % This means we have to rerun too many times.
                %
                % Also, we don't have the true irradiance level, just a
                % noise-free irradiance.  So, we should aim to set the
                % irradiance to a reasonable level here.
                %
                % load('radianceDataFiles');
                for i=1:length(radianceDataFiles)
                    % chdir(fullfile(nnGenRootPath,'local'));
                    % save('radianceDataFiles','radianceDataFiles');
                    
                    radianceData = load(radianceDataFiles{i});
                    
                    % Create an oi and set the parameters
                    clear oiParams;
                    oiParams.optics_name = lensType{lt};
                    oiParams.fov = fov;
                    switch ieParamFormat(lensType{lt})
                        case 'pinhole'
                            oiParams.optics_fnumber = 999;
                        otherwise
                            oiParams.optics_fnumber = fNumber(i);
                    end
                    oiParams.optics_focalLength = filmDistanceVec(i)*1e-3; % In meters
                    [~, label] = fileparts(radianceDataFiles{i});
                    oiParams.name = label;
                    
                    oi = buildOi(radianceData.multispectralImage, [], oiParams);

                    ieAddObject(oi);
                    oiWindow;
                    
                end
            end
        end
    end
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




