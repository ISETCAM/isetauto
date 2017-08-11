function hints = nnHintsInit(varargin)
% Initialize the rendering hints with the defaults
%
% 
% We should set up input parameters for this
%
% BW/HB SCIEN Stanford, 2017

%%
p = inputParser;

% Some can be reset.
p.addParameter('imageWidth',640,@isnumeric);
p.addParameter('imageHeight',480,@isnumeric);
p.addParameter('recipeName',tempname,@ischar);
p.addParameter('tokenPath','',@(x)(exist(x,'file')));
p.addParameter('gcloud',false,@islogical);

% Local docker image
p.addParameter('dockerImage','',@ischar);

p.parse(varargin{:});

%%
hints.imageWidth  =  p.Results.imageWidth;
hints.imageHeight = p.Results.imageHeight;
hints.recipeName  = p.Results.recipeName;  % Name of the directory output
if p.Results.gcloud
    hints.renderer = 'PBRTCloud';          % We're only using PBRT right now
else
    hints.renderer = 'PBRT';               % We're only using PBRT right now
end
hints.copyResources = 1;                   % Is this a logical?? (BW)
hints.isParallel = false;

% Change the docker container
hints.batchRenderStrategy = RtbAssimpStrategy(hints);

hints.batchRenderStrategy.remodelPerConditionAfterFunction = @MexximpRemodellerMoveCar;
hints.batchRenderStrategy.converter = RtbAssimpPBRTConverter(hints);
hints.batchRenderStrategy.converter.remodelAfterMappingsFunction = @PBRTRemodeller;
hints.batchRenderStrategy.converter.rewriteMeshData = false;
if p.Results.gcloud
    % Google cloud run
    hints.tokenPath = p.Results.tokenPath; 
    hints.batchRenderStrategy.renderer = RtbPBRTCloudRenderer(hints);
    if isempty(p.Results.dockerImage)
        dockerImage = 'gcr.io/primal-surfer-140120/pbrt-v2-spectral-gcloud';
    end
    hints.batchRenderStrategy.renderer.pbrt.dockerImage = dockerImage;
    hints.batchRenderStrategy.renderer.cloudFolder = fullfile('gs://primal-surfer-140120.appspot.com',hints.recipeName);
else
    % Local run
    hints.batchRenderStrategy.renderer = RtbPBRTRenderer(hints);
    if isempty(p.Results.dockerImage)
        dockerImage = 'pbrt-v2-spectral';
    end
    hints.batchRenderStrategy.renderer.pbrt.dockerImage = dockerImage;
end

end