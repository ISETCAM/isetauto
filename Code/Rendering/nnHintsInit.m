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
p.parse(varargin{:});

%%
hints.imageWidth =  p.Results.imageWidth;
hints.imageHeight = p.Results.imageHeight;
hints.recipeName = 'Car-Different-Lenses'; % Name of the render
hints.renderer = 'PBRT'; % We're only using PBRT right now
hints.copyResources = 1;
hints.isParallel = false;

% Change the docker container
hints.batchRenderStrategy = RtbAssimpStrategy(hints);

hints.batchRenderStrategy.remodelPerConditionAfterFunction = @MexximpRemodellerMoveCar;
hints.batchRenderStrategy.converter = RtbAssimpPBRTConverter(hints);
hints.batchRenderStrategy.converter.remodelAfterMappingsFunction = @PBRTRemodeller;
hints.batchRenderStrategy.converter.rewriteMeshData = false;
hints.batchRenderStrategy.renderer = RtbPBRTRenderer(hints);
hints.batchRenderStrategy.renderer.pbrt.dockerImage = 'vistalab/pbrt-v2-spectral';


end