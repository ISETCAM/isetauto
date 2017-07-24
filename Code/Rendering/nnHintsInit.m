function hints = nnHintsInit()
% Initialize the rendering hints with the defaults
%
% 
% We should set up input parameters for this
%
% BW/HB SCIEN Stanford, 2017

hints.imageWidth = 640;
hints.imageHeight = 480;
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