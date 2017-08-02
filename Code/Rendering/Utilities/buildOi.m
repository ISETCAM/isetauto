function oi = buildOi(photons, depthMap, oiParams)
% Build and ISET Optical Image
% Required parameters
%   photons  - h x w x 31 multi-spectral radiance matrix
%   depthMap - h x w depth matrix
%   oiParams - struct oiSet() params
%
% HB, Scien Team, 2017

%% New optical image with photon and depth data.
oi = oiCreate;

oi = initDefaultSpectrum(oi);

oi = oiSet(oi, 'photons', single(photons));
oi = oiSet(oi, 'depthmap', depthMap);

%% This should be only parameters consistent with oiSet()

% And the parameters should be one at a time, not conjoined as below.
% It looks like we set optics name, optics fnumber, optics focal length (m)
%    and fov
%
 names = fieldnames(oiParams);
 for ii=1:length(names)
    oi = oiSet(oi,names{ii},oiParams.(names{ii}));
 end


%% Transfer various lens parameters to optical image.
% if strcmp(oiParams.lensType, 'realisticDiffraction')
%     oi = oiSet(oi, 'optics name', oiParams.specFile);
%     oi = oiSet(oi, 'optics fnumber', oiParams.filmDistance / oiParams.apertureDiameter);
%     
% elseif strcmp(oiParams.lensType, 'PinholeLens') % TODO: What are the real names for this in PBRT?
%     % Pinholes have no real aperture size.  So, we set the f-number
%     % really big.
%     oi = oiSet(oi, 'optics name', 'pinhole');
%     oi = oiSet(oi, 'optics fnumber', 999);
%     
% elseif strcmp(oiParams.lensType, 'IdealLens') % TODO: What are the real names for this in PBRT?
%     % This case is a diffraction limited lens but with an aperture of a
%     % real size.
%     oi = oiSet(oi, 'optics name', 'diffraction limited');
%     oi = oiSet(oi, 'optics fnumber', oiParams.filmDistance / oiParams.apertureDiameter);
%     
% else
%     fprintf('Could not find lens type in hints! Setting arbitrary values.\n')
%     TODO: What sort of arbitrary values should go here?
%     
% end
% 
% oi = oiSet(oi, 'optics focal length', oiParams.filmDistance * 1e-3);

% Funny, I computed this elsewhere for a 480 x 640 (height, weight) parameter.
% Same formula.  We should just send in the wAngular parameter, computing this
% in the main routine

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
