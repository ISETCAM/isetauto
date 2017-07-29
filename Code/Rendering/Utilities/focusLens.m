function [ sensorDistance ] = focusLens( lensFileName, distance )
% Compute the distance from the lens to the sensor
%
%   sensorDistance = focusLens( lensFileName, distance )
%
% This function uses CISET autofocus function to compute the distance
% between the lens and the sensor (sensorDistance) so that the point at
% distance away is in focus
%
% See also: t_autofocus.m (in ciset)
%
% Copyright, Henryk Blasinski 2017

%%  Initialize a point and a camera

point{1} = [0 0 -distance];

% This human eye model has a focal length of 16.5 mm, which we confirm when
% running ray trace in PBRT and ray trace in CISET. See -
lens = lensC('fileName',lensFileName);

film = filmC;

camera = psfCameraC('lens',lens,'film',film,'pointsource',point);

%%  Find the film focal length for this wavelength

% Call autofocus, setting the indices of refraction of air and water
camera.autofocus(550,'nm',1,1);

% Show adjusted position for focus
sensorDistance = camera.film.position(3);


end

