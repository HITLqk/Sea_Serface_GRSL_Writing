function cfg = default_wind_components_config()
%DEFAULT_WIND_COMPONENTS_CONFIG Directional short wind-wave parameters.

rootDirectory = fileparts(mfilename('fullpath'));
cfg.backgroundMatFile = fullfile(rootDirectory,'..','NonLiner', ...
    'output','nonlinear_lie_elfouhaily_surface.mat');
cfg.randomSeed = 20260901;

% Resolved short-wave band. With dx=0.25 m, 1.5 m still has six samples.
cfg.wind.minimumWavelength = 1.5;
cfg.wind.maximumWavelength = 5.0;
cfg.wind.propagationDirectionDeg = 0.0;
cfg.wind.angularSpreadStdDeg = 8.0;
cfg.wind.radialTransitionFraction = 0.15;
cfg.wind.replacementEnergyRatio = 1.0;
cfg.wind.previewTimeStep = 0.20;
cfg.wind.gravity = 9.81;

% Crest detection prepares a target for the later Curl stage, but this
% module does not deform or overturn the surface.
cfg.detection.minimumHeightQuantile = 0.80;
cfg.detection.maximumAlongSlopeFraction = 0.25;
cfg.detection.heightWeight = 0.40;
cfg.detection.curvatureWeight = 0.35;
cfg.detection.componentWeight = 0.25;
cfg.detection.edgeMargin = 14.0;

cfg.output.figureVisible = 'on';
cfg.output.outputDirectory = fullfile(rootDirectory,'output');
cfg.output.saveSurfaceMat = true;
end
