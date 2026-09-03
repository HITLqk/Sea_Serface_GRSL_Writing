function cfg = default_nonlinear_lie_config()
%DEFAULT_NONLINEAR_LIE_CONFIG Configuration for nonlinear Elfouhaily sea.

cfg.randomSeed = 20260825;

cfg.domain.Lx = 128.0;
cfg.domain.Ly = 128.0;
cfg.domain.dx = 0.25;
cfg.domain.dy = 0.25;

cfg.sea.U10 = 10.0;
cfg.sea.inverseWaveAge = 0.84;
cfg.sea.windDirectionDeg = 0.0;
% Optional physical synthesis band. The defaults preserve the original
% full resolved Elfouhaily spectrum. Validation runs can set finite limits
% to match an instrument's observable wave scales.
cfg.sea.minimumWavenumber = 0.0;
cfg.sea.maximumWavenumber = inf;

% Second-order Lie/Creamer-type bound-wave correction. The Elfouhaily
% spectrum already carries the wind dependence, so the default extra gain
% is unity and wind speed is never inserted into a wavenumber operator.
cfg.lie.baseGain = 1.0;
cfg.lie.referenceWindSpeed = 10.0;
cfg.lie.windExponent = 0.0;
cfg.lie.minimumWindGain = 1.0;
cfg.lie.maximumWindGain = 1.0;
cfg.lie.directionStrength = 0.35;
cfg.lie.directionExponent = 2.0;

% Fixed physical bandwidth for the quadratic bound-wave operation. The
% input contains wavelengths >= about 1 m; its second harmonics remain
% below the 0.25 m grid Nyquist limit. A 3/2 padded product removes aliasing.
cfg.lie.nonlinearInputCutoff = 6.0;
cfg.lie.nonlinearOutputCutoff = 12.0;
cfg.lie.filterTransitionFraction = 0.15;
cfg.lie.dealiasExpansion = 1.5;

% Horizontal Riesz displacement sharpens crests. Its amplitude is reduced
% automatically until the horizontal mapping remains one-to-one.
cfg.lie.horizontalDisplacement = 0.75;
cfg.lie.minimumJacobian = 0.20;
cfg.lie.backoffFactor = 0.85;
cfg.lie.maximumBackoffSteps = 30;

cfg.output.figureVisible = 'on';
cfg.output.outputDirectory = fullfile(fileparts(mfilename('fullpath')), ...
    'output');
cfg.output.saveSurfaceMat = true;
end
