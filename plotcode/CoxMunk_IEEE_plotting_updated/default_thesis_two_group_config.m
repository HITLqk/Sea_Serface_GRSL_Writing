function cfg = default_thesis_two_group_config()
%DEFAULT_THESIS_TWO_GROUP_CONFIG Configuration for the clean two-group test.

cfg.windSpeeds = (1:10)';
cfg.realizationSeeds = (1:20)';
cfg.inverseWaveAge = 0.84;
cfg.windDirectionDeg = 0;

% The primary-wave grid uses dk = kp/12, satisfying the thesis condition
% dk <= kp/10. Its nonlinear input/output bands prevent quadratic aliasing.
cfg.primaryGridSize = 256;
cfg.primaryPeakSamples = 12;
cfg.primaryMaximumPeakMultiple = 8;
cfg.lieInputPeakMultiple = 4;
cfg.lieOutputPeakMultiple = 8;

% Independent octave tiles represent the short-wave contribution without
% requiring one prohibitively large millimetre-resolution domain.
cfg.shortWaveTileSize = 64;
cfg.shortWaveModesBelowBand = 8;
cfg.maximumOpticalWavenumber = pi*1000;

% Dimensionless strength of the second-order Lie/Creamer correction.  U10 is
% deliberately NOT multiplied into the quadratic term (that is dimensionally
% inconsistent); wind dependence already enters through the input spectrum.
cfg.modifiedLieScale = 1.0;

% Moment dressing.  The empirical target is the component-wise mean of the
% available Cox-Munk, Hu/TGRS and Guerin laws.  residualFraction=0 gives an
% exact MSS match; a small positive value retains realistic realization
% scatter while contracting the log-error toward the empirical target.
cfg.enableEmpiricalMssDressing = true;
cfg.mssResidualFraction = 0.15;
cfg.empiricalTargetMode = 'fusion';
cfg.requireNonlinearImprovement = true;

cfg.outputDirectory = fullfile(fileparts(mfilename('fullpath')),'output');
cfg.figureVisible = 'off';
cfg.exportResolution = 450;
end
