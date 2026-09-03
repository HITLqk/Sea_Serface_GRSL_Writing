function cfg = default_high_sea_state_config()
%DEFAULT_HIGH_SEA_STATE_CONFIG Configuration for 10-40 m/s validation.

cfg = default_thesis_two_group_config();
cfg.windSpeeds = (10:2:40)';
cfg.realizationSeeds = (1:20)';

% Above 13.3 m/s use the logarithmic Hu/TGRS branch.  Unlike Cox-Munk's
% unrestricted linear extrapolation, this branch has decreasing d(MSS)/dU.
cfg.empiricalTargetMode = 'hu-high-wind';
cfg.mssResidualFraction = 0.05;
cfg.requireNonlinearImprovement = true;

cfg.highWindTransition = 13.3;
cfg.tailSlopeStart = 24;
cfg.maximumTailToEarlySlopeRatio = 0.65;
cfg.minimumConcaveFraction = 0.65;
cfg.outputDirectory = fullfile(fileparts(mfilename('fullpath')), ...
    'output_high_sea_state');
cfg.figureVisible = 'off';
cfg.exportResolution = 450;
end
