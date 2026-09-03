function cfg = default_elfouhaily_ideal_curl_config()
%DEFAULT_ELFOUHAILY_IDEAL_CURL_CONFIG Configuration for a local 2-D curl.

cfg.randomSeed = 20260825;

cfg.domain.Lx = 32.0;
cfg.domain.Ly = 32.0;
cfg.domain.dx = 0.025;
cfg.domain.dy = 0.05;

cfg.sea.U10 = 5.0;
cfg.sea.inverseWaveAge = 0.84;
cfg.sea.windDirectionDeg = 0.0;
cfg.sea.targetHs = 0.35;

% A static surface supplies a geometric breaking-eligible crest, not a
% universal hydrodynamic breaking-onset prediction.
cfg.detection.propagationDirectionDeg = 0.0;
cfg.detection.smoothingLength = 0.25;
cfg.detection.heightQuantile = 0.88;
cfg.detection.slopeToleranceFactor = 0.32;
cfg.detection.curvatureQuantile = 0.70;
cfg.detection.forwardSlopeDistance = 0.80;
cfg.detection.edgeMargin = 4.0;
cfg.detection.refineRadius = 0.30;
cfg.detection.heightWeight = 0.50;
cfg.detection.curvatureWeight = 0.35;
cfg.detection.forwardSlopeWeight = 0.15;

% The material patch follows the detected crest ridge. Lengths are small
% relative to the 32 m sea and to the old 0.915 m pivot depth.
cfg.curl.crestHalfLength = 2.20;
cfg.curl.ridgeSearchHalfWidth = 0.55;
cfg.curl.ridgeSmoothSamples = 9;
cfg.curl.coreHalfWidth = 0.42;
cfg.curl.transitionWidth = 0.58;
cfg.curl.profileRear = -0.46;
cfg.curl.noseStart = 0.10;
cfg.curl.noseEnd = 0.42;
cfg.curl.lowerEnd = 0.86;
cfg.curl.upperAdvanceFraction = 1.45;
cfg.curl.noseRadiusFraction = 1.05;
cfg.curl.lowerAdvanceFactor = 0.92;
cfg.curl.lipTopDropFraction = 0.16;

cfg.curl.localHeightRadius = 1.25;

% Low-amplitude pre-shaping changes the neighboring noncurling surface too.
cfg.curl.crestLiftFraction = 0.015;
cfg.curl.evolutionLiftFraction = 0.005;
cfg.curl.forwardLeanFraction = 0.08;
cfg.curl.plateauFraction = 0.48;
cfg.curl.lowerBranchDropFraction = 0.82;
cfg.curl.rearSigmaFraction = 0.72;
cfg.curl.frontSigmaFraction = 0.42;
cfg.curl.evolutionSigmaFraction = 1.20;
cfg.curl.coreMaskThreshold = 0.08;

cfg.output.figureVisible = 'on';
cfg.output.saveMatFile = false;
cfg.output.outputDirectory = fullfile(fileparts(mfilename('fullpath')), ...
    'output');
end
