function cfg = default_breaker_morphology_validation_config()
%DEFAULT_BREAKER_MORPHOLOGY_VALIDATION_CONFIG Publication-style setup for
% breaker morphology validation.
%
% This configuration supports a compact GRSL-ready validation product:
%   1) one 1x2 scatter figure for joint morphology validation;
%   2) one small summary table exported as CSV and LaTeX.

patternsDirectory = fileparts(mfilename('fullpath'));
cfg.curlDirectory = fullfile(patternsDirectory,'..','Curl');
cfg.outputDirectory = fullfile(patternsDirectory,'output');
cfg.figureVisible = 'off';

% Use the same representative pre-impact curl family as the accepted Curl
% implementation. These values are not redrawn here; the Monte Carlo ranges
% below define the two compared groups.
cfg.curlOverrides.amplitudeCurl = 0.21;
cfg.curlOverrides.forwardGain = 1.05;
cfg.curlOverrides.verticalAngleRatio = 1.59;
cfg.curlOverrides.pivotDepth = 0.915;
cfg.curlOverrides.curlMultiplier = 1.16;

% Erinin et al. (JFM 2023), weak-to-strong plunging-breaker envelope.
cfg.reference.frontFaceAngleDeg = [65.0 70.0];
cfg.reference.rxOverLambda = [0.049 0.069];
cfg.reference.ryOverLambda = [0.045 0.059];

% The default normalization wavelength is the Elfouhaily spectral-peak
% wavelength. Set mode='fixed' to override.
cfg.normalization.mode = 'spectralPeak';
cfg.normalization.fixedWavelength = 1.1806;

cfg.extraction.centerlineHalfWidthCells = 0.55;
cfg.extraction.activeCurlFraction = 0.08;
cfg.extraction.frontForwardFraction = 0.35;
cfg.extraction.frontVerticalFraction = [0.15 0.85];

cfg.monteCarlo.randomSeed = 20260901;
cfg.monteCarlo.nPerGroup = 100;
cfg.monteCarlo.original = parameter_ranges( ...
    [0.17 0.23],[0.45 0.65],[0.80 1.10],[1.00 1.30],[0.20 0.40]);
cfg.monteCarlo.corrected = parameter_ranges( ...
    [0.19 0.23],[1.12 1.20],[0.85 0.98],[0.98 1.12],[1.54 1.64]);

% Output names.
cfg.output.figureBaseName = 'Fig_breaker_morphology_validation';
cfg.output.tableBaseName = 'Table_breaker_morphology_statistics';
cfg.output.exportPdf = true;
cfg.output.exportPng = true;
cfg.output.pngResolution = 450;

% Publication style.
cfg.style.figureWidthIn = 7.10;
cfg.style.figureHeightIn = 3.30;
cfg.style.fontName = 'Times New Roman';
cfg.style.fontSize = 8.5;
cfg.style.axisLineWidth = 0.75;
cfg.style.gridAlpha = 0.14;
cfg.style.envelopeFaceColor = [0.87 0.92 0.98];
cfg.style.envelopeEdgeColor = [0.45 0.56 0.73];
cfg.style.previousColor = [0.45 0.45 0.45];
cfg.style.proposedColor = [0.05 0.38 0.70];
cfg.style.previousMedianColor = [0.15 0.15 0.15];
cfg.style.proposedMedianColor = [0.85 0.33 0.10];
cfg.style.passColor = [0.10 0.55 0.25];
cfg.style.markerSize = 17;
cfg.style.medianMarkerSize = 40;
cfg.style.legendTokenSize = [16 8];
cfg.style.legendFontSize = 8.0;

cfg.labels.previous = 'Previous--Shallow-curl configuration';
cfg.labels.proposed = 'Proposed--Curling geometry';
cfg.labels.reference = 'Erinin et al.--Experimental envelope';
end

function ranges = parameter_ranges(amplitudeCurl,curlMultiplier, ...
    pivotDepth,forwardGain,verticalAngleRatio)
ranges.amplitudeCurl = amplitudeCurl;
ranges.curlMultiplier = curlMultiplier;
ranges.pivotDepth = pivotDepth;
ranges.forwardGain = forwardGain;
ranges.verticalAngleRatio = verticalAngleRatio;
end
