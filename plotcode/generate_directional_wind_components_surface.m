function surfaceData = generate_directional_wind_components_surface(cfg)
%GENERATE_DIRECTIONAL_WIND_COMPONENTS_SURFACE Replace a short-wave band.
%   The nonlinear Elfouhaily background is separated into a target short
%   band and its residual. The target band is replaced, not duplicated, by
%   an equal-energy random wind-wave field with a prescribed traveling
%   direction. A steep crest candidate is reported for the later Curl
%   module; no curl deformation is applied here.

arguments
    cfg (1,1) struct = default_wind_components_config()
end

validate_config(cfg);
background = load_background(cfg.backgroundMatFile);
rng(cfg.randomSeed,'twister');

[Ny,Nx] = size(background.Z);
dx = infer_spacing(background.XLinear,2);
dy = infer_spacing(background.YLinear,1);
Lx = Nx*dx;
Ly = Ny*dy;
dkx = 2*pi/Lx;
dky = 2*pi/Ly;
kx = ifftshift((-Nx/2:Nx/2-1)*dkx);
ky = ifftshift((-Ny/2:Ny/2-1)*dky);
[KX,KY] = meshgrid(kx,ky);
K = hypot(KX,KY);
Ksafe = max(K,eps);

kLow = 2*pi/cfg.wind.maximumWavelength;
kHigh = 2*pi/cfg.wind.minimumWavelength;
nyquistRadius = min(pi/dx,pi/dy);
assert(kHigh < 0.90*nyquistRadius, ...
    ['The requested minimum wavelength is too close to the grid ', ...
    'Nyquist limit. Increase minimumWavelength or reduce dx/dy.']);

bandWindow = smooth_radial_bandpass(K,kLow,kHigh, ...
    cfg.wind.radialTransitionFraction);
backgroundHat = fft2(background.Z);
removedBand = real(ifft2(backgroundHat.*bandWindow));
backgroundResidual = background.Z-removedBand;

[componentHat0,spectrumShape,componentMeta] = ...
    generate_traveling_component_hat(K,KX,KY,bandWindow,cfg.wind);
componentT0 = real(ifft2(componentHat0));
targetRms = std(removedBand,0,'all')* ...
    sqrt(cfg.wind.replacementEnergyRatio);
rawRms = std(componentT0,0,'all');
assert(rawRms > eps,'The generated short-wave component has zero energy.');
componentScale = targetRms/rawRms;
componentHat0 = componentScale*componentHat0;
componentT0 = componentScale*componentT0;

omega = sqrt(cfg.wind.gravity*K);
psi = deg2rad(cfg.wind.propagationDirectionDeg);
signedDirection = sign(KX*cos(psi)+KY*sin(psi));
phaseT1 = exp(-1i*omega.*signedDirection*cfg.wind.previewTimeStep);
componentT1 = real(ifft2(componentHat0.*phaseT1));

Z = backgroundResidual+componentT0;
ZPreview = backgroundResidual+componentT1;
dZdt = real(ifft2(-1i*omega.*signedDirection.*componentHat0));

detection = detect_steep_directional_crest( ...
    background.X,background.Y,background.XLinear,background.YLinear, ...
    Z,componentT0,dx,dy,cfg);
metrics = calculate_metrics(background.Z,backgroundResidual,removedBand, ...
    componentT0,Z,ZPreview,KX,KY,componentHat0,dx,dy,cfg, ...
    componentMeta,detection);

surfaceData = struct();
surfaceData.X = background.X;
surfaceData.Y = background.Y;
surfaceData.Z = Z;
surfaceData.ZPreview = ZPreview;
surfaceData.XLinear = background.XLinear;
surfaceData.YLinear = background.YLinear;
surfaceData.ZOriginalBackground = background.Z;
surfaceData.ZBackgroundResidual = backgroundResidual;
surfaceData.ZRemovedBand = removedBand;
surfaceData.ZWindComponent = componentT0;
surfaceData.dZWindDt = dZdt;
surfaceData.componentHat0 = componentHat0;
surfaceData.spectrumShape = spectrumShape;
surfaceData.KX = KX;
surfaceData.KY = KY;
surfaceData.K = K;
surfaceData.bandWindow = bandWindow;
surfaceData.detection = detection;
surfaceData.metrics = metrics;
surfaceData.cfg = cfg;
end

function background = load_background(fileName)
assert(isfile(fileName),'Nonlinear background MAT file was not found: %s',fileName);
data = load(fileName,'X','Y','Z','XLinear','YLinear','ZLinear','cfg');
required = {'X','Y','Z','XLinear','YLinear'};
for i = 1:numel(required)
    assert(isfield(data,required{i}), ...
        'Background MAT file is missing variable %s.',required{i});
end
background = data;
assert(isequal(size(background.X),size(background.Y),size(background.Z)), ...
    'Background X, Y and Z arrays must have identical sizes.');
end

function spacing = infer_spacing(grid,dimension)
if dimension == 2
    differences = diff(grid(1,:),1,2);
else
    differences = diff(grid(:,1),1,1);
end
spacing = median(abs(differences),'all');
assert(isfinite(spacing) && spacing > 0,'Unable to infer grid spacing.');
end

function response = smooth_radial_bandpass(K,kLow,kHigh,transitionFraction)
lowStop = max(0,(1-transitionFraction)*kLow);
lowPass = kLow;
highPass = kHigh;
highStop = (1+transitionFraction)*kHigh;
response = zeros(size(K));
response(K >= lowPass & K <= highPass) = 1;
lower = K > lowStop & K < lowPass;
qLower = (K(lower)-lowStop)/max(lowPass-lowStop,eps);
response(lower) = 0.5*(1-cos(pi*qLower));
upper = K > highPass & K < highStop;
qUpper = (K(upper)-highPass)/max(highStop-highPass,eps);
response(upper) = 0.5*(1+cos(pi*qUpper));
end

function [componentHat,shape,meta] = generate_traveling_component_hat( ...
    K,KX,KY,bandWindow,wind)
[Ny,Nx] = size(K);
Ksafe = max(K,eps);
kLow = 2*pi/wind.maximumWavelength;
kHigh = 2*pi/wind.minimumWavelength;
kCenter = sqrt(kLow*kHigh);
sigmaLogK = log(kHigh/kLow)/4;
radialShape = exp(-0.5*(log(Ksafe/kCenter)/sigmaLogK).^2);

psi = deg2rad(wind.propagationDirectionDeg);
directionCosine = (KX*cos(psi)+KY*sin(psi))./Ksafe;
axisAngle = acos(min(max(abs(directionCosine),0),1));
sigmaAngle = deg2rad(wind.angularSpreadStdDeg);
directionShape = exp(-0.5*(axisAngle/max(sigmaAngle,eps)).^2);
shape = bandWindow.*radialShape.*directionShape;
shape(K == 0) = 0;

white = (randn(Ny,Nx)+1i*randn(Ny,Nx))/sqrt(2);
negX = [1,Nx:-1:2];
negY = [1,Ny:-1:2];
white = (white+conj(white(negY,negX)))/sqrt(2);
componentHat = sqrt(shape).*white;

% Enforce exact Hermitian symmetry after the directional shaping.
componentHat = (componentHat+conj(componentHat(negY,negX)))/2;
componentHat(1,1) = 0;
meta.centralWavenumber = kCenter;
meta.centralWavelength = 2*pi/kCenter;
end

function detection = detect_steep_directional_crest( ...
    X,Y,XParameter,YParameter,Z,component,dx,dy,cfg)
psi = deg2rad(cfg.wind.propagationDirectionDeg);
[Zx,Zy] = gradient(Z,dx,dy);
alongSlope = cos(psi)*Zx+sin(psi)*Zy;
[slopeX,slopeY] = gradient(alongSlope,dx,dy);
alongCurvature = cos(psi)*slopeX+sin(psi)*slopeY;

heightThreshold = array_quantile(Z,cfg.detection.minimumHeightQuantile);
maximumCrestSlope = cfg.detection.maximumAlongSlopeFraction* ...
    std(alongSlope,0,'all');
valid = Z >= heightThreshold & alongCurvature < 0 & ...
    abs(alongSlope) <= maximumCrestSlope;
valid = valid & X >= min(X,[],'all')+cfg.detection.edgeMargin & ...
    X <= max(X,[],'all')-cfg.detection.edgeMargin & ...
    Y >= min(Y,[],'all')+cfg.detection.edgeMargin & ...
    Y <= max(Y,[],'all')-cfg.detection.edgeMargin;
assert(any(valid,'all'),'No steep directional crest satisfies the detection gate.');

heightScore = normalize_positive(Z-heightThreshold);
curvatureScore = normalize_positive(-alongCurvature);
componentScore = normalize_positive(component);
score = cfg.detection.heightWeight*heightScore + ...
    cfg.detection.curvatureWeight*curvatureScore + ...
    cfg.detection.componentWeight*componentScore;
score(~valid) = -Inf;
[~,index] = max(score,[],'all');

detection.linearIndex = index;
detection.x = X(index);
detection.y = Y(index);
detection.parameterX = XParameter(index);
detection.parameterY = YParameter(index);
detection.z = Z(index);
detection.windComponentElevation = component(index);
detection.alongSlope = alongSlope(index);
detection.alongCurvature = alongCurvature(index);
detection.heightThreshold = heightThreshold;
detection.candidateCount = nnz(valid);
detection.score = score(index);
detection.propagationDirectionDeg = ...
    cfg.wind.propagationDirectionDeg;
detection.localSteepnessProxy = ...
    (2*pi/sqrt(cfg.wind.minimumWavelength* ...
    cfg.wind.maximumWavelength))*max(component(index),0);
detection.alongSlopeMap = alongSlope;
detection.alongCurvatureMap = alongCurvature;
end

function metrics = calculate_metrics(original,residual,removed,component, ...
    combined,preview,KX,KY,componentHat,dx,dy,cfg,componentMeta,detection)
psi = deg2rad(cfg.wind.propagationDirectionDeg);
energy = abs(componentHat).^2;
positiveTravel = (KX*cos(psi)+KY*sin(psi)) > 0;
positiveEnergy = energy.*positiveTravel;
directionX = sum(positiveEnergy.*KX./max(hypot(KX,KY),eps),'all');
directionY = sum(positiveEnergy.*KY./max(hypot(KX,KY),eps),'all');
estimatedDirection = mod(rad2deg(atan2(directionY,directionX)),360);

[originalX,originalY] = gradient(original,dx,dy);
[combinedX,combinedY] = gradient(combined,dx,dy);
metrics.originalRms = std(original,0,'all');
metrics.removedBandRms = std(removed,0,'all');
metrics.windComponentRms = std(component,0,'all');
metrics.combinedRms = std(combined,0,'all');
metrics.relativeRmsChange = (metrics.combinedRms-metrics.originalRms)/ ...
    max(metrics.originalRms,eps);
metrics.originalMss = mean(originalX.^2+originalY.^2,'all');
metrics.combinedMss = mean(combinedX.^2+combinedY.^2,'all');
metrics.residualReconstructionError = max(abs( ...
    residual+removed-original),[],'all');
metrics.previewChangeRms = std(preview-combined,0,'all');
metrics.estimatedPropagationDirectionDeg = estimatedDirection;
metrics.directionErrorDeg = angular_difference_deg( ...
    estimatedDirection,cfg.wind.propagationDirectionDeg);
metrics.centralWavelength = componentMeta.centralWavelength;
metrics.selectedCrestSteepnessProxy = detection.localSteepnessProxy;
end

function value = array_quantile(array,q)
values = sort(array(:));
index = 1+(numel(values)-1)*q;
lower = floor(index);
upper = ceil(index);
fraction = index-lower;
value = (1-fraction)*values(lower)+fraction*values(upper);
end

function normalized = normalize_positive(array)
array = max(array,0);
scale = array_quantile(array,0.99);
normalized = min(array/max(scale,eps),1);
end

function difference = angular_difference_deg(a,b)
difference = abs(mod(a-b+180,360)-180);
end

function validate_config(cfg)
assert(cfg.wind.minimumWavelength > 0 && ...
    cfg.wind.maximumWavelength > cfg.wind.minimumWavelength, ...
    'The wind-wave wavelength band is invalid.');
assert(cfg.wind.angularSpreadStdDeg > 0 && ...
    cfg.wind.angularSpreadStdDeg < 90, ...
    'angularSpreadStdDeg must lie in (0,90).');
assert(cfg.wind.replacementEnergyRatio > 0, ...
    'replacementEnergyRatio must be positive.');
assert(cfg.wind.radialTransitionFraction > 0 && ...
    cfg.wind.radialTransitionFraction < 0.5, ...
    'radialTransitionFraction must lie in (0,0.5).');
end
