function result = extract_breaker_morphology_metrics(surfaceData,cfg)
%EXTRACT_BREAKER_MORPHOLOGY_METRICS Measure the center-plane curl geometry.

dx = surfaceData.cfg.domain.Lx/size(surfaceData.X,2);
dy = surfaceData.cfg.domain.Ly/size(surfaceData.Y,1);
stripHalfWidth = cfg.extraction.centerlineHalfWidthCells*max(dx,dy);
centerMask = abs(surfaceData.localV) <= stripHalfWidth;

u0 = surfaceData.localU(centerMask);
u = surfaceData.localUFinal(centerMask);
z0 = surfaceData.Z0(centerMask);
z = surfaceData.Z(centerMask);
thetaCurl = surfaceData.thetaCurl(centerMask);

[u0,order] = sort(u0);
u = u(order);
z0 = z0(order);
z = z(order);
thetaCurl = thetaCurl(order);

% Average duplicate propagation coordinates if a rotated center strip is used.
[u0,~,group] = unique(round(u0/max(dx,eps))*max(dx,eps),'stable');
u = accumarray(group,u,[],@mean);
z0 = accumarray(group,z0,[],@mean);
z = accumarray(group,z,[],@mean);
thetaCurl = accumarray(group,thetaCurl,[],@mean);

active = thetaCurl >= cfg.extraction.activeCurlFraction*max(thetaCurl);
assert(nnz(active) >= 5,'Too few center-plane points belong to the curl.');

uActive = u(active);
zActive = z(active);
rx = max(uActive)-min(uActive);
ry = max(zActive)-min(zActive);

uLimit = max(uActive)-cfg.extraction.frontForwardFraction*rx;
zLimits = min(zActive)+cfg.extraction.frontVerticalFraction*ry;
front = active & u >= uLimit & z >= zLimits(1) & z <= zLimits(2);

% Retain the steep descending face when the forward window includes both
% sides of the folded parametric curve.
du = gradient(u,u0);
dz = gradient(z,u0);
descending = dz.*du < 0;
if nnz(front & descending) >= 3
    front = front & descending;
end
assert(nnz(front) >= 3,'Front-face extraction produced fewer than 3 points.');

frontU = u(front);
frontZ = z(front);
centered = [frontU-mean(frontU),frontZ-mean(frontZ)];
[~,~,vectors] = svd(centered,0);
direction = vectors(:,1);
frontFaceAngleDeg = atan2d(abs(direction(2)),abs(direction(1)));

lambdaPeak = peak_wavelength(surfaceData.cfg.sea);
switch lower(cfg.normalization.mode)
    case 'spectralpeak'
        wavelength = lambdaPeak;
    case 'fixed'
        wavelength = cfg.normalization.fixedWavelength;
    otherwise
        error('Unknown wavelength mode: %s',cfg.normalization.mode);
end
assert(isfinite(wavelength) && wavelength > 0, ...
    'The normalization wavelength must be positive.');

result.frontFaceAngleDeg = frontFaceAngleDeg;
result.rx = rx;
result.ry = ry;
result.wavelength = wavelength;
result.spectralPeakWavelength = lambdaPeak;
result.rxOverLambda = rx/wavelength;
result.ryOverLambda = ry/wavelength;
result.anglePass = in_range(frontFaceAngleDeg, ...
    cfg.reference.frontFaceAngleDeg);
result.rxPass = in_range(result.rxOverLambda, ...
    cfg.reference.rxOverLambda);
result.ryPass = in_range(result.ryOverLambda, ...
    cfg.reference.ryOverLambda);
result.allPass = result.anglePass && result.rxPass && result.ryPass;
result.profile.u0 = u0;
result.profile.u = u;
result.profile.z0 = z0;
result.profile.z = z;
result.profile.active = active;
result.profile.front = front;
result.frontFit.mean = [mean(frontU),mean(frontZ)];
result.frontFit.direction = direction;
end

function lambda = peak_wavelength(sea)
g = 9.81;
kp = g/sea.U10^2*sea.inverseWaveAge^2;
lambda = 2*pi/kp;
end

function tf = in_range(value,bounds)
tf = value >= bounds(1) && value <= bounds(2);
end
