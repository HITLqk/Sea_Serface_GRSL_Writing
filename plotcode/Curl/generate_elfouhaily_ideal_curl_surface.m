function surfaceData = generate_elfouhaily_ideal_curl_surface(cfg)
%GENERATE_ELFOUHAILY_IDEAL_CURL_SURFACE Add a local parametric curl.
%   The event is placed on a breaking-eligible 2-D Elfouhaily crest using
%   elevation, propagation-direction slope, and curvature. The original
%   structured mesh is deformed in both horizontal and vertical position.

arguments
    cfg (1,1) struct = default_elfouhaily_ideal_curl_config()
end

rng(cfg.randomSeed, 'twister');

Nx = round(cfg.domain.Lx/cfg.domain.dx);
Ny = round(cfg.domain.Ly/cfg.domain.dy);
dx = cfg.domain.Lx/Nx;
dy = cfg.domain.Ly/Ny;
x = (0:Nx-1)*dx;
y = (0:Ny-1)*dy;
[X0,Y0] = meshgrid(x,y);

Z0 = synthesize_elfouhaily_surface(Nx,Ny,cfg.domain.Lx, ...
    cfg.domain.Ly,cfg.sea);
detection = detect_breaking_eligible_crest(X0,Y0,Z0,cfg,dx,dy);

psi = deg2rad(cfg.detection.propagationDirectionDeg);
u = cos(psi).*(X0-detection.x) + sin(psi).*(Y0-detection.y);
v = -sin(psi).*(X0-detection.x) + cos(psi).*(Y0-detection.y);

ridge = track_crest_ridge(u,v,Z0,detection,cfg,dx,dy);
vQuery = min(max(v,ridge.v(1)),ridge.v(end));
uCenter = interp1(ridge.v,ridge.u,vQuery,'pchip');
zCrest = interp1(ridge.v,ridge.z,vQuery,'pchip');
ur = u-uCenter;

localWaveHeight = estimate_local_wave_height(u,v,Z0,detection,cfg,dx,dy);

a = cfg.curl.crestHalfLength;
b = cfg.curl.coreHalfWidth;
Wt = compact_cosine(ur/(b+cfg.curl.transitionWidth)).* ...
    compact_cosine(v/a);
Wc = compact_cosine(ur/b).*compact_cosine(v/a);

sigmaRear = cfg.curl.rearSigmaFraction*b;
sigmaFront = cfg.curl.frontSigmaFraction*b;
sigma = sigmaFront*ones(size(ur));
sigma(ur < 0) = sigmaRear;
sigmaEvolution = cfg.curl.evolutionSigmaFraction* ...
    (b+cfg.curl.transitionWidth);

crestShape = exp(-0.5*(ur./sigma).^2);
evolutionShape = exp(-0.5*(ur/sigmaEvolution).^2);
crestLift = cfg.curl.crestLiftFraction*localWaveHeight;
evolutionLift = cfg.curl.evolutionLiftFraction*localWaveHeight;
forwardLean = cfg.curl.forwardLeanFraction*localWaveHeight;

zPre = Z0 + (crestLift*crestShape + ...
    evolutionLift*evolutionShape).*Wt;
crestwiseWeight = compact_cosine(v/a);
uPre = u + forwardLean*evolutionShape.*Wt;

% A plunging breaker is multi-valued in z(u), so its core must be built as
% a material-parametric curve. The upper branch advances, a half-ellipse
% turns through the nose, and the lower branch continues forward at depth.
[uTemplate,zTemplate,profileMask] = plunging_profile(ur,uCenter, ...
    zCrest,zPre,localWaveHeight,cfg);
profileWeight = profileMask.*crestwiseWeight;
uFinal = uPre + profileWeight.*(uTemplate-uPre);
Z = zPre + profileWeight.*(zTemplate-zPre);
X = detection.x + cos(psi).*uFinal - sin(psi).*v;
Y = detection.y + sin(psi).*uFinal + cos(psi).*v;

[duFinalDx,duFinalDy] = gradient(uFinal,dx,dy);
propagationJacobian = cos(psi).*duFinalDx + ...
    sin(psi).*duFinalDy;
coreMask = Wc >= cfg.curl.coreMaskThreshold;
overturningMask = coreMask & propagationJacobian < 0;
transitionMask = Wt > 0 & ~coreMask;

faces = structured_triangles(Nx,Ny);
curlFacetMask = any(reshape(coreMask(faces),size(faces)),2);

surfaceData = struct();
surfaceData.X0 = X0;
surfaceData.Y0 = Y0;
surfaceData.Z0 = Z0;
surfaceData.X = X;
surfaceData.Y = Y;
surfaceData.Z = Z;
surfaceData.localU = u;
surfaceData.localV = v;
surfaceData.localURidge = uCenter;
surfaceData.localURelative = ur;
surfaceData.localUFinal = uFinal;
surfaceData.zPre = zPre;
surfaceData.thetaCurl = zeros(size(Z));
surfaceData.coreWeight = Wc;
surfaceData.transitionWeight = Wt;
surfaceData.curlMask = coreMask;
surfaceData.transitionMask = transitionMask;
surfaceData.propagationJacobian = propagationJacobian;
surfaceData.overturningMask = overturningMask;
surfaceData.faces = faces;
surfaceData.curlFacetMask = curlFacetMask;
surfaceData.verticesBaseline = [X0(:),Y0(:),Z0(:)];
surfaceData.vertices = [X(:),Y(:),Z(:)];
surfaceData.detection = detection;
surfaceData.ridge = ridge;
surfaceData.localWaveHeight = localWaveHeight;
surfaceData.cfg = cfg;

horizontalDisplacement = hypot(X-X0,Y-Y0);
forwardDisplacement = uFinal-u;
downwardDisplacement = Z0-Z;
surfaceData.metrics.maxElevationChange = max(abs(Z-Z0),[],'all');
surfaceData.metrics.maxUpwardDisplacement = max(Z-Z0,[],'all');
surfaceData.metrics.maxHorizontalDisplacement = ...
    max(horizontalDisplacement,[],'all');
surfaceData.metrics.maxForwardDisplacement = ...
    max(forwardDisplacement(coreMask));
surfaceData.metrics.maxDownwardDisplacement = ...
    max(downwardDisplacement(coreMask));
surfaceData.metrics.downwardToLocalHeight = ...
    surfaceData.metrics.maxDownwardDisplacement/localWaveHeight;
surfaceData.metrics.curlPointCount = nnz(coreMask);
surfaceData.metrics.transitionPointCount = nnz(transitionMask);
surfaceData.metrics.minimumPropagationJacobian = ...
    min(propagationJacobian(coreMask));
surfaceData.metrics.overturningPointCount = nnz(overturningMask);
surfaceData.metrics.forwardOverturningFraction = nnz( ...
    overturningMask & ur > 0)/max(nnz(overturningMask),1);
downwardCore = downwardDisplacement;
downwardCore(~coreMask) = -Inf;
[~,lowestIndex] = max(downwardCore,[],'all','linear');
surfaceData.metrics.lowestLipRelativeU = ur(lowestIndex);
outside = Wt == 0;
surfaceData.metrics.maxOutsideDisplacement = max([ ...
    abs(X(outside)-X0(outside)); abs(Y(outside)-Y0(outside)); ...
    abs(Z(outside)-Z0(outside))]);
surfaceData.metrics.crestwiseRidgeStd = std(ridge.u,0,'all');
end

function [uTarget,zTarget,mask] = plunging_profile(q,uC,zC,zBackground,H,cfg)
qRear = cfg.curl.profileRear;
qTop = cfg.curl.noseStart;
qBottom = cfg.curl.noseEnd;
qEnd = cfg.curl.lowerEnd;

uTarget = uC+q;
zTarget = zBackground;
mask = zeros(size(q));

upper = q >= qRear & q < qTop;
tUpper = smoothstep((q(upper)-qRear)/(qTop-qRear));
advance = cfg.curl.upperAdvanceFraction*H;
zLipTop = zC-cfg.curl.lipTopDropFraction*H;
uTarget(upper) = uC(upper)+q(upper)+advance.*tUpper;
zTarget(upper) = zTarget(upper)+(zLipTop(upper)-zTarget(upper)).* ...
    cfg.curl.plateauFraction.*tUpper;
mask(upper) = smoothstep((q(upper)-qRear)/(0.18*(qTop-qRear)));

uNoseBase = uC+qTop+advance;
zTop = zLipTop;
zBottom = zC-cfg.curl.lowerBranchDropFraction*H;
nose = q >= qTop & q <= qBottom;
tNose = (q(nose)-qTop)/(qBottom-qTop);
uTarget(nose) = uNoseBase(nose)+ ...
    cfg.curl.noseRadiusFraction*H.*sin(pi*tNose);
zTarget(nose) = zBottom(nose)+0.5.*(zTop(nose)-zBottom(nose)).* ...
    (1+cos(pi*tNose));
mask(nose) = 1;

lower = q > qBottom & q <= qEnd;
tLower = (q(lower)-qBottom)/(qEnd-qBottom);
uTarget(lower) = uNoseBase(lower)+ ...
    cfg.curl.lowerAdvanceFactor*(q(lower)-qBottom);
zBackgroundJoin = zBackground(lower);
blendUp = smoothstep(tLower);
zTarget(lower) = (1-blendUp).*zBottom(lower)+ ...
    blendUp.*zBackgroundJoin;
mask(lower) = 1-smoothstep(max((tLower-0.72)/0.28,0));

end

function y = smoothstep(x)
x = min(max(x,0),1);
y = x.^2.*(3-2*x);
end

function detection = detect_breaking_eligible_crest(X,Y,Z,cfg,dx,dy)
sigmaCells = cfg.detection.smoothingLength/max(0.5*(dx+dy),eps);
Zsmooth = gaussian_smooth(Z,sigmaCells);
[zx,zy] = gradient(Zsmooth,dx,dy);
psi = deg2rad(cfg.detection.propagationDirectionDeg);
etaU = cos(psi).*zx + sin(psi).*zy;
[etaUx,etaUy] = gradient(etaU,dx,dy);
etaUU = cos(psi).*etaUx + sin(psi).*etaUy;

heightThreshold = sample_quantile(Zsmooth,cfg.detection.heightQuantile);
slopeTolerance = cfg.detection.slopeToleranceFactor*std(etaU,0,'all');
negativeCurvature = max(-etaUU,0);
curvatureThreshold = sample_quantile(negativeCurvature( ...
    negativeCurvature > 0),cfg.detection.curvatureQuantile);

margin = cfg.detection.edgeMargin;
valid = X >= margin & X <= cfg.domain.Lx-margin & ...
    Y >= margin & Y <= cfg.domain.Ly-margin;
heightMask = valid & Zsmooth >= heightThreshold;
slopeMask = valid & abs(etaU) <= slopeTolerance;
curvatureMask = valid & etaUU <= -curvatureThreshold;
candidates = heightMask & slopeMask & curvatureMask;

if ~any(candidates,'all')
    curvatureMask = valid & etaUU < 0;
    candidates = heightMask & slopeMask & curvatureMask;
end
assert(any(candidates,'all'), ...
    'No crest satisfies the elevation, slope, and curvature conditions.');

x = X(1,:);
y = Y(:,1);
dForward = cfg.detection.forwardSlopeDistance;
forwardSlope = interp2(x,y,max(-etaU,0), ...
    X+dForward*cos(psi),Y+dForward*sin(psi),'linear',0);

heightScore = normalize_positive(Zsmooth-heightThreshold,valid);
curvatureScore = normalize_positive(negativeCurvature,valid);
forwardSlopeScore = normalize_positive(forwardSlope,valid);
score = cfg.detection.heightWeight*heightScore + ...
    cfg.detection.curvatureWeight*curvatureScore + ...
    cfg.detection.forwardSlopeWeight*forwardSlopeScore;
score(~candidates) = -Inf;
[~,coarseIndex] = max(score(:));

refine = valid & hypot(X-X(coarseIndex),Y-Y(coarseIndex)) <= ...
    cfg.detection.refineRadius & etaUU < 0;
indices = find(refine);
[~,localIndex] = max(Z(indices));
index = indices(localIndex);

detection.linearIndex = index;
detection.coarseIndex = coarseIndex;
detection.x = X(index);
detection.y = Y(index);
detection.z = Z(index);
detection.smoothedElevation = Zsmooth(index);
detection.heightThreshold = heightThreshold;
detection.slopeTolerance = slopeTolerance;
detection.curvatureThreshold = curvatureThreshold;
detection.candidateCount = nnz(candidates);
detection.Zsmooth = Zsmooth;
detection.directionalSlope = etaU;
detection.directionalCurvature = etaUU;
detection.forwardSlope = forwardSlope;
detection.heightMask = heightMask;
detection.slopeMask = slopeMask;
detection.curvatureMask = curvatureMask;
detection.candidateMask = candidates;
detection.score = score;
end

function ridge = track_crest_ridge(u,v,Z,detection,cfg,dx,dy)
a = cfg.curl.crestHalfLength;
dv = max(dx,dy);
nodeCount = max(31,2*ceil(a/dv)+1);
vNodes = linspace(-a,a,nodeCount);
uNodes = zeros(size(vNodes));
zNodes = zeros(size(vNodes));
for k = 1:nodeCount
    band = abs(v-vNodes(k)) <= 0.75*(vNodes(2)-vNodes(1)) & ...
        abs(u) <= cfg.curl.ridgeSearchHalfWidth;
    indices = find(band);
    if isempty(indices)
        uNodes(k) = 0;
        zNodes(k) = detection.z;
    else
        ridgeScore = detection.Zsmooth(indices) - ...
            0.06*abs(u(indices));
        [~,best] = max(ridgeScore);
        index = indices(best);
        uNodes(k) = u(index);
        zNodes(k) = Z(index);
    end
end

window = max(3,2*floor(cfg.curl.ridgeSmoothSamples/2)+1);
kernel = ones(1,window)/window;
uNodes = conv(uNodes,kernel,'same');
zNodes = conv(zNodes,kernel,'same');
edge = floor(window/2);
uNodes(1:edge) = uNodes(edge+1);
uNodes(end-edge+1:end) = uNodes(end-edge);
zNodes(1:edge) = zNodes(edge+1);
zNodes(end-edge+1:end) = zNodes(end-edge);
uNodes = uNodes-interp1(vNodes,uNodes,0,'linear');

ridge.v = vNodes;
ridge.u = uNodes;
ridge.z = zNodes;
end

function Hlocal = estimate_local_wave_height(u,v,Z,detection,cfg,dx,dy)
strip = abs(v) <= 1.5*max(dx,dy) & ...
    abs(u) <= cfg.curl.localHeightRadius;
uLine = u(strip);
zLine = Z(strip);
[uLine,order] = sort(uLine);
zLine = zLine(order);
guard = 0.25*cfg.curl.coreHalfWidth;
left = zLine(uLine < -guard);
right = zLine(uLine > guard);
if isempty(left) || isempty(right)
    Hlocal = 0.5*cfg.sea.targetHs;
else
    troughLevel = 0.5*(min(left)+min(right));
    Hlocal = detection.z-troughLevel;
end
Hlocal = min(max(Hlocal,0.35*cfg.sea.targetHs), ...
    1.25*cfg.sea.targetHs);
end

function values = normalize_positive(values,valid)
values = max(values,0);
scale = sample_quantile(values(valid),0.95);
values = min(values/max(scale,eps),1);
end

function value = sample_quantile(samples,p)
samples = sort(samples(isfinite(samples)));
assert(~isempty(samples),'Cannot calculate a quantile from empty data.');
position = 1+(numel(samples)-1)*min(max(p,0),1);
lowerIndex = floor(position);
upperIndex = ceil(position);
weight = position-lowerIndex;
value = (1-weight)*samples(lowerIndex)+weight*samples(upperIndex);
end

function Z = synthesize_elfouhaily_surface(Nx,Ny,Lx,Ly,sea)
g = 9.81;
km = 370.0;
cm = 0.23;
Omega = sea.inverseWaveAge;
uStar = sqrt(1e-3*(0.8+0.065*sea.U10))*sea.U10;
k0 = g/sea.U10^2;
kp = k0*Omega^2;
cp = sqrt(g/kp);

if Omega <= 1
    gamma = 1.7;
else
    gamma = 1.7+6*log10(Omega);
end
alphaP = 6e-3*sqrt(Omega);
alphaM = 0.01*(1+log(max(uStar/cm,eps)));

dkx = 2*pi/Lx;
dky = 2*pi/Ly;
kx = ifftshift((-floor(Nx/2):ceil(Nx/2)-1)*dkx);
ky = ifftshift((-floor(Ny/2):ceil(Ny/2)-1)*dky);
[KX,KY] = meshgrid(kx,ky);
K = hypot(KX,KY);
Ksafe = max(K,eps);

c = sqrt(g./Ksafe.*(1+(Ksafe/km).^2));
Lpm = exp(-1.25*(kp./Ksafe).^2);
Gamma = exp(-(sqrt(Ksafe/kp)-1).^2/(2*0.08^2));
Fp = Lpm.*gamma.^Gamma.*exp(-Omega/sqrt(10).* ...
    (sqrt(Ksafe/kp)-1));
Bl = 0.5*alphaP.*(cp./c).*Fp;
Fm = exp(-0.25*(Ksafe/km-1).^2);
Bh = 0.5*alphaM.*(cm./c).*Fm;
S = (Bl+Bh)./Ksafe.^3;

windDirection = deg2rad(sea.windDirectionDeg);
phi = atan2(KY,KX)-windDirection;
spread = 2/pi.*cos(phi).^2;
Psi = max(S.*spread./Ksafe,0);
Psi(K == 0) = 0;

white = (randn(Ny,Nx)+1i*randn(Ny,Nx))/sqrt(2);
negX = [1,Nx:-1:2];
negY = [1,Ny:-1:2];
white = (white+conj(white(negY,negX)))/sqrt(2);
H = sqrt(Psi*dkx*dky).*white;
Z = real(ifft2(H)*Nx*Ny);
Z = Z-mean(Z,'all');

if sea.targetHs > 0
    Z = Z*sea.targetHs/max(4*std(Z,0,'all'),eps);
end
end

function Zs = gaussian_smooth(Z,sigmaCells)
if sigmaCells <= 0
    Zs = Z;
    return;
end
radius = max(1,ceil(3*sigmaCells));
q = -radius:radius;
kernel = exp(-0.5*(q/sigmaCells).^2);
kernel = kernel/sum(kernel);
Zs = conv2(kernel,kernel,Z,'same');
end

function w = compact_cosine(q)
w = zeros(size(q));
inside = abs(q) < 1;
w(inside) = 0.5*(1+cos(pi*q(inside)));
end

function faces = structured_triangles(Nx,Ny)
[i,j] = meshgrid(1:Nx-1,1:Ny-1);
v00 = sub2ind([Ny,Nx],j(:),i(:));
v10 = sub2ind([Ny,Nx],j(:),i(:)+1);
v01 = sub2ind([Ny,Nx],j(:)+1,i(:));
v11 = sub2ind([Ny,Nx],j(:)+1,i(:)+1);
faces = [v00,v10,v11; v00,v11,v01];
end
