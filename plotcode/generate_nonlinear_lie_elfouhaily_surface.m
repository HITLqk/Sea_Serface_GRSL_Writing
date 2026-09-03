function surfaceData = generate_nonlinear_lie_elfouhaily_surface(cfg)
%GENERATE_NONLINEAR_LIE_ELFOUHAILY_SURFACE Build a nonlinear random sea.
%   A linear directional Elfouhaily realization is corrected by a
%   broadband second-order bound-wave term. For a monochromatic wave, the
%   operator 0.5*|grad|(eta^2) recovers the deep-water Stokes coefficient
%   0.5*k*a^2 at the second harmonic. Wind speed changes only the
%   dimensionless nonlinear gain, and wind direction enters through a
%   directional Fourier multiplier.

arguments
    cfg (1,1) struct = default_nonlinear_lie_config()
end

validate_config(cfg);
rng(cfg.randomSeed,'twister');

Nx = round(cfg.domain.Lx/cfg.domain.dx);
Ny = round(cfg.domain.Ly/cfg.domain.dy);
assert(mod(Nx,2) == 0 && mod(Ny,2) == 0, ...
    'The x and y sample counts must be even.');

dx = cfg.domain.Lx/Nx;
dy = cfg.domain.Ly/Ny;
x = (0:Nx-1)*dx;
y = (0:Ny-1)*dy;
[X0,Y0] = meshgrid(x,y);

dkx = 2*pi/cfg.domain.Lx;
dky = 2*pi/cfg.domain.Ly;
kx = ifftshift((-Nx/2:Nx/2-1)*dkx);
ky = ifftshift((-Ny/2:Ny/2-1)*dky);
[KX,KY] = meshgrid(kx,ky);
K = hypot(KX,KY);
Ksafe = max(K,eps);

[Psi,spectrumMeta] = elfouhaily_directional_spectrum( ...
    K,KX,KY,cfg.sea);
spectralBandMask = K >= cfg.sea.minimumWavenumber & ...
    K <= cfg.sea.maximumWavenumber;
Psi(~spectralBandMask) = 0;
etaLinear = synthesize_real_surface(Psi,dkx,dky,Nx,Ny);
etaLinear = etaLinear-mean(etaLinear,'all');

nyquistRadius = min(pi/dx,pi/dy);
assert(cfg.lie.nonlinearOutputCutoff <= 0.98*nyquistRadius, ...
    ['The fixed nonlinear output cutoff %.3f rad/m is not resolved by ', ...
    'the current grid (Nyquist %.3f rad/m). Reduce dx/dy.'], ...
    cfg.lie.nonlinearOutputCutoff,nyquistRadius);
inputPass = (1-cfg.lie.filterTransitionFraction)* ...
    cfg.lie.nonlinearInputCutoff;
outputPass = (1-cfg.lie.filterTransitionFraction)* ...
    cfg.lie.nonlinearOutputCutoff;
inputFilter = smooth_radial_lowpass( ...
    K,inputPass,cfg.lie.nonlinearInputCutoff);
outputFilter = smooth_radial_lowpass( ...
    K,outputPass,cfg.lie.nonlinearOutputCutoff);

psiWind = deg2rad(cfg.sea.windDirectionDeg);
windProjection = abs((KX*cos(psiWind)+KY*sin(psiWind))./Ksafe);
directionWeight = (1-cfg.lie.directionStrength) + ...
    cfg.lie.directionStrength*windProjection.^ ...
    cfg.lie.directionExponent;
directionWeight(K == 0) = 0;

windGain = cfg.lie.baseGain* ...
    (cfg.sea.U10/cfg.lie.referenceWindSpeed)^cfg.lie.windExponent;
windGain = min(max(windGain,cfg.lie.minimumWindGain), ...
    cfg.lie.maximumWindGain);

etaLinearHat = fft2(etaLinear);
etaForNonlinear = real(ifft2(etaLinearHat.*inputFilter));
etaSquaredHat = dealiased_square_hat( ...
    etaForNonlinear,cfg.lie.dealiasExpansion);
etaSquaredHat(1,1) = 0;
etaSecondHat = 0.5*K.*directionWeight.*outputFilter.*etaSquaredHat;
etaSecond = real(ifft2(etaSecondHat));
etaSecond = etaSecond-mean(etaSecond,'all');
etaNonlinear = etaLinear+windGain*etaSecond;
etaNonlinear = etaNonlinear-mean(etaNonlinear,'all');

% Riesz displacement gives the Lie variables a horizontal component. A
% line search prevents local folding of the x-y parameterization.
Dx = real(ifft2(-1i*(KX./Ksafe).*etaLinearHat.*inputFilter));
Dy = real(ifft2(-1i*(KY./Ksafe).*etaLinearHat.*inputFilter));

lambda = cfg.lie.horizontalDisplacement*windGain;
[lambda,minJacobian,backoffSteps] = stable_horizontal_gain( ...
    Dx,Dy,dx,dy,lambda,cfg.lie);
X = X0+lambda*Dx;
Y = Y0+lambda*Dy;
Z = etaNonlinear;

metrics = calculate_metrics(etaLinear,etaSecond,Z, ...
    minJacobian,lambda,windGain,dx,dy);

surfaceData = struct();
surfaceData.X0 = X0;
surfaceData.Y0 = Y0;
surfaceData.ZLinear = etaLinear;
surfaceData.X = X;
surfaceData.Y = Y;
surfaceData.Z = Z;
surfaceData.secondOrderElevation = etaSecond;
surfaceData.horizontalDisplacementX = lambda*Dx;
surfaceData.horizontalDisplacementY = lambda*Dy;
surfaceData.Psi = Psi;
surfaceData.KX = KX;
surfaceData.KY = KY;
surfaceData.K = K;
surfaceData.spectrumMeta = spectrumMeta;
surfaceData.metrics = metrics;
surfaceData.metrics.horizontalBackoffSteps = backoffSteps;
surfaceData.metrics.nonlinearInputCutoff = ...
    cfg.lie.nonlinearInputCutoff;
surfaceData.metrics.nonlinearOutputCutoff = ...
    cfg.lie.nonlinearOutputCutoff;
surfaceData.metrics.minimumSynthesisWavenumber = ...
    cfg.sea.minimumWavenumber;
surfaceData.metrics.maximumSynthesisWavenumber = ...
    cfg.sea.maximumWavenumber;
surfaceData.cfg = cfg;
end

function response = smooth_radial_lowpass(K,passWavenumber,stopWavenumber)
assert(passWavenumber >= 0 && stopWavenumber > passWavenumber, ...
    'Low-pass stop wavenumber must exceed its pass wavenumber.');
response = ones(size(K));
response(K >= stopWavenumber) = 0;
transition = K > passWavenumber & K < stopWavenumber;
q = (K(transition)-passWavenumber)/ ...
    (stopWavenumber-passWavenumber);
response(transition) = 0.5*(1+cos(pi*q));
end

function productHat = dealiased_square_hat(field,expansion)
% Compute field.^2 on an expanded Fourier grid, then crop spectrally.
[Ny,Nx] = size(field);
paddedNy = 2*ceil(expansion*Ny/2);
paddedNx = 2*ceil(expansion*Nx/2);

fieldHat = fftshift(fft2(field));
paddedHat = complex(zeros(paddedNy,paddedNx));
rowStart = (paddedNy-Ny)/2+1;
columnStart = (paddedNx-Nx)/2+1;
rows = rowStart:rowStart+Ny-1;
columns = columnStart:columnStart+Nx-1;
paddedHat(rows,columns) = fieldHat;

paddingScale = (paddedNy*paddedNx)/(Ny*Nx);
paddedField = real(ifft2(ifftshift(paddedHat)))*paddingScale;
paddedProductHat = fftshift(fft2(paddedField.^2));
cropped = paddedProductHat(rows,columns)/paddingScale;
productHat = ifftshift(cropped);
end

function [Psi,meta] = elfouhaily_directional_spectrum(K,KX,KY,sea)
% Unified Elfouhaily spectrum retained from the original project, with
% explicit zero-wavenumber handling and a normalized directional factor.
g = 9.81;
Cd10N = 0.00144;
km = 370.0;
cm = 0.23;
age = sea.inverseWaveAge;
U10 = sea.U10;
uStar = sqrt(Cd10N)*U10;
k0 = g/U10^2;
kp = k0*age^2;
cp = sqrt(g/kp);
sigma = 0.08*(1+4*age^(-3));
alphaP = 0.006*age^0.55;

if uStar <= cm
    alphaM = 0.01*(1+log(uStar/cm));
else
    alphaM = 0.01*(1+3*log(uStar/cm));
end
if age <= 1
    gamma = 1.7;
else
    gamma = 1.7+6*log(age);
end

Ksafe = max(K,eps);
c = sqrt((g./Ksafe).*(1+(Ksafe/km).^2));
Lpm = exp(-1.25*(kp./Ksafe).^2);
Gamma = exp(-0.5*((sqrt(Ksafe/kp)-1)/sigma).^2);
Jp = gamma.^Gamma;
Fp = Lpm.*Jp.*exp(-age/sqrt(10).*(sqrt(Ksafe/kp)-1));
Fm = Lpm.*Jp.*exp(-0.25*(Ksafe/km-1).^2);
Bl = 0.5*alphaP*(cp./c).*Fp;
Bh = 0.5*alphaM*(cm./c).*Fm;
S = (Bl+Bh)./Ksafe.^3;

a0 = log(2)/4;
ap = 4;
am = 0.13*uStar/cm;
delta = tanh(a0+ap*(c/cp).^2.5+am*(cm./c).^2.5);
phi = atan2(KY,KX);
phiWind = deg2rad(sea.windDirectionDeg);
spread = (1+delta.*cos(2*(phi-phiWind)))/(2*pi);
Psi = max(S.*spread./Ksafe,0);
Psi(K == 0) = 0;

meta.peakWavenumber = kp;
meta.peakWavelength = 2*pi/kp;
meta.frictionVelocity = uStar;
meta.alphaLong = alphaP;
meta.alphaShort = alphaM;
end

function eta = synthesize_real_surface(Psi,dkx,dky,Nx,Ny)
white = (randn(Ny,Nx)+1i*randn(Ny,Nx))/sqrt(2);
negX = [1,Nx:-1:2];
negY = [1,Ny:-1:2];
white = (white+conj(white(negY,negX)))/sqrt(2);
coefficients = sqrt(Psi*dkx*dky).*white;
etaComplex = ifft2(coefficients)*Nx*Ny;
imaginaryResidual = max(abs(imag(etaComplex)),[],'all');
assert(imaginaryResidual < 1e-10, ...
    'Hermitian spectral sampling did not produce a real surface.');
eta = real(etaComplex);
end

function [lambda,minJ,steps] = stable_horizontal_gain(Dx,Dy,dx,dy,lambda,lie)
[Dxx,Dxy] = gradient(Dx,dx,dy);
[Dyx,Dyy] = gradient(Dy,dx,dy);
for steps = 0:lie.maximumBackoffSteps
    J = (1+lambda*Dxx).*(1+lambda*Dyy) - ...
        lambda^2.*Dxy.*Dyx;
    minJ = min(J,[],'all');
    if minJ >= lie.minimumJacobian
        return;
    end
    lambda = lambda*lie.backoffFactor;
end
error('Unable to obtain a non-folding horizontal Lie mapping.');
end

function metrics = calculate_metrics(zLinear,zSecond,zNonlinear, ...
    minJacobian,lambda,windGain,dx,dy)
[sxLinear,syLinear] = gradient(zLinear,dx,dy);
[sxNonlinear,syNonlinear] = gradient(zNonlinear,dx,dy);
metrics.linearHs = 4*std(zLinear,0,'all');
metrics.nonlinearHs = 4*std(zNonlinear,0,'all');
metrics.linearSkewness = standardized_moment(zLinear,3);
metrics.nonlinearSkewness = standardized_moment(zNonlinear,3);
metrics.linearKurtosis = standardized_moment(zLinear,4);
metrics.nonlinearKurtosis = standardized_moment(zNonlinear,4);
metrics.linearMss = mean(sxLinear.^2+syLinear.^2,'all');
metrics.nonlinearMss = mean(sxNonlinear.^2+syNonlinear.^2,'all');
metrics.maximumCrest = max(zNonlinear,[],'all');
metrics.minimumTrough = min(zNonlinear,[],'all');
metrics.crestTroughRatio = metrics.maximumCrest/ ...
    max(abs(metrics.minimumTrough),eps);
metrics.secondOrderRms = sqrt(mean(zSecond.^2,'all'));
metrics.minimumHorizontalJacobian = minJacobian;
metrics.appliedHorizontalGain = lambda;
metrics.appliedWindGain = windGain;
end

function value = standardized_moment(z,order)
centered = z-mean(z,'all');
sigma = std(centered,0,'all');
value = mean((centered/max(sigma,eps)).^order,'all');
end

function validate_config(cfg)
assert(cfg.sea.U10 > 0,'U10 must be positive.');
assert(cfg.sea.inverseWaveAge > 0,'Inverse wave age must be positive.');
assert(cfg.sea.minimumWavenumber >= 0 && ...
    cfg.sea.maximumWavenumber > cfg.sea.minimumWavenumber, ...
    'The sea-surface synthesis wavenumber band is invalid.');
assert(cfg.domain.Lx > 0 && cfg.domain.Ly > 0, ...
    'Domain lengths must be positive.');
assert(cfg.lie.directionStrength >= 0 && ...
    cfg.lie.directionStrength <= 1, ...
    'directionStrength must lie in [0,1].');
assert(cfg.lie.nonlinearInputCutoff > 0 && ...
    cfg.lie.nonlinearOutputCutoff >= cfg.lie.nonlinearInputCutoff, ...
    'The fixed nonlinear wavenumber cutoffs are invalid.');
assert(cfg.lie.filterTransitionFraction > 0 && ...
    cfg.lie.filterTransitionFraction < 1, ...
    'filterTransitionFraction must lie in (0,1).');
assert(cfg.lie.dealiasExpansion >= 1.5, ...
    'Quadratic de-aliasing requires an expansion of at least 1.5.');
end
