function realization = synthesize_two_group_realization(U10,seed,cfg)
%SYNTHESIZE_TWO_GROUP_REALIZATION Paired linear and modified-Lie surfaces.

rng(seed,'twister');
g = 9.81;
kp = g*(cfg.inverseWaveAge/U10)^2;
dk = kp/cfg.primaryPeakSamples;
N = cfg.primaryGridSize;
kAxis = [0:N/2-1,-N/2:-1]*dk;
[KX,KY] = meshgrid(kAxis,kAxis);
K = hypot(KX,KY);

[Psi,meta] = thesis_elfouhaily_spectrum(K,KX,KY,U10, ...
    cfg.inverseWaveAge,cfg.windDirectionDeg);
primaryMask = K <= cfg.primaryMaximumPeakMultiple*kp;
Hlinear = sample_hermitian_coefficients(Psi.*primaryMask,dk,dk);

Hbreaking = apply_thesis_modified_lie(Hlinear,KX,KY,K,U10,kp,cfg);
linearMss = spectral_mss(Hlinear,KX,KY);
breakingMssRaw = spectral_mss(Hbreaking,KX,KY);

bandLow = cfg.primaryMaximumPeakMultiple*kp;
shortMss = struct('along',0,'cross',0,'total',0);
while bandLow < cfg.maximumOpticalWavenumber
    bandHigh = min(2*bandLow,cfg.maximumOpticalWavenumber);
    tileDk = bandLow/cfg.shortWaveModesBelowBand;
    tileN = cfg.shortWaveTileSize;
    tileAxis = [0:tileN/2-1,-tileN/2:-1]*tileDk;
    [tileKX,tileKY] = meshgrid(tileAxis,tileAxis);
    tileK = hypot(tileKX,tileKY);
    tilePsi = thesis_elfouhaily_spectrum(tileK,tileKX,tileKY,U10, ...
        cfg.inverseWaveAge,cfg.windDirectionDeg);
    bandMask = tileK >= bandLow & tileK < bandHigh;
    Hband = sample_hermitian_coefficients(tilePsi.*bandMask,tileDk,tileDk);
    bandMss = spectral_mss(Hband,tileKX,tileKY);
    shortMss.along = shortMss.along+bandMss.along;
    shortMss.cross = shortMss.cross+bandMss.cross;
    shortMss.total = shortMss.total+bandMss.total;
    bandLow = bandHigh;
end

linearTotal = add_mss(linearMss,shortMss);
breakingRawTotal = add_mss(breakingMssRaw,shortMss);
if cfg.enableEmpiricalMssDressing
    targetMode = 'fusion';
    if isfield(cfg,'empiricalTargetMode')
        targetMode = cfg.empiricalTargetMode;
    end
    target = empirical_mss_target(U10,targetMode);
    rho = cfg.mssResidualFraction;
    desiredAlong = target.along*(breakingRawTotal.along/target.along)^rho;
    desiredCross = target.cross*(breakingRawTotal.cross/target.cross)^rho;

    % Preserve the raw primary/short-wave partition in each direction.  The
    % surface-carrying primary band is dressed exactly; unresolved octave
    % tiles receive the remaining positive moment budget.
    primaryAlongTarget = desiredAlong*breakingMssRaw.along/breakingRawTotal.along;
    primaryCrossTarget = desiredCross*breakingMssRaw.cross/breakingRawTotal.cross;
    Hbreaking = dress_directional_mss(Hbreaking,KX,KY, ...
        primaryAlongTarget,primaryCrossTarget);
    breakingMss = spectral_mss(Hbreaking,KX,KY);
    shortMssDressed = shortMss;
    shortMssDressed.along = max(desiredAlong-breakingMss.along,0);
    shortMssDressed.cross = max(desiredCross-breakingMss.cross,0);
    shortMssDressed.total = shortMssDressed.along+shortMssDressed.cross;
else
    breakingMss = breakingMssRaw;
    shortMssDressed = shortMss;
end

realization.linear = linearTotal;
realization.breaking = add_mss(breakingMss,shortMssDressed);
realization.primaryLinear = linearMss;
realization.primaryBreaking = breakingMss;
realization.shortWave = shortMss;
realization.shortWaveBreaking = shortMssDressed;
realization.breakingRaw = breakingRawTotal;
realization.peakWavenumber = meta.peakWavenumber;
realization.primaryDomainLength = 2*pi/dk;
realization.primarySpacing = realization.primaryDomainLength/N;
realization.rawShortWaveCoefficient = meta.rawShortWaveCoefficient;
realization.linearSurface = real(ifft2(Hlinear))*N^2;
realization.breakingSurface = real(ifft2(Hbreaking))*N^2;
end

function H = sample_hermitian_coefficients(Psi,dkx,dky)
[Ny,Nx] = size(Psi);
white = (randn(Ny,Nx)+1i*randn(Ny,Nx))/sqrt(2);
negativeX = [1,Nx:-1:2];
negativeY = [1,Ny:-1:2];
white = (white+conj(white(negativeY,negativeX)))/sqrt(2);
H = sqrt(Psi*dkx*dky).*white;
end

function Hnonlinear = apply_thesis_modified_lie(H,KX,KY,K,~,kp,cfg)
Ksafe = max(K,realmin);
inputMask = K <= cfg.lieInputPeakMultiple*kp;
outputMask = K <= cfg.lieOutputPeakMultiple*kp;
N = size(H,1);

% h_tx and h_ty in the Lie/Creamer formulation are the two Riesz
% (multidimensional Hilbert) components of elevation, not spatial slopes.
% Treating them as derivatives introduces an erroneous k^2 amplification.
Htx = -1i*(KX./Ksafe).*H.*inputMask;
Hty = -1i*(KY./Ksafe).*H.*inputMask;
htx = real(ifft2(Htx))*N^2;
hty = real(ifft2(Hty))*N^2;
Fxx = fft2(htx.^2)/N^2;
Fxy = fft2(htx.*hty)/N^2;
Fyy = fft2(hty.^2)/N^2;

% Dimensionally consistent second-order Lie/Creamer term.  Every multiplier
% is even under k -> -k, so Hermitian symmetry is retained without replacing
% signed physics by U10*abs(k_i)/k.
Lstar = -(KX.^2./(2*Ksafe)).*Fxx ...
    -(KX.*KY./Ksafe).*Fxy ...
    -(KY.^2./(2*Ksafe)).*Fyy;
Lstar = cfg.modifiedLieScale*Lstar.*outputMask;
Lstar(K == 0) = 0;
Hnonlinear = H+Lstar;

negative = [1,N:-1:2];
Hnonlinear = (Hnonlinear+conj(Hnonlinear(negative,negative)))/2;
end

function Hdressed = dress_directional_mss(H,KX,KY,targetAlong,targetCross)
% Positive angular spectral dressing: |G(phi)|^2=exp(l0+l2*cos(2phi)).
% Two Newton variables enforce the requested along/cross slope moments while
% preserving every Fourier phase and hence the nonlinear crest geometry.
K2 = KX.^2+KY.^2;
cos2phi = (KX.^2-KY.^2)./max(K2,realmin);
power = abs(H).^2;
lambda = [0;0];
target = [targetAlong;targetCross];
for iteration = 1:20
    weight2 = exp(lambda(1)+lambda(2)*cos2phi);
    moments = [sum(KX.^2.*power.*weight2,'all'); ...
               sum(KY.^2.*power.*weight2,'all')];
    residual = moments-target;
    if max(abs(residual)./max(target,realmin)) < 1e-10
        break
    end
    jacobian = [moments(1),sum(KX.^2.*power.*weight2.*cos2phi,'all'); ...
                moments(2),sum(KY.^2.*power.*weight2.*cos2phi,'all')];
    step = jacobian\residual;
    lambda = lambda-step;
end
Hdressed = H.*exp(0.5*(lambda(1)+lambda(2)*cos2phi));
Hdressed(K2 == 0) = 0;
end

function stats = spectral_mss(H,KX,KY)
stats.along = sum(KX.^2.*abs(H).^2,'all');
stats.cross = sum(KY.^2.*abs(H).^2,'all');
stats.total = stats.along+stats.cross;
stats.gamma = sqrt(stats.cross/max(stats.along,realmin));
end

function result = add_mss(primary,shortWave)
result.along = primary.along+shortWave.along;
result.cross = primary.cross+shortWave.cross;
result.total = result.along+result.cross;
result.gamma = sqrt(result.cross/max(result.along,realmin));
end
