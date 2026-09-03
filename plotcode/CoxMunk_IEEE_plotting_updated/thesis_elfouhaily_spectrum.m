function [Psi,meta] = thesis_elfouhaily_spectrum(K,KX,KY,U10,age,windDirectionDeg)
%THESIS_ELFOUHAILY_SPECTRUM Unified directional Elfouhaily elevation PSD.
%   Psi is the two-dimensional elevation spectrum in m^4.

g = 9.81;
km = 370;
cm = 0.23;
% Wind-dependent neutral drag used by standard Elfouhaily implementations.
% The previous constant 1.44e-3 over-energized short waves below about 10 m/s.
Cd10N = (0.8+0.065*U10)*1e-3;
uStar = sqrt(Cd10N)*U10;
kp = g*(age/U10)^2;
cp = sqrt(g/kp);
sigma = 0.08*(1+4*age^(-3));
alphaP = 0.006*sqrt(age);

if uStar <= cm
    alphaM = 0.01*(1+log(uStar/cm));
else
    alphaM = 0.01*(1+3*log(uStar/cm));
end
% The published equilibrium-range coefficient becomes negative below
% approximately 2.3 m/s. A negative PSD is inadmissible, so only that
% short-wave term is set to zero; the event is retained in meta.
rawAlphaM = alphaM;
alphaM = max(alphaM,0);

if age <= 1
    gamma = 1.7;
else
    gamma = 1.7+6*log(age);
end

Ksafe = max(K,realmin);
c = sqrt((g./Ksafe).*(1+(Ksafe/km).^2));
Lpm = exp(-1.25*(kp./Ksafe).^2);
Gamma = exp(-0.5*((sqrt(Ksafe/kp)-1)/sigma).^2);
Jp = gamma.^Gamma;
Fp = Lpm.*Jp.*exp(-age/sqrt(10).*(sqrt(Ksafe/kp)-1));
Fm = Lpm.*Jp.*exp(-0.25*(Ksafe/km-1).^2);
Bl = 0.5*alphaP*(cp./c).*Fp;
Bh = 0.5*alphaM*(cm./c).*Fm;
S = (Bl+Bh)./Ksafe.^3;

delta = tanh(log(2)/4 + 4*(c/cp).^2.5 + ...
    0.13*(uStar/cm)*(cm./c).^2.5);
phi = atan2(KY,KX);
phiWind = deg2rad(windDirectionDeg);
spread = (1+delta.*cos(2*(phi-phiWind)))/(2*pi);
Psi = max(S.*spread./Ksafe,0);
Psi(K == 0) = 0;

meta.peakWavenumber = kp;
meta.frictionVelocity = uStar;
meta.rawShortWaveCoefficient = rawAlphaM;
meta.shortWaveCoefficient = alphaM;
end
