function target = empirical_mss_target(U10,mode)
%EMPIRICAL_MSS_TARGET Component-wise fusion of independent slope laws.
%   The target is not an additional observation.  It is the arithmetic mean
%   of available Cox-Munk, Hu/TGRS and Guerin along/cross MSS values.  Using
%   components (rather than total MSS alone) also constrains anisotropy.

if nargin < 2
    mode = 'fusion';
end

if U10 < 7
    huTotal = 14.6e-3*sqrt(U10);
elseif U10 < 13.3
    huTotal = 3e-3+5.12e-3*U10;
else
    huTotal = 138e-3*log10(U10)-84e-3;
end
huGamma = 0.864;
huAlong = huTotal/(1+huGamma^2);
huCross = huTotal-huAlong;

if strcmpi(mode,'hu-high-wind')
    target.along = huAlong;
    target.cross = huCross;
    target.total = huTotal;
    target.gamma = huGamma;
    return
elseif ~strcmpi(mode,'fusion')
    error('Unknown empirical target mode: %s',mode);
end

U12p5 = U10/0.98;
along = 3.16e-3*U12p5;
cross = 3e-3+1.92e-3*U12p5;
alongValues = [along,huAlong];
crossValues = [cross,huCross];

[guerinAlong,guerinCross] = local_guerin(U10);
if isfinite(guerinAlong)
    alongValues(end+1) = guerinAlong;
    crossValues(end+1) = guerinCross;
end

target.along = mean(alongValues);
target.cross = mean(crossValues);
target.total = target.along+target.cross;
target.gamma = sqrt(target.cross/target.along);
end

function [along,cross] = local_guerin(U10)
grid = (3:0.5:15)';
alongData = 0.01*[1.10 1.16 1.24 1.34 1.46 1.59 1.74 1.91 2.09 ...
    2.28 2.49 2.68 2.86 3.05 3.23 3.42 3.60 3.79 3.98 4.13 ...
    4.40 4.46 4.52 4.75 4.86]';
crossData = 0.01*[0.97 1.03 1.11 1.20 1.28 1.37 1.45 1.53 1.62 ...
    1.70 1.79 1.88 1.98 2.08 2.20 2.33 2.46 2.59 2.72 2.83 ...
    3.00 3.06 3.10 3.26 3.36]';
along = interp1(grid,alongData,U10,'linear',NaN);
cross = interp1(grid,crossData,U10,'linear',NaN);
end
