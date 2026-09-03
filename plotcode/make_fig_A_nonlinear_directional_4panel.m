%% make_fig_A_nonlinear_directional_4panel.m
% Paper-ready 4-panel figure for
% Section II.A "Nonlinear Background and Directional Wind Waves"
%
% Assumed function files in the same project:
%   - default_nonlinear_lie_config.m
%   - generate_nonlinear_lie_elfouhaily_surface.m
%   - default_wind_components_config.m
%   - generate_directional_wind_components_surface.m
%
% If your local files still carry suffixes such as "(1)" or "(2)",
% please rename them so that each filename matches its function name.
%
% The 4 panels are:
%   (a) Linear vs nonlinear background profile
%   (b) Removed short-wave band vs equal-energy directional replacement
%   (c) Final sea surface top view after directional short-wave replacement
%   (d) Local 3-D view around the selected steep crest
%
% No subplot titles are used, so the figure can be inserted into the paper
% directly and explained entirely in the figure caption.

clear; clc; close all;

%% -------------------- user-adjustable settings ----------------------- %%
outputDir = fullfile(pwd,'output');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

% -------- visual compression controls for paper figure --------
topViewColorExpand = 1.20;   % >1 enlarges color range, making relief look milder
localZAspect       = 0.30;   % smaller value -> flatter local 3D appearance
localViewAz        = 42;
localViewEl        = 24;

% -------- optional physical amplitude reduction --------------
shortWaveEnergyRatio = 0.75; % set to 1.0 to keep original energy


figureName = 'Fig_IIA_nonlinear_directional_4panel';
exportDPI  = 300;

% Local 3-D window half-width (metres)
detailHalfWidth = 14.0;

% Section length for the two line-profile panels (metres)
sectionHalfLength = 20.0;
sectionSamples    = 800;

%% -------------------- global plotting defaults ---------------------- %%
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',10.5);
set(groot,'defaultAxesLineWidth',0.8);
set(groot,'defaultLineLineWidth',1.4);
set(groot,'defaultAxesBox','on');
set(groot,'defaultAxesLayer','top');
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% -------------------- Step 1: nonlinear background ------------------ %%
cfgNL = default_nonlinear_lie_config();
cfgNL.output.figureVisible  = 'off';
cfgNL.output.outputDirectory = outputDir;
cfgNL.output.saveSurfaceMat = false;

% Optional: adjust these to match the exact case used in your manuscript.
% cfgNL.sea.U10 = 10.0;
% cfgNL.sea.windDirectionDeg = 0.0;
% cfgNL.domain.Lx = 128.0; cfgNL.domain.Ly = 128.0;
% cfgNL.domain.dx = 0.25;  cfgNL.domain.dy = 0.25;

nlData = generate_nonlinear_lie_elfouhaily_surface(cfgNL);

% Save a temporary background MAT file for the directional short-wave step.
backgroundMatFile = fullfile(outputDir,'background_for_directional_shortwaves.mat');
X       = nlData.X;
Y       = nlData.Y;
Z       = nlData.Z;
XLinear = nlData.X0;
YLinear = nlData.Y0;
ZLinear = nlData.ZLinear;
cfg     = nlData.cfg; %#ok<NASGU>
save(backgroundMatFile,'X','Y','Z','XLinear','YLinear','ZLinear','cfg','-v7.3');

%% -------------------- Step 2: directional short waves -------------- %%
cfgW = default_wind_components_config();
cfgW.backgroundMatFile      = backgroundMatFile;
cfgW.output.figureVisible   = 'off';
cfgW.output.outputDirectory = outputDir;
cfgW.output.saveSurfaceMat  = false;
cfgW.wind.replacementEnergyRatio = shortWaveEnergyRatio;

% Recommended defaults already match the manuscript text:
%   lambda_s,min = 1.5 m
%   lambda_s,max = 5.0 m
%   sigma_theta  = 8 deg
%   rho_s        = 1.0
% Adjust only if your paper uses another case.

windData = generate_directional_wind_components_surface(cfgW);

%% -------------------- common quantities ----------------------------- %%
psi = deg2rad(cfgW.wind.propagationDirectionDeg);

% Selected crest / local window centre from the directional short-wave step
xc = windData.detection.parameterX;
yc = windData.detection.parameterY;
zc = windData.detection.z;

% 1-D section along the propagation direction through the selected crest
s = linspace(-sectionHalfLength,sectionHalfLength,sectionSamples);
xSec = xc + s*cos(psi);
ySec = yc + s*sin(psi);

zLinearProfile = interp2(nlData.X0,nlData.Y0,nlData.ZLinear, ...
    xSec,ySec,'linear',NaN);
zNonlinearProfile = interp2(windData.XLinear,windData.YLinear, ...
    windData.ZOriginalBackground,xSec,ySec,'linear',NaN);
zRemovedProfile = interp2(windData.XLinear,windData.YLinear, ...
    windData.ZRemovedBand,xSec,ySec,'linear',NaN);
zComponentProfile = interp2(windData.XLinear,windData.YLinear, ...
    windData.ZWindComponent,xSec,ySec,'linear',NaN);

removedRMS   = std(windData.ZRemovedBand,0,'all');
componentRMS = std(windData.ZWindComponent,0,'all');

% Local 3-D region indices
xVec = windData.XLinear(1,:);
yVec = windData.YLinear(:,1);
colIdx = find(abs(xVec-xc) <= detailHalfWidth);
rowIdx = find(abs(yVec-yc) <= detailHalfWidth);

XLocal = windData.X(rowIdx,colIdx);
YLocal = windData.Y(rowIdx,colIdx);
ZLocal = windData.Z(rowIdx,colIdx);

% Common elevation colour scale for panels (c) and (d)
zMin = min(windData.Z,[],'all');
zMax = max(windData.Z,[],'all');
commonCLim = [zMin zMax];

%% -------------------- Step 3: 4-panel paper figure ----------------- %%
fig = figure('Color','w','Position',[80 80 1180 860]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% -------------------- (a) Linear vs nonlinear profile -----------------
ax1 = nexttile(tl,1);
plot(ax1,s,zLinearProfile,'Color',[0.25 0.25 0.25],'LineWidth',1.35); hold(ax1,'on');
plot(ax1,s,zNonlinearProfile,'Color',[0.80 0.10 0.10],'LineWidth',1.55);
grid(ax1,'on');
box(ax1,'on');
xlabel(ax1,'Propagation coordinate $u$ (m)');
ylabel(ax1,'Surface elevation $\eta$ (m)');
legend(ax1,{'Linear background $\eta_L$','Nonlinear background $\eta_N$'}, ...
    'Location','northwest','Box','off');
add_panel_label(ax1,'(a)');

% -------- (b) Removed band vs equal-energy directional replacement ----
ax2 = nexttile(tl,2);
plot(ax2,s,zRemovedProfile,'Color',[0.10 0.35 0.75],'LineWidth',1.35); hold(ax2,'on');
plot(ax2,s,zComponentProfile,'Color',[0.85 0.35 0.10],'LineWidth',1.55);
grid(ax2,'on');
box(ax2,'on');
xlabel(ax2,'Propagation coordinate $u$ (m)');
ylabel(ax2,'Short-wave elevation (m)');
legend(ax2,{sprintf('Removed band, RMS = %.4f m',removedRMS), ...
            sprintf('Directional replacement, RMS = %.4f m',componentRMS)}, ...
    'Location','northwest','Box','off');
text(ax2,0.98,0.06,'Equal-energy replacement ($\rho_s=1$)', ...
    'Units','normalized','HorizontalAlignment','right', ...
    'VerticalAlignment','bottom','FontSize',9.5);
add_panel_label(ax2,'(b)');

% -------------------- (c) Final surface top view ----------------------
ax3 = nexttile(tl,3);
surf(ax3,windData.XLinear,windData.YLinear,zeros(size(windData.Z)), ...
    windData.Z,'EdgeColor','none');
view(ax3,2);
axis(ax3,'image'); axis(ax3,'tight');
xlabel(ax3,'$x$ (m)');
ylabel(ax3,'$y$ (m)');
colormap(ax3,turbo);
clim(ax3,commonCLim);
hold(ax3,'on');
plot(ax3,xc,yc,'kp','MarkerFaceColor','w','MarkerSize',9,'LineWidth',1.0);
rectangle(ax3,'Position',[xc-detailHalfWidth,yc-detailHalfWidth, ...
    2*detailHalfWidth,2*detailHalfWidth], ...
    'EdgeColor','k','LineStyle','--','LineWidth',0.95);
arrowLength = 10.0;
quiver(ax3,xc,yc,arrowLength*cos(psi),arrowLength*sin(psi),0, ...
    'k','LineWidth',1.2,'MaxHeadSize',0.8);
cb3 = colorbar(ax3);
cb3.Label.String = 'Surface elevation (m)';
cb3.Label.Interpreter = 'latex';
add_panel_label(ax3,'(c)');

% -------------------- (d) Local 3-D view ------------------------------
ax4 = nexttile(tl,4);
surf(ax4,XLocal,YLocal,ZLocal,ZLocal,'EdgeColor','none');
shading(ax4,'interp');
axis(ax4,'tight');
view(ax4,localViewAz,localViewEl);
pbaspect(ax4,[1 1 localZAspect]);
grid(ax4,'on');
camproj(ax4,'orthographic');
camlight(ax4,'headlight'); lighting(ax4,'gouraud');
material(ax4,'dull');
xlabel(ax4,'$x$ (m)');
ylabel(ax4,'$y$ (m)');
zlabel(ax4,'$z$ (m)');
colormap(ax4,turbo);
clim(ax4,commonCLim);
hold(ax4,'on');
plot3(ax4,windData.detection.x,windData.detection.y,zc, ...
    'kp','MarkerFaceColor','w','MarkerSize',9,'LineWidth',1.0);
quiver3(ax4,windData.detection.x,windData.detection.y,zc+0.08, ...
    4*cos(psi),4*sin(psi),0,0,'k','LineWidth',1.2,'MaxHeadSize',0.7);
add_panel_label(ax4,'(d)');

%% -------------------- export ---------------------------------------- %%
exportgraphics(fig,fullfile(outputDir,[figureName '.png']), ...
    'Resolution',exportDPI,'BackgroundColor','white');
exportgraphics(fig,fullfile(outputDir,[figureName '.pdf']), ...
    'ContentType','vector','BackgroundColor','white');

fprintf('4-panel figure exported to:\n');
fprintf('  %s\n',fullfile(outputDir,[figureName '.png']));
fprintf('  %s\n',fullfile(outputDir,[figureName '.pdf']));

%% -------------------- local helper functions ------------------------ %%
function add_panel_label(ax,labelText)
text(ax,0.02,0.98,labelText,'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontWeight','bold','FontSize',11, ...
    'BackgroundColor','w','Margin',1.0);
end
