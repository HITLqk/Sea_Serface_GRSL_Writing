clear; close all; clc;

cfg = default_elfouhaily_ideal_curl_config();
surfaceData = generate_elfouhaily_ideal_curl_surface(cfg);

assert(surfaceData.metrics.overturningPointCount > 0, ...
    'The local patch is deformed but does not geometrically overturn.');
assert(surfaceData.metrics.maxOutsideDisplacement < 1e-12, ...
    'The deformation leaks outside the compact transition support.');
assert(surfaceData.metrics.downwardToLocalHeight < 0.90, ...
    'The downward displacement is too large for the local wave height.');
assert(surfaceData.metrics.maxUpwardDisplacement < ...
    0.10*surfaceData.localWaveHeight, ...
    'The curl incorrectly raises the crest instead of extending forward.');
assert(surfaceData.metrics.forwardOverturningFraction > 0.80, ...
    'The overturned material is not concentrated on the forward side.');
assert(surfaceData.metrics.lowestLipRelativeU > 0, ...
    'The lowest part of the lip is not in front of the detected crest.');
assert(surfaceData.metrics.crestwiseRidgeStd > 0, ...
    'The detected crest ridge has collapsed to a straight extrusion.');

if ~exist(cfg.output.outputDirectory,'dir')
    mkdir(cfg.output.outputDirectory);
end

fprintf('Local curled Elfouhaily surface generated.\n');
fprintf('  height quantile threshold : %.4f m\n', ...
    surfaceData.detection.heightThreshold);
fprintf('  slope tolerance           : %.5f\n', ...
    surfaceData.detection.slopeTolerance);
fprintf('  curvature threshold       : %.5f 1/m\n', ...
    surfaceData.detection.curvatureThreshold);
fprintf('  joint candidates          : %d\n', ...
    surfaceData.detection.candidateCount);
fprintf('  selected crest            : (%.3f, %.3f, %.3f) m\n', ...
    surfaceData.detection.x,surfaceData.detection.y, ...
    surfaceData.detection.z);
fprintf('  estimated local wave H    : %.4f m\n', ...
    surfaceData.localWaveHeight);

style = publication_style();
psi = deg2rad(cfg.detection.propagationDirectionDeg);

%% Common quantities
zAll = [surfaceData.Z0(:); surfaceData.Z(:)];
zAbs = max(abs(zAll));
surfaceCLim = 1.02*[-zAbs zAbs];

% Local crop used in the 3-D close-up panel
uAxis = [-1.15 1.25];
vAxis = [-0.95 0.95];
localCrop = surfaceData.localU >= uAxis(1) & surfaceData.localU <= uAxis(2) & ...
    surfaceData.localV >= vAxis(1) & surfaceData.localV <= vAxis(2);
Xloc = surfaceData.X; Yloc = surfaceData.Y; Zloc = surfaceData.Z;
Xloc(~localCrop) = NaN; Yloc(~localCrop) = NaN; Zloc(~localCrop) = NaN;
localZVals = surfaceData.Z(localCrop);

% Crop rectangle corners for global 3-D context
cornersU = [uAxis(1) uAxis(2) uAxis(2) uAxis(1) uAxis(1)];
cornersV = [vAxis(1) vAxis(1) vAxis(2) vAxis(2) vAxis(1)];
cornersX = surfaceData.detection.x + cos(psi).*cornersU - sin(psi).*cornersV;
cornersY = surfaceData.detection.y + sin(psi).*cornersU + cos(psi).*cornersV;
rectZ = max(surfaceData.Z0(localCrop)) + 0.01;

%% Figure 1: crest-placement conditions (2x2)
fig1 = new_pub_figure(cfg.output.figureVisible,7.10,5.75);
tl1 = tiledlayout(fig1,2,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl1,1);
plot_condition_mask(ax1,surfaceData.X0,surfaceData.Y0, ...
    surfaceData.detection.heightMask,surfaceData,style,'(a) Height');

ax2 = nexttile(tl1,2);
plot_condition_mask(ax2,surfaceData.X0,surfaceData.Y0, ...
    surfaceData.detection.slopeMask,surfaceData,style,'(b) Slope');

ax3 = nexttile(tl1,3);
plot_condition_mask(ax3,surfaceData.X0,surfaceData.Y0, ...
    surfaceData.detection.curvatureMask,surfaceData,style,'(c) Curvature');

ax4 = nexttile(tl1,4);
scoreMap = surfaceData.detection.score;
validScore = isfinite(scoreMap);
imagesc(ax4,surfaceData.X0(1,:),surfaceData.Y0(:,1),scoreMap);
set(ax4,'YDir','normal'); axis(ax4,'image'); hold(ax4,'on');
colormap(ax4,style.scoreMap);
contour(ax4,surfaceData.X0,surfaceData.Y0,surfaceData.detection.candidateMask, ...
    [0.5 0.5],'Color',style.candidateColor,'LineWidth',1.1);
contour(ax4,surfaceData.X0,surfaceData.Y0,surfaceData.transitionWeight, ...
    [0.02 0.02],'Color',style.supportColor,'LineWidth',1.1);
contour(ax4,surfaceData.X0,surfaceData.Y0,surfaceData.detection.Zsmooth,6, ...
    'Color',style.contourColor,'LineWidth',0.40);
plot(ax4,surfaceData.detection.x,surfaceData.detection.y,'p', ...
    'MarkerFaceColor',style.markerFill,'MarkerEdgeColor',style.markerEdge, ...
    'MarkerSize',11,'LineWidth',1.0);
xlabel(ax4,'$x$ (m)'); ylabel(ax4,'$y$ (m)');
add_panel_label(ax4,'(d) Joint score');
if any(validScore(:))
    caxis(ax4,[min(scoreMap(validScore)) max(scoreMap(validScore))]);
end

% Create centered legend beneath the tiledlayout.
hSel = plot(ax4,nan,nan,'p','MarkerFaceColor',style.markerFill, ...
    'MarkerEdgeColor',style.markerEdge,'MarkerSize',10,'LineWidth',1.0);
hCand = plot(ax4,nan,nan,'-','Color',style.candidateColor,'LineWidth',1.1);
hSupp = plot(ax4,nan,nan,'-','Color',style.supportColor,'LineWidth',1.1);
lgd1 = legend(ax4,[hSel,hCand,hSupp], ...
    {'Selected crest $oldsymbol{x}_c$','Joint-candidate boundary', ...
     'Curl-support boundary'}, ...
    'Orientation','horizontal','Box','off');
lgd1.Layout.Tile = 'south';
lgd1.ItemTokenSize = [18 8];

export_pub_figure(fig1,fullfile(cfg.output.outputDirectory, ...
    '01_curl_detection_conditions'),true,cfg.output);

%% Figure 2: global/local 3-D two-panel view
fig2 = new_pub_figure(cfg.output.figureVisible,7.10,3.45);
tl2 = tiledlayout(fig2,1,2,'TileSpacing','compact','Padding','compact');

ax5 = nexttile(tl2,1);
surf(ax5,surfaceData.X,surfaceData.Y,surfaceData.Z, ...
    surfaceData.Z,'EdgeColor','none');
hold(ax5,'on');
plot3(ax5,surfaceData.detection.x,surfaceData.detection.y, ...
    surfaceData.detection.z,'p','MarkerFaceColor',style.markerFill, ...
    'MarkerEdgeColor',style.markerEdge,'MarkerSize',11,'LineWidth',1.0);
plot3(ax5,cornersX,cornersY,rectZ*ones(size(cornersX)),'-', ...
    'Color',style.supportColor,'LineWidth',1.0);
axis(ax5,'tight');
pbaspect(ax5,[1 1 0.30]);
view(ax5,37,26); grid(ax5,'on');
box(ax5,'on');
colormap(ax5,style.surfaceMap); caxis(ax5,surfaceCLim);
xlabel(ax5,'$x$ (m)'); ylabel(ax5,'$y$ (m)'); zlabel(ax5,'$z$ (m)');
add_panel_label(ax5,'(a) Global sea surface');

ax6 = nexttile(tl2,2);
surf(ax6,Xloc,Yloc,Zloc,Zloc,'EdgeColor','none');
hold(ax6,'on');
mesh(ax6,Xloc,Yloc,Zloc,'EdgeColor',style.meshColor, ...
    'FaceAlpha',0,'EdgeAlpha',0.18);
plot3(ax6,surfaceData.detection.x,surfaceData.detection.y, ...
    surfaceData.detection.z,'p','MarkerFaceColor',style.markerFill, ...
    'MarkerEdgeColor',style.markerEdge,'MarkerSize',11,'LineWidth',1.0);
axis(ax6,'tight');
pbaspect(ax6,[1.6 1.0 0.75]);
view(ax6,-18,14); grid(ax6,'on');
box(ax6,'on');
colormap(ax6,style.surfaceMap); caxis(ax6,surfaceCLim);
xlabel(ax6,'$x$ (m)'); ylabel(ax6,'$y$ (m)'); zlabel(ax6,'$z$ (m)');
add_panel_label(ax6,'(b) Local curled geometry');
zlim(ax6,[min(localZVals)-0.02, max(localZVals)+0.02]);
cb2 = colorbar(ax6);
cb2.Label.String = 'Surface elevation (m)';
cb2.Label.Interpreter = 'latex';
cb2.Box = 'off';

export_pub_figure(fig2,fullfile(cfg.output.outputDirectory, ...
    '02_curl_global_local_3d'),false,cfg.output);

%% Figure 3: center propagation-direction section
fig3 = new_pub_figure(cfg.output.figureVisible,3.45,2.70);
ax7 = axes(fig3); hold(ax7,'on'); box(ax7,'on'); grid(ax7,'on');

strip = abs(surfaceData.localV) <= 0.55*max(cfg.domain.dx,cfg.domain.dy);
[uBase,order] = sort(surfaceData.localU(strip));
zBase = surfaceData.Z0(strip);
zPre = surfaceData.zPre(strip);
uCurl = surfaceData.localUFinal(strip);
zCurl = surfaceData.Z(strip);

plot(ax7,uBase,zBase(order),'-','Color',style.baseColor,'LineWidth',1.25);
plot(ax7,uBase,zPre(order),'--','Color',style.preColor,'LineWidth',1.20);
plot(ax7,uCurl(order),zCurl(order),'-','Color',style.curlColor,'LineWidth',1.70);
plot(ax7,0,surfaceData.detection.z,'p','MarkerFaceColor',style.markerFill, ...
    'MarkerEdgeColor',style.markerEdge,'MarkerSize',8,'LineWidth',0.9);

xlim(ax7,[-1.35 1.45]);
yData = [zBase(:); zPre(:); zCurl(:); surfaceData.detection.z];
yPad = 0.08*max(range(yData),0.10);
ylim(ax7,[min(yData)-yPad, max(yData)+yPad]);
xlabel(ax7,'Propagation coordinate $u$ (m)');
ylabel(ax7,'Elevation $z$ (m)');
lgd3 = legend(ax7,{'Original crest','Pre-shaped crest', ...
    'Curled crest','Selected crest'}, ...
    'Location','southoutside','Orientation','horizontal','Box','off');
lgd3.ItemTokenSize = [16 8];

export_pub_figure(fig3,fullfile(cfg.output.outputDirectory, ...
    '03_curl_center_profile'),true,cfg.output);

if cfg.output.saveMatFile
    save(fullfile(cfg.output.outputDirectory,'elfouhaily_local_curl.mat'), ...
        'surfaceData','-v7.3');
end

%% -------------------------------------------------------------------
function plot_condition_mask(ax,X,Y,mask,surfaceData,style,panelText)
imagesc(ax,X(1,:),Y(:,1),double(mask));
set(ax,'YDir','normal'); axis(ax,'image'); hold(ax,'on');
colormap(ax,style.maskMap);
contour(ax,X,Y,surfaceData.detection.Zsmooth,6, ...
    'Color',style.contourColor,'LineWidth',0.40);
plot(ax,surfaceData.detection.x,surfaceData.detection.y,'p', ...
    'MarkerFaceColor',style.markerFill,'MarkerEdgeColor',style.markerEdge, ...
    'MarkerSize',11,'LineWidth',1.0);
xlabel(ax,'$x$ (m)'); ylabel(ax,'$y$ (m)');
add_panel_label(ax,panelText);
end

function fig = new_pub_figure(visibility,widthIn,heightIn)
fig = figure('Visible',visibility,'Color','w');
fig.Units = 'inches';
fig.Position = [0.8 0.8 widthIn heightIn];
set(fig,'Renderer','painters');
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',8.5);
set(groot,'defaultAxesLineWidth',0.75);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
end

function export_pub_figure(fig,baseFile,vectorPreferred,outCfg)
if nargin < 4
    outCfg.exportPdf = true;
    outCfg.exportPng = true;
end
if outCfg.exportPdf
    if vectorPreferred
        exportgraphics(fig,[baseFile '.pdf'],'ContentType','vector');
    else
        exportgraphics(fig,[baseFile '.pdf'],'ContentType','image','Resolution',450);
    end
end
if outCfg.exportPng
    exportgraphics(fig,[baseFile '.png'],'Resolution',450);
end
end

function add_panel_label(ax,labelText)
text(ax,0.02,0.98,labelText,'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontWeight','bold','FontSize',9, ...
    'BackgroundColor','w','Margin',1.0);
end

function style = publication_style()
style.maskMap = [0.97 0.97 0.98; 0.25 0.49 0.77];
style.scoreMap = parula(256);
style.surfaceMap = parula(256);
style.baseColor = [0.30 0.30 0.30];
style.preColor = [0.12 0.45 0.72];
style.curlColor = [0.84 0.22 0.14];
style.candidateColor = [0.08 0.38 0.68];
style.supportColor = [0.86 0.56 0.08];
style.contourColor = [0.62 0.62 0.62];
style.markerFill = [0.99 0.93 0.35];
style.markerEdge = [0.10 0.10 0.10];
style.meshColor = [0.18 0.18 0.18];
end
