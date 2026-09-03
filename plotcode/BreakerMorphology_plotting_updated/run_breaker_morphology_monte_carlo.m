clear; close all; clc;

cfg = default_breaker_morphology_validation_config();
assert(isfolder(cfg.curlDirectory), ...
    'Curl directory not found: %s',cfg.curlDirectory);
addpath(cfg.curlDirectory);
cleanupPath = onCleanup(@() rmpath(cfg.curlDirectory)); %#ok<NASGU>

rng(cfg.monteCarlo.randomSeed,'twister');
groupKeys = ["Original","Corrected"];
groupLabels = [string(cfg.labels.previous), string(cfg.labels.proposed)];
rangeSets = {cfg.monteCarlo.original,cfg.monteCarlo.corrected};
nPerGroup = cfg.monteCarlo.nPerGroup;
nRows = numel(groupKeys)*nPerGroup;

% ---------------------- storage -------------------------------------- %
group = strings(nRows,1);
groupLabel = strings(nRows,1);
iteration = zeros(nRows,1);
randomSeed = zeros(nRows,1);
amplitudeCurl = zeros(nRows,1);
curlMultiplier = zeros(nRows,1);
pivotDepth = zeros(nRows,1);
forwardGain = zeros(nRows,1);
verticalAngleRatio = zeros(nRows,1);
frontFaceAngleDeg = nan(nRows,1);
rxOverLambda = nan(nRows,1);
ryOverLambda = nan(nRows,1);
valid = false(nRows,1);
pass = false(nRows,1);
errorMessage = strings(nRows,1);

% ---------------------- Monte Carlo generation ----------------------- %
row = 0;
for g = 1:numel(groupKeys)
    ranges = rangeSets{g};
    for k = 1:nPerGroup
        row = row+1;
        group(row) = groupKeys(g);
        groupLabel(row) = groupLabels(g);
        iteration(row) = k;
        randomSeed(row) = randi([1 2^31-1]);
        amplitudeCurl(row) = draw_uniform(ranges.amplitudeCurl);
        curlMultiplier(row) = draw_uniform(ranges.curlMultiplier);
        pivotDepth(row) = draw_uniform(ranges.pivotDepth);
        forwardGain(row) = draw_uniform(ranges.forwardGain);
        verticalAngleRatio(row) = draw_uniform(ranges.verticalAngleRatio);

        curlCfg = default_elfouhaily_ideal_curl_config();
        curlCfg.randomSeed = randomSeed(row);
        curlCfg.curl.amplitudeCurl = amplitudeCurl(row);
        curlCfg.curl.curlMultiplier = curlMultiplier(row);
        curlCfg.curl.pivotDepth = pivotDepth(row);
        curlCfg.curl.forwardGain = forwardGain(row);
        curlCfg.curl.verticalAngleRatio = verticalAngleRatio(row);

        try
            surfaceData = generate_elfouhaily_ideal_curl_surface(curlCfg);
            metrics = extract_breaker_morphology_metrics(surfaceData,cfg);
            frontFaceAngleDeg(row) = metrics.frontFaceAngleDeg;
            rxOverLambda(row) = metrics.rxOverLambda;
            ryOverLambda(row) = metrics.ryOverLambda;
            valid(row) = true;
            pass(row) = metrics.allPass;
        catch ME
            errorMessage(row) = string(ME.message);
        end

        if mod(k,10) == 0
            fprintf('%s: %d/%d complete\n',groupKeys(g),k,nPerGroup);
        end
    end
end

raw = table(group,groupLabel,iteration,randomSeed,amplitudeCurl, ...
    curlMultiplier,pivotDepth,forwardGain,verticalAngleRatio, ...
    frontFaceAngleDeg,rxOverLambda,ryOverLambda,valid,pass,errorMessage);

summary = build_summary(raw,groupKeys,groupLabels);
disp(summary);

if ~isfolder(cfg.outputDirectory)
    mkdir(cfg.outputDirectory);
end
writetable(raw,fullfile(cfg.outputDirectory, ...
    'breaker_morphology_monte_carlo_raw.csv'));
writetable(summary,fullfile(cfg.outputDirectory, ...
    'breaker_morphology_monte_carlo_summary.csv'));
save(fullfile(cfg.outputDirectory, ...
    'breaker_morphology_monte_carlo.mat'),'cfg','raw','summary');
write_summary_table_tex(summary,cfg);

% ---------------------- publication figure --------------------------- %
fig = new_pub_figure(cfg.figureVisible,cfg.style.figureWidthIn, ...
    cfg.style.figureHeightIn,cfg.style);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl,1);
plot_joint_panel(ax1,raw,summary,cfg,'rx');
add_panel_label(ax1,'(a)');

ax2 = nexttile(tl,2);
plot_joint_panel(ax2,raw,summary,cfg,'ry');
add_panel_label(ax2,'(b)');

% Shared legend centered below both panels.
hRef = patch(nan,nan,cfg.style.envelopeFaceColor,'EdgeColor', ...
    cfg.style.envelopeEdgeColor,'LineWidth',0.8,'FaceAlpha',0.75);
hPrev = scatter(nan,nan,cfg.style.markerSize, ...
    cfg.style.previousColor,'o','LineWidth',0.9);
hProp = scatter(nan,nan,cfg.style.markerSize, ...
    cfg.style.proposedColor,'^','LineWidth',0.9);
hPass = scatter(nan,nan,cfg.style.markerSize+6, ...
    cfg.style.passColor,'o','filled','MarkerEdgeColor',cfg.style.passColor);
lgd = legend(tl,[hRef,hPrev,hProp,hPass], ...
    {cfg.labels.reference,cfg.labels.previous,cfg.labels.proposed, ...
    'Proposed--Joint-pass realization'}, ...
    'Orientation','horizontal','Box','off');
lgd.Layout.Tile = 'south';
lgd.ItemTokenSize = cfg.style.legendTokenSize;
lgd.FontSize = cfg.style.legendFontSize;

baseFile = fullfile(cfg.outputDirectory,cfg.output.figureBaseName);
export_pub_figure(fig,baseFile,cfg.output);

fprintf('Publication figure written to:\n  %s.pdf\n  %s.png\n', ...
    baseFile,baseFile);
fprintf('Publication table written to:\n  %s.tex\n', ...
    fullfile(cfg.outputDirectory,cfg.output.tableBaseName));

% ==================================================================== %
function value = draw_uniform(bounds)
value = bounds(1)+(bounds(2)-bounds(1))*rand();
end

function summary = build_summary(raw,groupKeys,groupLabels)
n = numel(groupKeys);
group = groupKeys(:);
groupLabel = groupLabels(:);
requestedCount = zeros(n,1);
validCount = zeros(n,1);
passCount = zeros(n,1);
passRate = zeros(n,1);
angleMean = nan(n,1); angleCi95 = nan(n,1);
rxMean = nan(n,1); rxCi95 = nan(n,1);
ryMean = nan(n,1); ryCi95 = nan(n,1);
angleMedian = nan(n,1); rxMedian = nan(n,1); ryMedian = nan(n,1);
for k = 1:n
    selected = raw.group == group(k);
    accepted = selected & raw.valid;
    requestedCount(k) = nnz(selected);
    validCount(k) = nnz(accepted);
    passCount(k) = nnz(accepted & raw.pass);
    passRate(k) = passCount(k)/max(validCount(k),1);
    [angleMean(k),angleCi95(k)] = mean_ci(raw.frontFaceAngleDeg(accepted));
    [rxMean(k),rxCi95(k)] = mean_ci(raw.rxOverLambda(accepted));
    [ryMean(k),ryCi95(k)] = mean_ci(raw.ryOverLambda(accepted));
    angleMedian(k) = median(raw.frontFaceAngleDeg(accepted),'omitnan');
    rxMedian(k) = median(raw.rxOverLambda(accepted),'omitnan');
    ryMedian(k) = median(raw.ryOverLambda(accepted),'omitnan');
end
summary = table(group,groupLabel,requestedCount,validCount,passCount, ...
    passRate,angleMean,angleCi95,rxMean,rxCi95,ryMean,ryCi95, ...
    angleMedian,rxMedian,ryMedian);
end

function [average,halfWidth] = mean_ci(values)
values = values(isfinite(values));
assert(~isempty(values),'No finite values are available.');
average = mean(values);
halfWidth = 1.96*std(values)/sqrt(max(numel(values),1));
end

function plot_joint_panel(ax,raw,summary,cfg,whichScale)
hold(ax,'on');
set(ax,'FontName',cfg.style.fontName, ...
    'FontSize',cfg.style.fontSize, ...
    'LineWidth',cfg.style.axisLineWidth, ...
    'TickDir','out','Box','on','Layer','top');
grid(ax,'on');
ax.GridAlpha = cfg.style.gridAlpha;
ax.MinorGridAlpha = 0.5*cfg.style.gridAlpha;

angleBounds = cfg.reference.frontFaceAngleDeg;
switch lower(whichScale)
    case 'rx'
        scaleBounds = cfg.reference.rxOverLambda;
        yLabelText = '$r_x/\lambda_p$';
        yData = raw.rxOverLambda;
        panelMedianField = 'rxMedian';
        yPadding = 0.003;
    case 'ry'
        scaleBounds = cfg.reference.ryOverLambda;
        yLabelText = '$r_y/\lambda_p$';
        yData = raw.ryOverLambda;
        panelMedianField = 'ryMedian';
        yPadding = 0.003;
    otherwise
        error('Unknown panel selector: %s',whichScale);
end

% Reference envelope.
patch(ax,angleBounds([1 2 2 1]),scaleBounds([1 1 2 2]), ...
    cfg.style.envelopeFaceColor,'EdgeColor',cfg.style.envelopeEdgeColor, ...
    'LineStyle','-','LineWidth',0.8,'FaceAlpha',0.75);

% Previous group.
prev = raw.group == "Original" & raw.valid;
scatter(ax,raw.frontFaceAngleDeg(prev),yData(prev),cfg.style.markerSize, ...
    cfg.style.previousColor,'o','LineWidth',0.85, ...
    'MarkerFaceColor','none');

% Proposed group: all valid realizations.
prop = raw.group == "Corrected" & raw.valid;
scatter(ax,raw.frontFaceAngleDeg(prop),yData(prop),cfg.style.markerSize, ...
    cfg.style.proposedColor,'^','LineWidth',0.85, ...
    'MarkerFaceColor','none');

% Highlight proposed realizations that jointly pass all three constraints.
propPass = prop & raw.pass;
scatter(ax,raw.frontFaceAngleDeg(propPass),yData(propPass), ...
    cfg.style.markerSize+6,cfg.style.passColor,'o','filled', ...
    'MarkerEdgeColor',cfg.style.passColor);

% Median markers (no separate legend entry).
prevRow = find(summary.group == "Original",1);
propRow = find(summary.group == "Corrected",1);
plot(ax,summary.angleMedian(prevRow),summary.(panelMedianField)(prevRow), ...
    's','MarkerSize',5.6,'LineWidth',0.9, ...
    'MarkerEdgeColor',cfg.style.previousMedianColor, ...
    'MarkerFaceColor','w');
plot(ax,summary.angleMedian(propRow),summary.(panelMedianField)(propRow), ...
    'd','MarkerSize',5.6,'LineWidth',0.9, ...
    'MarkerEdgeColor',cfg.style.proposedMedianColor, ...
    'MarkerFaceColor','w');

xlabel(ax,'Front-face angle $\theta_f$ (deg)');
ylabel(ax,yLabelText);

% Focus the limits tightly around the observed range and reference box.
validY = yData(raw.valid & isfinite(yData));
yMin = min([validY(:); scaleBounds(:)]) - yPadding;
yMax = max([validY(:); scaleBounds(:)]) + yPadding;
xMin = min([raw.frontFaceAngleDeg(raw.valid); angleBounds(:)]) - 3.5;
xMax = max([raw.frontFaceAngleDeg(raw.valid); angleBounds(:)]) + 3.5;
xlim(ax,[xMin xMax]);
ylim(ax,[max(0,yMin) yMax]);
end

function fig = new_pub_figure(visibility,widthIn,heightIn,style)
fig = figure('Visible',visibility,'Color','w');
fig.Units = 'inches';
fig.Position = [0.8 0.8 widthIn heightIn];
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 widthIn heightIn];
set(fig,'Renderer','painters');
set(groot,'defaultAxesFontName',style.fontName);
set(groot,'defaultTextFontName',style.fontName);
set(groot,'defaultAxesFontSize',style.fontSize);
set(groot,'defaultAxesLineWidth',style.axisLineWidth);
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
end

function export_pub_figure(fig,baseFile,outCfg)
if outCfg.exportPdf
    exportgraphics(fig,[baseFile '.pdf'],'ContentType','vector', ...
        'BackgroundColor','white');
end
if outCfg.exportPng
    exportgraphics(fig,[baseFile '.png'],'Resolution', ...
        outCfg.pngResolution,'BackgroundColor','white');
end
end

function add_panel_label(ax,labelText)
text(ax,0.03,0.97,labelText,'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontWeight','bold','FontSize',9, ...
    'BackgroundColor','w','Margin',1.0);
end

function write_summary_table_tex(summary,cfg)
outFile = fullfile(cfg.outputDirectory,[cfg.output.tableBaseName '.tex']);
fid = fopen(outFile,'w');
assert(fid>0,'Could not open table file for writing: %s',outFile);
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid,'%% Auto-generated by run_breaker_morphology_monte_carlo.m\n');
fprintf(fid,'\\begin{table}[!t]\n');
fprintf(fid,'    \\centering\n');
fprintf(fid,'    \\caption{Monte Carlo statistics of curling-wave morphology.}\n');
fprintf(fid,'    \\label{tab:breaker_morphology_statistics}\n');
fprintf(fid,'    \\begin{tabular}{lcccc}\n');
fprintf(fid,'        \\toprule\n');
fprintf(fid,'        Group & $\\theta_f$ (deg) & $r_x/\\lambda_p$ & $r_y/\\lambda_p$ & Joint pass \\\\ \n');
fprintf(fid,'        \\midrule\n');
for k = 1:height(summary)
    fprintf(fid,['        %s & %.2f $\\pm$ %.2f & %.5f $\\pm$ %.5f & ' ...
        '%.5f $\\pm$ %.5f & %.0f\\%% \\\\ \n'], ...
        escape_tex(summary.groupLabel(k)), ...
        summary.angleMean(k), summary.angleCi95(k), ...
        summary.rxMean(k), summary.rxCi95(k), ...
        summary.ryMean(k), summary.ryCi95(k), ...
        100*summary.passRate(k));
end
fprintf(fid,['        %s & %g--%g & %g--%g & %g--%g & --- \\\\ \n'], ...
    escape_tex(cfg.labels.reference), ...
    cfg.reference.frontFaceAngleDeg(1), cfg.reference.frontFaceAngleDeg(2), ...
    cfg.reference.rxOverLambda(1), cfg.reference.rxOverLambda(2), ...
    cfg.reference.ryOverLambda(1), cfg.reference.ryOverLambda(2));
fprintf(fid,'        \\bottomrule\n');
fprintf(fid,'    \\end{tabular}\n');
fprintf(fid,'\\end{table}\n');
end

function out = escape_tex(in)
out = char(in);
out = strrep(out,'\\','\\textbackslash{}');
out = strrep(out,'_','\\_');
out = strrep(out,'--','--');
end
