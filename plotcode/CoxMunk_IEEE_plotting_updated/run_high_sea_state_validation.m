function [raw,summary,reference,assessment] = run_high_sea_state_validation(cfg)
%RUN_HIGH_SEA_STATE_VALIDATION Validate nonlinear surfaces from 10 to 40 m/s.
%   The corrected curve is required to improve upon both the simulated linear
%   Elfouhaily surface and the integrated Elfouhaily spectrum against the
%   Hu/TGRS high-wind law.  It must also remain monotone and flatten at high U.

arguments
    cfg (1,1) struct = default_high_sea_state_config()
end
validate_high_config(cfg);
if ~isfolder(cfg.outputDirectory)
    mkdir(cfg.outputDirectory);
end

reference = two_group_mss_references(cfg.windSpeeds,cfg);
groups = ["Linear Elfouhaily","Raw Modified Lie", ...
    "MSS-Constrained Modified Lie"];
nRows = numel(cfg.windSpeeds)*numel(cfg.realizationSeeds)*numel(groups);
U10 = zeros(nRows,1); Seed = zeros(nRows,1); Group = strings(nRows,1);
MssAlong = zeros(nRows,1); MssCross = zeros(nRows,1);
MssTotal = zeros(nRows,1); Gamma = zeros(nRows,1);
row = 0;

for windIndex = 1:numel(cfg.windSpeeds)
    wind = cfg.windSpeeds(windIndex);
    for seedIndex = 1:numel(cfg.realizationSeeds)
        seed = cfg.realizationSeeds(seedIndex);
        result = synthesize_two_group_realization(wind,seed,cfg);
        append(groups(1),result.linear);
        append(groups(2),result.breakingRaw);
        append(groups(3),result.breaking);
    end
end

raw = table(U10,Seed,Group,MssAlong,MssCross,MssTotal,Gamma);
summary = summarize_high_results(raw,groups,cfg.windSpeeds);
assessment = assess_high_results(summary,reference,groups,cfg);

fig = plot_high_wind_mss(summary,reference,groups,cfg);
export_ieee_figure(fig,fullfile(cfg.outputDirectory, ...
    'Fig_high_sea_state_mss_validation'),cfg.exportResolution);
close(fig);
writetable(raw,fullfile(cfg.outputDirectory,'high_sea_state_raw.csv'));
writetable(summary,fullfile(cfg.outputDirectory,'high_sea_state_summary.csv'));
writetable(reference,fullfile(cfg.outputDirectory,'high_sea_state_reference.csv'));
writetable(assessment,fullfile(cfg.outputDirectory,'high_sea_state_assessment.csv'));
save(fullfile(cfg.outputDirectory,'high_sea_state_validation.mat'), ...
    'raw','summary','reference','assessment','cfg','-v7.3');
disp(assessment);

    function append(groupName,mss)
        row = row+1;
        U10(row) = wind;
        Seed(row) = seed;
        Group(row) = groupName;
        MssAlong(row) = mss.along;
        MssCross(row) = mss.cross;
        MssTotal(row) = mss.total;
        Gamma(row) = mss.gamma;
    end
end

function summary = summarize_high_results(raw,groups,winds)
n = numel(groups)*numel(winds);
Group = strings(n,1); U10 = zeros(n,1); AlongMedian = zeros(n,1);
CrossMedian = zeros(n,1); TotalMedian = zeros(n,1);
TotalQ05 = zeros(n,1); TotalQ95 = zeros(n,1); GammaMedian = zeros(n,1);
row = 0;
for groupIndex = 1:numel(groups)
    for windIndex = 1:numel(winds)
        row = row+1;
        use = raw.Group == groups(groupIndex) & raw.U10 == winds(windIndex);
        values = sort(raw.MssTotal(use));
        Group(row) = groups(groupIndex); U10(row) = winds(windIndex);
        AlongMedian(row) = median(raw.MssAlong(use));
        CrossMedian(row) = median(raw.MssCross(use));
        TotalMedian(row) = median(values);
        TotalQ05(row) = local_quantile(values,0.05);
        TotalQ95(row) = local_quantile(values,0.95);
        GammaMedian(row) = median(raw.Gamma(use));
    end
end
summary = table(Group,U10,AlongMedian,CrossMedian,TotalMedian, ...
    TotalQ05,TotalQ95,GammaMedian);
end

function assessment = assess_high_results(summary,reference,groups,cfg)
Group = groups(:);
RMSE_HuTGRS = zeros(numel(groups),1);
RMSE_EmpiricalTarget = zeros(numel(groups),1);
TailToEarlySlopeRatio = zeros(numel(groups),1);
ConcaveFraction = zeros(numel(groups),1);
for index = 1:numel(groups)
    selected = sortrows(summary(summary.Group == groups(index),:),'U10');
    y = selected.TotalMedian;
    slopes = diff(y)./diff(selected.U10);
    early = selected.U10(1:end-1) < 16;
    tail = selected.U10(1:end-1) >= cfg.tailSlopeStart;
    RMSE_HuTGRS(index) = sqrt(mean((y-reference.TgrsHuTotal).^2));
    RMSE_EmpiricalTarget(index) = sqrt(mean((y-reference.EmpiricalTotal).^2));
    TailToEarlySlopeRatio(index) = mean(slopes(tail))/mean(slopes(early));
    ConcaveFraction(index) = mean(diff(slopes) <= 0);
end

elfRMSE = sqrt(mean((reference.ElfouhailyTotal-reference.TgrsHuTotal).^2));
ElfouhailyIntegralRMSE_HuTGRS = repmat(elfRMSE,numel(groups),1);
assessment = table(Group,RMSE_HuTGRS,RMSE_EmpiricalTarget, ...
    ElfouhailyIntegralRMSE_HuTGRS,TailToEarlySlopeRatio,ConcaveFraction);

if cfg.requireNonlinearImprovement
    linear = assessment(assessment.Group == groups(1),:);
    corrected = assessment(assessment.Group == groups(3),:);
    correctedCurve = sortrows(summary(summary.Group == groups(3),:),'U10');
    assert(corrected.RMSE_HuTGRS < linear.RMSE_HuTGRS, ...
        'Corrected nonlinear MSS did not outperform the linear realization.');
    assert(corrected.RMSE_HuTGRS < elfRMSE, ...
        'Corrected nonlinear MSS did not outperform integrated Elfouhaily.');
    assert(all(diff(correctedCurve.TotalMedian) > 0), ...
        'Corrected high-wind MSS must be strictly increasing.');
    assert(corrected.TailToEarlySlopeRatio < cfg.maximumTailToEarlySlopeRatio, ...
        'The corrected curve does not flatten sufficiently at high wind.');
    assert(corrected.ConcaveFraction >= cfg.minimumConcaveFraction, ...
        'Too few corrected-curve intervals are concave.');
end
end

function fig = plot_high_wind_mss(summary,reference,groups,cfg)
% Publication-ready high-sea-state figure for IEEE two-column use.
% Panel (a): total MSS. Panel (b): finite-difference MSS growth rate.

style = ieee_style();
fig = figure('Visible',cfg.figureVisible,'Color','w','Units','inches', ...
    'Position',[0.6 0.6 7.16 3.15],'Renderer','painters');
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

linear = sortrows(summary(summary.Group == groups(1),:),'U10');
rawLie = sortrows(summary(summary.Group == groups(2),:),'U10');
proposed = sortrows(summary(summary.Group == groups(3),:),'U10');

% (a) High-wind total MSS
ax1 = nexttile(tl,1); hold(ax1,'on');
fill(ax1,[proposed.U10;flipud(proposed.U10)], ...
    [proposed.TotalQ05;flipud(proposed.TotalQ95)],style.proposed, ...
    'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
h1 = plot(ax1,linear.U10,linear.TotalMedian,'s-', ...
    'Color',style.linear,'MarkerFaceColor','w','MarkerSize',3.8, ...
    'LineWidth',1.15,'DisplayName','Elfouhaily et al.--Linear realization');
h2 = plot(ax1,rawLie.U10,rawLie.TotalMedian,'^-', ...
    'Color',style.raw,'MarkerFaceColor','w','MarkerSize',3.8, ...
    'LineWidth',1.10,'DisplayName','Proposed--Raw modified Lie');
h3 = plot(ax1,proposed.U10,proposed.TotalMedian,'o-', ...
    'Color',style.proposed,'MarkerFaceColor','w','MarkerSize',4.0, ...
    'LineWidth',1.50,'DisplayName','Proposed--MSS-constrained Lie');
h4 = plot(ax1,reference.U10,reference.CoxMunkTotal,'--', ...
    'Color',style.cox,'LineWidth',1.15,'DisplayName','Cox and Munk--Empirical MSS');
h5 = plot(ax1,reference.U10,reference.TgrsHuTotal,'-.', ...
    'Color',style.hu,'LineWidth',1.35,'DisplayName','Hu et al.--High-wind MSS');
h6 = plot(ax1,reference.U10,reference.ElfouhailyTotal,':', ...
    'Color',style.elf,'LineWidth',1.25,'DisplayName','Elfouhaily et al.--Unified spectrum');
xline(ax1,cfg.highWindTransition,':','Color',[0.68 0.68 0.68], ...
    'LineWidth',0.85,'HandleVisibility','off');
xlabel(ax1,'$U_{10}$ (m/s)','Interpreter','latex');
ylabel(ax1,'Total MSS','Interpreter','latex');
text(ax1,0.035,0.96,'(a)','Units','normalized','FontWeight','bold', ...
    'FontSize',8.4,'VerticalAlignment','top','BackgroundColor','w','Margin',0.5);
grid(ax1,'on'); box(ax1,'on'); xlim(ax1,[10 40]);

% (b) Growth rate: finite-difference derivative.
ax2 = nexttile(tl,2); hold(ax2,'on');
plot_growth(ax2,linear.U10,linear.TotalMedian,style.linear,'s-');
plot_growth(ax2,rawLie.U10,rawLie.TotalMedian,style.raw,'^-');
plot_growth(ax2,proposed.U10,proposed.TotalMedian,style.proposed,'o-');
plot_growth(ax2,reference.U10,reference.CoxMunkTotal,style.cox,'--');
plot_growth(ax2,reference.U10,reference.TgrsHuTotal,style.hu,'-.');
plot_growth(ax2,reference.U10,reference.ElfouhailyTotal,style.elf,':');
xline(ax2,cfg.highWindTransition,':','Color',[0.68 0.68 0.68], ...
    'LineWidth',0.85,'HandleVisibility','off');
yline(ax2,0,':','Color',[0.55 0.55 0.55],'HandleVisibility','off');
xlabel(ax2,'$U_{10}$ (m/s)','Interpreter','latex');
ylabel(ax2,'$\Delta$MSS/$\Delta U_{10}$','Interpreter','latex');
text(ax2,0.035,0.96,'(b)','Units','normalized','FontWeight','bold', ...
    'FontSize',8.4,'VerticalAlignment','top','BackgroundColor','w','Margin',0.5);
grid(ax2,'on'); box(ax2,'on'); xlim(ax2,[10 40]);

lgd = legend(ax1,[h1 h2 h3 h4 h5 h6], ...
    'Location','southoutside','Orientation','horizontal','NumColumns',3, ...
    'Box','off','FontSize',7.2,'Interpreter','latex');
lgd.Layout.Tile = 'south';
lgd.ItemTokenSize = [15 7];
apply_ieee_axes([ax1 ax2]);
end

function plot_growth(ax,U,y,color,lineSpec)
midU = 0.5*(U(1:end-1)+U(2:end));
slope = diff(y)./diff(U);
if contains(lineSpec,'o') || contains(lineSpec,'s') || contains(lineSpec,'^')
    plot(ax,midU,slope,lineSpec,'Color',color,'MarkerFaceColor','w', ...
        'MarkerSize',3.4,'LineWidth',1.10,'HandleVisibility','off');
else
    plot(ax,midU,slope,lineSpec,'Color',color,'LineWidth',1.15, ...
        'HandleVisibility','off');
end
end

function apply_ieee_axes(ax)
set(ax,'FontName','Times New Roman','FontSize',7.8,'LineWidth',0.65, ...
    'TickDir','out','TickLength',[0.016 0.016],'Layer','top', ...
    'GridAlpha',0.15,'MinorGridAlpha',0.08, ...
    'TickLabelInterpreter','latex');
end

function style = ieee_style()
style.linear = [0.35 0.35 0.35];
style.raw = [0.78 0.45 0.12];
style.proposed = [0.00 0.36 0.62];
style.cox = [0.05 0.05 0.05];
style.hu = [0.72 0.24 0.18];
style.elf = [0.22 0.52 0.30];
end

function export_ieee_figure(fig,baseName,resolution)
exportgraphics(fig,[baseName '.pdf'],'ContentType','vector', ...
    'BackgroundColor','white');
exportgraphics(fig,[baseName '.png'],'Resolution',resolution, ...
    'BackgroundColor','white');
end

function value = local_quantile(values,p)
position = 1+(numel(values)-1)*p;
lower = floor(position); upper = ceil(position); weight = position-lower;
value = values(lower)*(1-weight)+values(upper)*weight;
end

function validate_high_config(cfg)
assert(all(diff(cfg.windSpeeds) > 0) && min(cfg.windSpeeds) >= 10, ...
    'High-sea-state wind speeds must be strictly increasing and >=10 m/s.');
assert(max(cfg.windSpeeds) > 20, ...
    'Include winds above 20 m/s to test flattening.');
assert(strcmpi(cfg.empiricalTargetMode,'hu-high-wind'), ...
    'High-wind validation must use the Hu/TGRS logarithmic target mode.');
assert(cfg.mssResidualFraction >= 0 && cfg.mssResidualFraction < 0.2, ...
    'Use residualFraction in [0,0.2) for stable high-wind moment dressing.');
end
