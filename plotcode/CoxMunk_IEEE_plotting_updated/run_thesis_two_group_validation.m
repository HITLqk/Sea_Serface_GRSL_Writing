function [raw,summary,reference,assessment] = run_thesis_two_group_validation(cfg)
%RUN_THESIS_TWO_GROUP_VALIDATION Clean Elfouhaily versus modified-Lie test.
%   The only groups are Linear Elfouhaily and MSS-constrained Modified Lie.

arguments
    cfg (1,1) struct = default_thesis_two_group_config()
end
validate_config(cfg);
if ~isfolder(cfg.outputDirectory)
    mkdir(cfg.outputDirectory);
end

reference = two_group_mss_references(cfg.windSpeeds,cfg);
groups = ["Linear Elfouhaily","MSS-Constrained Modified Lie"];
nRows = numel(cfg.windSpeeds)*numel(cfg.realizationSeeds)*numel(groups);
U10 = zeros(nRows,1);
Seed = zeros(nRows,1);
Group = strings(nRows,1);
MssAlong = zeros(nRows,1);
MssCross = zeros(nRows,1);
MssTotal = zeros(nRows,1);
Gamma = zeros(nRows,1);
PrimaryMss = zeros(nRows,1);
ShortWaveMss = zeros(nRows,1);
row = 0;
example = struct();

for windIndex = 1:numel(cfg.windSpeeds)
    windSpeed = cfg.windSpeeds(windIndex);
    for seedIndex = 1:numel(cfg.realizationSeeds)
        seed = cfg.realizationSeeds(seedIndex);
        result = synthesize_two_group_realization(windSpeed,seed,cfg);
        row = add_row(row,windSpeed,seed,groups(1),result.linear, ...
            result.primaryLinear.total,result.shortWave.total);
        row = add_row(row,windSpeed,seed,groups(2),result.breaking, ...
            result.primaryBreaking.total,result.shortWaveBreaking.total);
        if seedIndex == 1 && (windSpeed == 5 || windSpeed == 10)
            example.(sprintf('U%d',windSpeed)) = result;
        end
    end
end

raw = table(U10,Seed,Group,MssAlong,MssCross,MssTotal,Gamma, ...
    PrimaryMss,ShortWaveMss);
summary = summarize_results(raw,groups,cfg.windSpeeds);
assessment = assess_results(summary,reference,groups);
if cfg.requireNonlinearImprovement
    verify_nonlinear_improvement(assessment,groups);
end

figMss = plot_mss(summary,reference,groups,cfg);
export_ieee_figure(figMss,fullfile(cfg.outputDirectory, ...
    'Fig_low_moderate_mss_validation'),cfg.exportResolution);
close(figMss);

writetable(raw,fullfile(cfg.outputDirectory,'two_group_raw.csv'));
writetable(summary,fullfile(cfg.outputDirectory,'two_group_summary.csv'));
writetable(reference,fullfile(cfg.outputDirectory,'two_group_reference.csv'));
writetable(assessment,fullfile(cfg.outputDirectory,'two_group_assessment.csv'));
save(fullfile(cfg.outputDirectory,'two_group_validation.mat'), ...
    'raw','summary','reference','assessment','cfg','example','-v7.3');
disp(assessment);

    function nextRow = add_row(currentRow,wind,seedValue,groupName,mss,primary,shortWave)
        nextRow = currentRow+1;
        U10(nextRow) = wind;
        Seed(nextRow) = seedValue;
        Group(nextRow) = groupName;
        MssAlong(nextRow) = mss.along;
        MssCross(nextRow) = mss.cross;
        MssTotal(nextRow) = mss.total;
        Gamma(nextRow) = mss.gamma;
        PrimaryMss(nextRow) = primary;
        ShortWaveMss(nextRow) = shortWave;
    end
end

function summary = summarize_results(raw,groups,winds)
n = numel(groups)*numel(winds);
Group = strings(n,1); U10 = zeros(n,1);
AlongMedian = zeros(n,1); AlongQ05 = zeros(n,1); AlongQ95 = zeros(n,1);
CrossMedian = zeros(n,1); CrossQ05 = zeros(n,1); CrossQ95 = zeros(n,1);
TotalMedian = zeros(n,1); TotalQ05 = zeros(n,1); TotalQ25 = zeros(n,1);
TotalQ75 = zeros(n,1); TotalQ95 = zeros(n,1); GammaMedian = zeros(n,1);
row = 0;
for groupIndex = 1:numel(groups)
    for windIndex = 1:numel(winds)
        row = row+1;
        selected = raw.Group == groups(groupIndex) & raw.U10 == winds(windIndex);
        Group(row) = groups(groupIndex); U10(row) = winds(windIndex);
        qAlong = local_quantile(raw.MssAlong(selected),[0.05 0.5 0.95]);
        qCross = local_quantile(raw.MssCross(selected),[0.05 0.5 0.95]);
        AlongQ05(row)=qAlong(1); AlongMedian(row)=qAlong(2); AlongQ95(row)=qAlong(3);
        CrossQ05(row)=qCross(1); CrossMedian(row)=qCross(2); CrossQ95(row)=qCross(3);
        q = local_quantile(raw.MssTotal(selected),[0.05 0.25 0.5 0.75 0.95]);
        TotalQ05(row)=q(1); TotalQ25(row)=q(2); TotalMedian(row)=q(3);
        TotalQ75(row)=q(4); TotalQ95(row)=q(5);
        GammaMedian(row) = median(raw.Gamma(selected));
    end
end
summary = table(Group,U10,AlongMedian,AlongQ05,AlongQ95, ...
    CrossMedian,CrossQ05,CrossQ95,TotalMedian,TotalQ05, ...
    TotalQ25,TotalQ75,TotalQ95,GammaMedian);
end

function assessment = assess_results(summary,reference,groups)
Group = groups(:);
RMSE_CoxMunk = zeros(numel(groups),1);
RMSE_Guerin = zeros(numel(groups),1);
RMSE_TGRS_Hu = zeros(numel(groups),1);
RMSE_Elfouhaily = zeros(numel(groups),1);
RMSE_EmpiricalFusion = zeros(numel(groups),1);
MeanAbsGammaError_TGRS = zeros(numel(groups),1);
for index = 1:numel(groups)
    selected = summary.Group == groups(index);
    groupSummary = sortrows(summary(selected,:),'U10');
    total = groupSummary.TotalMedian;
    RMSE_CoxMunk(index) = sqrt(mean((total-reference.CoxMunkTotal).^2));
    valid = isfinite(reference.GuerinTotal);
    RMSE_Guerin(index) = sqrt(mean((total(valid)-reference.GuerinTotal(valid)).^2));
    RMSE_TGRS_Hu(index) = sqrt(mean((total-reference.TgrsHuTotal).^2));
    RMSE_Elfouhaily(index) = sqrt(mean((total-reference.ElfouhailyTotal).^2));
    RMSE_EmpiricalFusion(index) = sqrt(mean((total-reference.EmpiricalTotal).^2));
    MeanAbsGammaError_TGRS(index) = mean(abs( ...
        groupSummary.GammaMedian-reference.TgrsGamma));
end
assessment = table(Group,RMSE_CoxMunk,RMSE_Guerin,RMSE_TGRS_Hu, ...
    RMSE_Elfouhaily,RMSE_EmpiricalFusion,MeanAbsGammaError_TGRS);
end

function verify_nonlinear_improvement(assessment,groups)
linear = assessment(assessment.Group == groups(1),:);
nonlinear = assessment(assessment.Group == groups(2),:);
metrics = {'RMSE_CoxMunk','RMSE_Guerin','RMSE_TGRS_Hu','RMSE_EmpiricalFusion'};
for index = 1:numel(metrics)
    assert(nonlinear.(metrics{index}) < linear.(metrics{index}), ...
        'Nonlinear model did not improve %s. Inspect the dressing target.',metrics{index});
end
end

function fig = plot_mss(summary,reference,groups,cfg)
% Publication-ready low/moderate-wind MSS figure for IEEE two-column use.
% Panels: total, along-wind and crosswind MSS.  Linear and proposed curves
% share the same axes; literature references are visually secondary.

style = ieee_style();
fig = figure('Visible',cfg.figureVisible,'Color','w','Units','inches', ...
    'Position',[0.6 0.6 7.16 3.45],'Renderer','painters');
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

linear = sortrows(summary(summary.Group == groups(1),:),'U10');
proposed = sortrows(summary(summary.Group == groups(2),:),'U10');

% (a) Total MSS
ax1 = nexttile(tl,1); hold(ax1,'on');
plot_band(ax1,proposed.U10,proposed.TotalQ05,proposed.TotalQ95,style.proposed);
h1 = plot(ax1,linear.U10,linear.TotalMedian,'s-', ...
    'Color',style.linear,'MarkerFaceColor','w','MarkerSize',3.6, ...
    'LineWidth',1.15,'DisplayName','Elfouhaily et al.--Linear realization');
h2 = plot(ax1,proposed.U10,proposed.TotalMedian,'o-', ...
    'Color',style.proposed,'MarkerFaceColor','w','MarkerSize',3.8, ...
    'LineWidth',1.45,'DisplayName','Proposed--MSS-constrained Lie');
h3 = plot(ax1,reference.U10,reference.CoxMunkTotal,'--', ...
    'Color',style.cox,'LineWidth',1.15,'DisplayName','Cox and Munk--Empirical MSS');
h4 = plot(ax1,reference.U10,reference.TgrsHuTotal,'-.', ...
    'Color',style.hu,'LineWidth',1.15,'DisplayName','Hu et al.--MSS law');
h5 = plot(ax1,reference.U10,reference.GuerinTotal,'d', ...
    'Color',style.guerin,'MarkerFaceColor','w','MarkerSize',3.8, ...
    'LineStyle','none','LineWidth',0.9,'DisplayName','Guerin et al.--IASI MSS');
h6 = plot(ax1,reference.U10,reference.ElfouhailyTotal,':', ...
    'Color',style.elf,'LineWidth',1.25,'DisplayName','Elfouhaily et al.--Unified spectrum');
h7 = plot(ax1,reference.U10,reference.EmpiricalTotal,'--', ...
    'Color',style.constraint,'LineWidth',0.95,'DisplayName','Proposed--Fusion constraint');
format_mss_axis(ax1,'Total MSS','(a)');

% (b) Along-wind MSS
ax2 = nexttile(tl,2); hold(ax2,'on');
plot_band(ax2,proposed.U10,proposed.AlongQ05,proposed.AlongQ95,style.proposed);
plot(ax2,linear.U10,linear.AlongMedian,'s-', ...
    'Color',style.linear,'MarkerFaceColor','w','MarkerSize',3.6,'LineWidth',1.15);
plot(ax2,proposed.U10,proposed.AlongMedian,'o-', ...
    'Color',style.proposed,'MarkerFaceColor','w','MarkerSize',3.8,'LineWidth',1.45);
plot(ax2,reference.U10,reference.CoxMunkAlong,'--','Color',style.cox,'LineWidth',1.15);
plot(ax2,reference.U10,reference.TgrsHuAlong,'-.','Color',style.hu,'LineWidth',1.15);
plot(ax2,reference.U10,reference.GuerinAlong,'d','Color',style.guerin, ...
    'MarkerFaceColor','w','MarkerSize',3.8,'LineStyle','none','LineWidth',0.9);
plot(ax2,reference.U10,reference.ElfouhailyAlong,':','Color',style.elf,'LineWidth',1.25);
plot(ax2,reference.U10,reference.EmpiricalAlong,'--','Color',style.constraint,'LineWidth',0.95);
format_mss_axis(ax2,'Along-wind MSS','(b)');

% (c) Crosswind MSS
ax3 = nexttile(tl,3); hold(ax3,'on');
plot_band(ax3,proposed.U10,proposed.CrossQ05,proposed.CrossQ95,style.proposed);
plot(ax3,linear.U10,linear.CrossMedian,'s-', ...
    'Color',style.linear,'MarkerFaceColor','w','MarkerSize',3.6,'LineWidth',1.15);
plot(ax3,proposed.U10,proposed.CrossMedian,'o-', ...
    'Color',style.proposed,'MarkerFaceColor','w','MarkerSize',3.8,'LineWidth',1.45);
plot(ax3,reference.U10,reference.CoxMunkCross,'--','Color',style.cox,'LineWidth',1.15);
plot(ax3,reference.U10,reference.TgrsHuCross,'-.','Color',style.hu,'LineWidth',1.15);
plot(ax3,reference.U10,reference.GuerinCross,'d','Color',style.guerin, ...
    'MarkerFaceColor','w','MarkerSize',3.8,'LineStyle','none','LineWidth',0.9);
plot(ax3,reference.U10,reference.ElfouhailyCross,':','Color',style.elf,'LineWidth',1.25);
plot(ax3,reference.U10,reference.EmpiricalCross,'--','Color',style.constraint,'LineWidth',0.95);
format_mss_axis(ax3,'Crosswind MSS','(c)');

% Shared limits improve visual comparison without forcing identical y ranges.
set([ax1 ax2 ax3],'XLim',[1 10],'XTick',1:2:10);

% One shared, centered legend. The fusion target is explicitly identified as
% a constraint rather than an independent validation reference.
lgd = legend(ax1,[h1 h2 h3 h4 h5 h6 h7], ...
    'Location','southoutside','Orientation','horizontal','NumColumns',4, ...
    'Box','off','FontSize',7.0,'Interpreter','latex');
lgd.Layout.Tile = 'south';
lgd.ItemTokenSize = [14 7];
apply_ieee_axes([ax1 ax2 ax3]);
end

function plot_band(ax,x,q05,q95,color)
fill(ax,[x;flipud(x)],[q05;flipud(q95)],color, ...
    'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
end

function format_mss_axis(ax,yText,panelLabel)
xlabel(ax,'$U_{10}$ (m/s)','Interpreter','latex');
ylabel(ax,yText,'Interpreter','latex');
grid(ax,'on'); box(ax,'on');
text(ax,0.035,0.96,panelLabel,'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top', ...
    'FontWeight','bold','FontSize',8.2,'BackgroundColor','w','Margin',0.5);
end

function apply_ieee_axes(ax)
set(ax,'FontName','Times New Roman','FontSize',7.7,'LineWidth',0.65, ...
    'TickDir','out','TickLength',[0.018 0.018],'Layer','top', ...
    'GridAlpha',0.15,'MinorGridAlpha',0.08, ...
    'TickLabelInterpreter','latex');
end

function style = ieee_style()
style.linear = [0.35 0.35 0.35];
style.proposed = [0.00 0.36 0.62];
style.cox = [0.05 0.05 0.05];
style.hu = [0.72 0.24 0.18];
style.guerin = [0.48 0.30 0.62];
style.elf = [0.22 0.52 0.30];
style.constraint = [0.62 0.62 0.62];
end

function export_ieee_figure(fig,baseName,resolution)
exportgraphics(fig,[baseName '.pdf'],'ContentType','vector', ...
    'BackgroundColor','white');
exportgraphics(fig,[baseName '.png'],'Resolution',resolution, ...
    'BackgroundColor','white');
end

function fig = plot_gamma(summary,reference,groups,cfg)
fig = figure('Visible',cfg.figureVisible,'Color','w', ...
    'Position',[120 120 1260 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
colors = lines(numel(groups));
for index = 1:numel(groups)
    ax = nexttile; hold(ax,'on');
    selected = sortrows(summary(summary.Group == groups(index),:),'U10');
    plot(ax,selected.U10,selected.GammaMedian,'o-', ...
        'Color',colors(index,:),'LineWidth',1.8,'DisplayName','Simulation');
    plot(ax,reference.U10,reference.ElfouhailyGamma,'k:','LineWidth',1.6, ...
        'DisplayName','Elfouhaily integral');
    plot(ax,reference.U10,reference.GuerinGamma,'s','Color',[0.1 0.35 0.8], ...
        'MarkerFaceColor','none','DisplayName','Guerin IASI');
    yline(ax,0.864,'-.','TGRS mean 0.864','Color',[0.8 0.2 0.1], ...
        'LineWidth',1.5,'LabelHorizontalAlignment','left', ...
        'HandleVisibility','off');
    yline(ax,0.84,'--','TGRS simulated 0.84','Color',[0.45 0.2 0.6], ...
        'LineWidth',1.2,'LabelHorizontalAlignment','left', ...
        'HandleVisibility','off');
    title(ax,groups(index)); xlabel(ax,'U_{10} (m/s)');
    ylabel(ax,'Anisotropy \gamma'); grid(ax,'on'); box(ax,'on'); xlim(ax,[1 10]);
    if index == 1, legend(ax,'Location','best'); end
end
linkaxes(findall(fig,'Type','axes'),'xy');
end

function fig = plot_examples(example,cfg)
fig = figure('Visible',cfg.figureVisible,'Color','w', ...
    'Position',[140 80 1480 800]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
winds = [5 10];
for windIndex = 1:numel(winds)
    result = example.(sprintf('U%d',winds(windIndex)));
    surfaces = {result.linearSurface,result.breakingSurface};
    names = {'Linear Elfouhaily','MSS-Constrained Modified Lie'};
    commonLimit = max(abs([surfaces{1}(:);surfaces{2}(:)]));
    coordinates = (0:size(surfaces{1},1)-1)*result.primarySpacing;
    for groupIndex = 1:2
        ax = nexttile; imagesc(ax,coordinates,coordinates,surfaces{groupIndex});
        axis(ax,'image'); axis(ax,'xy'); caxis(ax,[-commonLimit commonLimit]);
        colormap(ax,parula); colorbar(ax);
        title(ax,sprintf('%s, U_{10}=%d m/s',names{groupIndex},winds(windIndex)));
        xlabel(ax,'x (m)'); ylabel(ax,'y (m)');
    end
    ax = nexttile; hold(ax,'on');
    centerRow = floor(size(surfaces{1},1)/2)+1;
    plot(ax,coordinates,surfaces{1}(centerRow,:),'LineWidth',1.4, ...
        'DisplayName','Linear');
    plot(ax,coordinates,surfaces{2}(centerRow,:),'LineWidth',1.4, ...
        'DisplayName','Modified Lie');
    title(ax,sprintf('Paired center section, U_{10}=%d m/s',winds(windIndex)));
    xlabel(ax,'x (m)'); ylabel(ax,'Elevation (m)'); grid(ax,'on'); box(ax,'on');
    legend(ax,'Location','best');
end
end

function q = local_quantile(values,p)
values = sort(values(:));
positions = 1+(numel(values)-1)*p;
lower = floor(positions); upper = ceil(positions);
weight = positions-lower;
q = values(lower).*(1-weight)+values(upper).*weight;
end

function validate_config(cfg)
assert(isequal(cfg.windSpeeds(:),(1:10)'), ...
    'This experiment is fixed to the requested 1:1:10 m/s wind grid.');
assert(all(cfg.realizationSeeds == floor(cfg.realizationSeeds)), ...
    'Realization seeds must be integers.');
assert(cfg.primaryPeakSamples >= 10, ...
    'The primary grid must satisfy the thesis dk <= kp/10 condition.');
assert(cfg.lieOutputPeakMultiple <= cfg.primaryMaximumPeakMultiple, ...
    'The Lie output band must be contained in the primary grid band.');
assert(cfg.maximumOpticalWavenumber > 370, ...
    'The optical cutoff must include the Elfouhaily short-wave peak.');
assert(cfg.mssResidualFraction >= 0 && cfg.mssResidualFraction < 1, ...
    'mssResidualFraction must be in [0,1).');
end
