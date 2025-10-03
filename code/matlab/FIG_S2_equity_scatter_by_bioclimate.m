%% figureS1_typology_scatter.m
% ------------------------------------------------------------------------------
% SI Figure S1: Typology scatter — ΔNDVI (equity change) vs ΔLST, colored by Biome/Köppen
%
% Inputs (preferred, from equity_build_and_analyze.m):
%   T_NDVI : table with City, Slope_Bottom, Slope_Top, (optionally) Delta_TopMinusBottom
%   T_LST  : table with City, Slope_Bottom, Slope_Top, Delta_TopMinusBottom
%   A_HYPO_MEANS : table with City and a typology column (BIOMES_cat or KOPPEN_cat)
% If T_NDVI/T_LST are not in the workspace, the script will try to load
% them from equity_panel.mat.
% ------------------------------------------------------------------------------

%% ===================== USER CONFIG =====================
pointSize      = 55;      % scatter marker size
alphaFace      = 0.85;    % face alpha
alphaEdge      = 0.25;    % edge alpha
labelExtremesN = 3;       % how many extremes to auto-label on each axis

% --- output toggles ---
WritePNG = false;          % <- set false to skip PNG
WritePDF = false;          % <- set false to skip PDF
pngOut   = 'SI_typology_scatter_DeltaNDVI_vs_DeltaLST.png';
pdfOut   = 'SI_typology_scatter_DeltaNDVI_vs_DeltaLST.pdf';

% --- typology toggle: 'BIOME' or 'KOPPEN' ---
TypologyMode = 'BIOME';   % options: 'BIOME' (BIOMES_cat) | 'KOPPEN' (KOPPEN_cat)
%% ======================================================

% ---- Pull LME per-city slopes (preferred) ----
if ~(exist('T_NDVI','var')==1 && istable(T_NDVI)) || ~(exist('T_LST','var')==1 && istable(T_LST))
    assert(exist('equity_panel.mat','file')==2, ...
        'T_NDVI/T_LST not in workspace. Run the build script or load equity_panel.mat.');
    Sload = load('equity_panel.mat','T_NDVI','T_LST');
    T_NDVI = Sload.T_NDVI;
    T_LST  = Sload.T_LST;
end

% ---- Choose typology column on A_HYPO_MEANS ----
assert(exist('A_HYPO_MEANS','var')==1 && istable(A_HYPO_MEANS), ...
       'A_HYPO_MEANS not found. Load the predictors with BIOMES_cat/KOPPEN_cat.');

v = string(A_HYPO_MEANS.Properties.VariableNames);
switch upper(strtrim(TypologyMode))
    case 'BIOME'
        wantCol = "BIOMES_cat";
        errMsg  = 'A_HYPO_MEANS must contain BIOMES_cat for TypologyMode=''BIOME''.';
    case 'KOPPEN'
        wantCol = "KOPPEN_cat";
        errMsg  = 'A_HYPO_MEANS must contain KOPPEN_cat for TypologyMode=''KOPPEN''.';
    otherwise
        error('TypologyMode must be ''BIOME'' or ''KOPPEN''.');
end
ixTyp = find(v == wantCol, 1);
assert(~isempty(ixTyp), errMsg);

% ---- column guards on LME tables ----
needN_base = {'City','Slope_Bottom','Slope_Top'};
needL      = {'City','Slope_Bottom','Slope_Top','Delta_TopMinusBottom'};
assert(all(ismember(needN_base, T_NDVI.Properties.VariableNames)), 'T_NDVI missing required columns (City, slopes).');
assert(all(ismember(needL,      T_LST .Properties.VariableNames)), 'T_LST missing required columns.');

% ---- compute features on copies (don’t mutate originals) ----
tN = T_NDVI(:, needN_base);
tL = T_LST (  :, {'City','Slope_Bottom','Slope_Top','Delta_TopMinusBottom'});

% Ensure robust string joins
tN.City = string(tN.City);
tL.City = string(tL.City);

% Gather typology from A_HYPO_MEANS and normalize the name to 'Typology'
meta = A_HYPO_MEANS(:, ["City", wantCol]);
meta.Properties.VariableNames = {'City','Typology'};
meta.City = string(meta.City);

% ---- X-axis: ΔNDVI (Top − Bottom) per decade ----
% Prefer precomputed Delta_TopMinusBottom if present; otherwise compute from slopes
if ismember('Delta_TopMinusBottom', T_NDVI.Properties.VariableNames)
    tN.DeltaNDVI = T_NDVI.Delta_TopMinusBottom;
else
    tN.DeltaNDVI = tN.Slope_Top - tN.Slope_Bottom;
end

% ---- Y-axis: ΔLST (Top − Bottom) per decade (from LME table) ----
tL.DeltaLST = tL.Delta_TopMinusBottom;

% Merge NDVI + LST + typology
JL = innerjoin(innerjoin(tN, tL, 'Keys','City'), meta, 'Keys','City');

% Filter valid rows and well-defined typologies
grp = categorical(string(JL.Typology));
ok  = isfinite(JL.DeltaNDVI) & isfinite(JL.DeltaLST) & ~isundefined(grp);
JL  = JL(ok,:);
grp = removecats(categorical(string(JL.Typology)));

% (soft check) ΔLST should be >= 0 in current panel; warn if not
if any(JL.DeltaLST < -1e-8)
    warning('Some ΔLST values are negative; verify model sign conventions.');
end

% ---- color palette by typology ----
cats = categories(grp);
K    = numel(cats);
palette = lines(max(K,4));                     % distinct, colorblind-friendly baseline
Cmap = containers.Map(cats, num2cell(palette(1:K,:),2));

% ---- plot ----
figure('Color','w','Position',[100 100 900 650]); hold on; box on; grid on;

% scatter by typology
for k = 1:K
    idx = grp == cats{k};
    scatter(JL.DeltaNDVI(idx), JL.DeltaLST(idx), pointSize, ...
        'MarkerFaceColor', Cmap(cats{k}), 'MarkerEdgeColor','k', ...
        'MarkerFaceAlpha', alphaFace, 'MarkerEdgeAlpha', alphaEdge);
end

% reference lines at zero
plot([0 0], ylim, 'k:','LineWidth',1);
plot(xlim, [0 0], 'k:','LineWidth',1);

% labels, legend, title
xlabel('\DeltaNDVI = slope_{Top} − slope_{Bottom}  (per decade)');  % more negative = equity gain
ylabel('\DeltaLST = slope_{Top} − slope_{Bottom}  (°C per decade)');
title(sprintf('Typology of equity pathways by %s', ...
    ternary(strcmpi(TypologyMode,'KOPPEN'),'Köppen','biome')));

lgd = legend(cats, 'Location', 'bestoutside');
try
    lgd.Title.String = ternary(strcmpi(TypologyMode,'KOPPEN'),'Köppen','Biome');
catch
    try
        title(lgd, ternary(strcmpi(TypologyMode,'KOPPEN'),'Köppen','Biome'));
    catch
        pos = lgd.Position; % normalized figure units
        annotation('textbox', [pos(1) pos(2)+pos(4)  pos(3) 0.03], ...
            'String', ternary(strcmpi(TypologyMode,'KOPPEN'),'Köppen','Biome'), ...
            'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold');
    end
end

% quadrant annotations (reworded for ΔNDVI rather than balance)
xl = xlim; yl = ylim;
dx = 0.03*(xl(2)-xl(1)); dy = 0.04*(yl(2)-yl(1));
text(xl(2)-dx, yl(2)-dy, 'Top-tail warming \uparrow', 'HorizontalAlignment','right','Color',[0.25 0.25 0.25]);
text(xl(2)-dx, yl(1)+dy, 'Bottom warms faster \downarrow', 'HorizontalAlignment','right','Color',[0.25 0.25 0.25]);
text(xl(1)+dx, yl(2)-3*dy, '\leftarrow NDVI equity gain (gap closes)', 'HorizontalAlignment','left','Color',[0.1 0.35 0.9]);
text(xl(2)-dx, yl(2)-3*dy, 'NDVI equity loss (gap widens) \rightarrow', 'HorizontalAlignment','right','Color',[0.85 0.2 0.2]);

% optional: label a few extremes for readability
[~,ix_hi] = maxk(abs(JL.DeltaNDVI), labelExtremesN);
[~,iy_hi] = maxk(abs(JL.DeltaLST),  labelExtremesN);
labIdx = unique([ix_hi; iy_hi]);
for ii = labIdx(:)'
    text(JL.DeltaNDVI(ii), JL.DeltaLST(ii), " " + string(JL.City(ii)), ...
        'FontSize',8, 'Color',[0.2 0.2 0.2], 'HorizontalAlignment','left', ...
        'VerticalAlignment','middle');
end

% ---- final axis housekeeping ----
% y-axis must start at 0; choose a modest headroom for the top
ymax = max(JL.DeltaLST(:));
if ~isfinite(ymax) || isempty(ymax), ymax = 1; end
pad  = 0.05 * max(1e-6, ymax);         % small 5% headroom
ylim([0, ymax + pad]);

% keep x-lims auto, but nudge a little pad
xl = xlim;
xpad = 0.03 * max(1e-6, (xl(2)-xl(1)));
xlim([xl(1)-xpad, xl(2)+xpad]);

set(gca,'Layer','top','TickDir','out');

% ---- Save to disk (toggled) ----
if WritePNG
    exportgraphics(gcf, pngOut, 'Resolution', 300);
end
if WritePDF
    try, exportgraphics(gcf, pdfOut); catch, end
end
if WritePNG || WritePDF
    fprintf('Saved figure to:\n');
    if WritePNG, fprintf('  %s\n', pngOut); end
    if WritePDF, fprintf('  %s\n', pdfOut); end
else
    fprintf('Save toggles are OFF (WritePNG=false, WritePDF=false). No files written.\n');
end

% -------------- Local tiny helper --------------
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
