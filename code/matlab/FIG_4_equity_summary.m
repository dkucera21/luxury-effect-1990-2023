%% figure4_equity_summary.m
% ------------------------------------------------------------------------------
% 4-panel equity summary figure
%  Top row  : interaction plots with 95% CIs (NDVI, LST) from stacked-tail LMEs
%  Bottom   : dumbbells of Top-N ΔLST and ΔNDVI (NDVI spines colored by balance)
% Requirements in workspace (or in equity_panel.mat):
%   - T_NDVI (table with City, Slope_Bottom, Slope_Top, Delta_TopMinusBottom)
%   - T_LST  (table with City, Slope_Bottom, Slope_Top, Delta_TopMinusBottom)
%   - lmeN   (LinearMixedModel for NDVI stacked tails)
%   - lmeL   (LinearMixedModel for LST  stacked tails)
% ------------------------------------------------------------------------------

%% ===================== USER CONFIG =====================
nShow     = 15;                     % show top-N cities in dumbbells
clrBottom = [0.16 0.45 0.90];       % low-income marker color (blue)
clrTop    = [0.95 0.49 0.20];       % high-income marker color (orange)

SaveFig   = true;                   % write PNG/PDF to disk?
OutBase   = 'equity_panel_figure';  % base filename (no extension)
DPI       = 300;                    % PNG resolution
%% =======================================================

%% ---------- Try to pull required vars; fall back to MAT file ----------
needLoad = ~(exist('T_NDVI','var')==1 && istable(T_NDVI) && ...
             exist('T_LST','var')==1  && istable(T_LST)  && ...
             exist('lmeN','var')==1   && isa(lmeN,'LinearMixedModel') && ...
             exist('lmeL','var')==1   && isa(lmeL,'LinearMixedModel'));

if needLoad && exist('equity_panel.mat','file')==2
    S = load('equity_panel.mat');
    if ~exist('T_NDVI','var') && isfield(S,'T_NDVI'), T_NDVI = S.T_NDVI; end
    if ~exist('T_LST','var')  && isfield(S,'T_LST'),  T_LST  = S.T_LST;  end
    if ~exist('lmeN','var')   && isfield(S,'lmeN'),   lmeN   = S.lmeN;   end
    if ~exist('lmeL','var')   && isfield(S,'lmeL'),   lmeL   = S.lmeL;   end
end

assert(exist('T_NDVI','var')==1 && istable(T_NDVI), 'T_NDVI not found.');
assert(exist('T_LST','var')==1  && istable(T_LST),  'T_LST not found.');
assert(exist('lmeN','var')==1   && isa(lmeN,'LinearMixedModel'), 'lmeN not found.');
assert(exist('lmeL','var')==1   && isa(lmeL,'LinearMixedModel'), 'lmeL not found.');

%% ---------- Prep dumbbell tables ----------
TopL = sortrows(T_LST,  'Delta_TopMinusBottom', 'descend');
TopL = TopL(1:min(nShow, height(TopL)), :);

if ~ismember('Balance', T_NDVI.Properties.VariableNames)
    T_NDVI.Balance = abs(T_NDVI.Slope_Bottom) - abs(T_NDVI.Slope_Top); % + = bottom-up, - = top-down
end
TopN = sortrows(T_NDVI, 'Delta_TopMinusBottom', 'ascend');             % most negative Δ = biggest NDVI equity gain
TopN = TopN(1:min(nShow, height(TopN)), :);

B0 = max(abs(TopN.Balance)); if isempty(B0) || B0==0, B0 = 0.05; end

%% ---------- Layout ----------
f = figure('Color','w','Position',[100 80 1400 900]);
tlo = tiledlayout(f,2,2,'TileSpacing','compact','Padding','loose');

% ===================== A) NDVI interaction (top-left) =====================
axA = nexttile(tlo,1); hold(axA,'on');
plot_tail_panel(axA, lmeN, 'NDVI: Top vs Bottom tails', 'NDVI (unitless)');

% ===================== B) LST interaction (top-right) =====================
axB = nexttile(tlo,2); hold(axB,'on');
plot_tail_panel(axB, lmeL, 'LST: Top vs Bottom tails', 'LST (°C)');

% ===================== C) ΔLST dumbbells (bottom-left) ====================
axC = nexttile(tlo,4); hold(axC,'on'); box(axC,'off');
y = (1:height(TopL))';
for i = 1:height(TopL)
    b = TopL.Slope_Bottom(i);
    t = TopL.Slope_Top(i);

    % grey spine
    plot(axC, [b t], [y(i) y(i)], '-', 'Color', 0.65*[1 1 1], 'LineWidth', 4);

    % endpoints
    plot(axC, b, y(i), 'o', 'MarkerFaceColor', clrBottom, 'MarkerEdgeColor','none', 'MarkerSize', 7);
    plot(axC, t, y(i), 'o', 'MarkerFaceColor', clrTop,    'MarkerEdgeColor','none', 'MarkerSize', 7);

    % annotate delta
    d = TopL.Delta_TopMinusBottom(i);
    xpad = 0.04*range([TopL.Slope_Bottom;TopL.Slope_Top]);
    text(axC, max(b,t)+xpad, y(i), sprintf('\\Delta = %+0.2f', d), ...
        'FontSize', 10, 'Color', [0.30 0.30 0.30], 'VerticalAlignment','middle');
end
set(axC,'YTick',y,'YTickLabel',TopL.City,'YDir','reverse','TickDir','out');
xlabel(axC,'LST slope (°C / decade)'); ylabel(axC,'City');
title(axC,'C) Top 10 \DeltaLST equity gains');
% legend
% legend (use explicit handles so the grey spine isn't used)
hB = plot(axC, NaN, NaN, 'o', 'MarkerFaceColor', clrBottom, 'MarkerEdgeColor','none','MarkerSize',7);
hT = plot(axC, NaN, NaN, 'o', 'MarkerFaceColor', clrTop,    'MarkerEdgeColor','none','MarkerSize',7);
legend(axC, [hB hT], {'Bottom (low-income)','Top (high-income)'}, ...
       'Location','southoutside', 'Box','off');
% x-lims with small pad
xlC = [min([TopL.Slope_Bottom;TopL.Slope_Top]) max([TopL.Slope_Bottom;TopL.Slope_Top])];
pad = 0.05*range(xlC); if pad==0, pad=0.1; end
xlim(axC, xlC + [-pad pad]);

% ===================== D) ΔNDVI dumbbells (bottom-right) ==================
axD = nexttile(tlo,3); hold(axD,'on'); box(axD,'off');
colormap(axD, balanceCmap_posBlue(256));  % diverging colormap for balance
caxis(axD, [-B0 B0]);

y2 = (1:height(TopN))';
for i = 1:height(TopN)
    b = TopN.Slope_Bottom(i);
    t = TopN.Slope_Top(i);
    bal = TopN.Balance(i);  % |bottom| - |top|

    % colored spine by balance index
    spineColor = colorFromBalance(axD, bal, B0);
    plot(axD, [b t], [y2(i) y2(i)], '-', 'Color', spineColor, 'LineWidth', 4);

    % endpoints
    plot(axD, b, y2(i), 'o', 'MarkerFaceColor', clrBottom, 'MarkerEdgeColor','none', 'MarkerSize', 7);
    plot(axD, t, y2(i), 'o', 'MarkerFaceColor', clrTop,    'MarkerEdgeColor','none', 'MarkerSize', 7);

    % annotate Δ and balance
    d = TopN.Delta_TopMinusBottom(i);
    xpad = 0.04*range([TopN.Slope_Bottom;TopN.Slope_Top]);
	text(axD, max(b,t)+xpad, y2(i), sprintf('\\Delta = %+0.4f', d), ...
        'FontSize', 10, 'Color', [0.25 0.25 0.25], 'VerticalAlignment','middle');
end
set(axD,'YTick',y2,'YTickLabel',TopN.City,'YDir','reverse','TickDir','out');
xlabel(axD,'NDVI slope (per decade)'); ylabel(axD,'City');
title(axD,'D) Top 10 \DeltaNDVI equity gains');

% add NDVI balance colorbar UNDER the NDVI panel (no overlap)
cb = colorbar(axD, 'Location','southoutside');
cb.Label.String = 'Balance index (|bottom| − |top|)';
cb.Limits = [-B0 B0];
cb.Ticks  = [-B0 0 B0];
cb.TickLabels = compose('%+.2f', cb.Ticks);

% tighten x-lims with a pad
xlD = [min([TopN.Slope_Bottom;TopN.Slope_Top]) max([TopN.Slope_Bottom;TopN.Slope_Top])];
padD = 0.05*range(xlD); if padD==0, padD=0.02; end
xlim(axD, xlD + [-padD padD]);

%% ---------- Numeric sanity (fixed-effects only lines) ----------
yrs = [1990 2023]'; time = (yrs - 1990)/10;
mk = @(tail) table(yrs, time, categorical(repmat(tail,numel(yrs),1),["Bottom","Top"]), ...
                   categorical(repmat("AnyCity",numel(yrs),1)), ...
                   'VariableNames', {'Year','time_dec','tail','City'});

[yB_N,~] = predict(lmeN, mk("Bottom"), 'Conditional', false);
[yT_N,~] = predict(lmeN, mk("Top"),    'Conditional', false);
[yB_L,~] = predict(lmeL, mk("Bottom"), 'Conditional', false);
[yT_L,~] = predict(lmeL, mk("Top"),    'Conditional', false);

fprintf('\nNDVI FE lines (1990→2023): Bottom %.4f→%.4f (Δ=%+0.4f) | Top %.4f→%.4f (Δ=%+0.4f) | Gap Δ (Top-Bot) = %+0.4f\n', ...
    yB_N(1), yB_N(2), yB_N(2)-yB_N(1), yT_N(1), yT_N(2), yT_N(2)-yT_N(1), (yT_N(2)-yB_N(2)) - (yT_N(1)-yB_N(1)));
fprintf('LST  FE lines (1990→2023): Bottom %.4f→%.4f (Δ=%+0.4f) | Top %.4f→%.4f (Δ=%+0.4f) | Gap Δ (Bot-Top) = %+0.4f\n\n', ...
    yB_L(1), yB_L(2), yB_L(2)-yB_L(1), yT_L(1), yT_L(2), yT_L(2)-yT_L(1), (yB_L(2)-yT_L(2)) - (yB_L(1)-yT_L(1)));

%% ---------- Optional save ----------
if SaveFig
    exportgraphics(f, [OutBase '.png'], 'Resolution', DPI);
    try, exportgraphics(f, [OutBase '.pdf']); catch, end
    fprintf('Figure written: %s.{png,pdf}\n', OutBase);
end

%% ======================================================================
%% ------------------------- Helper functions ---------------------------
%% ======================================================================
function plot_tail_panel(ax, lme, panelTitle, yLabel)
    years   = (1990:2023)'; timeDec = (years - 1990)/10;

    anyCity = lme.Variables.City(find(~ismissing(lme.Variables.City),1));
    mkNew = @(tl) table(years, timeDec, categorical(repmat(string(tl),numel(years),1),["Bottom","Top"]), ...
                        repmat(anyCity,numel(years),1), 'VariableNames',{'Year','time_dec','tail','City'});
    newB = mkNew('Bottom'); newT = mkNew('Top');

    [yB, ciB] = predict(lme, newB, 'Conditional', false, 'Alpha', 0.05);
    [yT, ciT] = predict(lme, newT, 'Conditional', false, 'Alpha', 0.05);

    CF   = lme.Coefficients; cn = string(CF.Name);
    idxI = find(contains(cn,"time_dec:tail"), 1);
    bInt = CF.Estimate(idxI); pInt = CF.pValue(idxI);

    fillCI(ax, years, ciB(:,1), ciB(:,2), [0 0.45 0.74], 0.15);
    fillCI(ax, years, ciT(:,1), ciT(:,2), [0.85 0.33 0.10], 0.15);
    p1 = plot(ax, years, yB, 'LineWidth', 2.5, 'Color',[0 0.45 0.74]);
    p2 = plot(ax, years, yT, 'LineWidth', 2.5, 'Color',[0.85 0.33 0.10],'LineStyle','--');

    grid(ax,'on'); box(ax,'off');
    xlabel(ax,'Year'); ylabel(ax,yLabel); title(ax,panelTitle,'FontWeight','bold');
    legend(ax,[p1 p2],{'Bottom (low-income)','Top (high-income)'},'Location','best','Box','off');

    yl = ylim(ax);
    txt = sprintf('time×tail = %.3f (p = %s)', bInt, ternary(pInt<1e-3,'<0.001',num2str(pInt,'%.3f')));
    text(ax, years(1)+0.5, yl(2)-0.05*(yl(2)-yl(1)), txt, 'VerticalAlignment','top', ...
         'FontSize',10,'BackgroundColor',[1 1 1 0.6],'Margin',6);
end

function fillCI(axh, x, lo, hi, baseColor, alphaVal)
    xv = [x; flipud(x)]; yv = [lo; flipud(hi)];
    fc = baseColor + (1-baseColor)*0.6;
    patch('Parent',axh,'XData',xv,'YData',yv,'FaceColor',fc, ...
          'EdgeColor','none','FaceAlpha',alphaVal);
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function C = balanceCmap_posBlue(n)
% Diverging map: negative -> red, zero -> white, positive -> blue
if nargin<1, n = 256; end
half = floor(n/2);
red  = [0.84 0.15 0.16];
blue = [0.16 0.45 0.90];

% left half (negatives): red -> white
C1 = [linspace(red(1),1,half)' linspace(red(2),1,half)' linspace(red(3),1,half)'];
% right half (positives): white -> blue
C2 = [linspace(1,blue(1),n-half)' linspace(1,blue(2),n-half)' linspace(1,blue(3),n-half)'];
C  = [C1; C2];
end

function col = colorFromBalance(ax, bal, B0)
    C = colormap(ax); n = size(C,1);
    if B0 <= 0, col = [0.5 0.5 0.5]; return; end
    u = max(-B0, min(B0, bal));
    t = (u + B0) / (2*B0);
    idx = max(1, min(n, 1 + round(t*(n-1))));
    col = C(idx,:);
end
