%% figure_climate_equity_boxes.m
% Boxplots (with jitter) of equity-change metrics by climate.
% Produces 2 figures: ΔLST and ΔNDVI, each with BIOME & Köppen panels.

%% ---------- Guards & join (reuse your existing data) ----------
assert(exist('CityEquity','var')==1 && istable(CityEquity), 'CityEquity not found.');
assert(exist('A_HYPO_MEANS','var')==1 && istable(A_HYPO_MEANS), 'A_HYPO_MEANS not found.');

CityEquity.City   = string(CityEquity.City);
A_HYPO_MEANS.City = string(A_HYPO_MEANS.City);

v = string(A_HYPO_MEANS.Properties.VariableNames);
ixBiome  = find(lower(v)=="biomes_cat", 1);
ixKoppen = find(contains(lower(v),"koppen") & endsWith(lower(v),"_cat"), 1);
if isempty(ixBiome) && isempty(ixKoppen)
    error('A_HYPO_MEANS must contain BIOMES_cat and/or KOPPEN_cat.');
end
keep = "City"; if ~isempty(ixBiome), keep(end+1)=v(ixBiome); end
if ~isempty(ixKoppen), keep(end+1)=v(ixKoppen); end
Meta = A_HYPO_MEANS(:, keep);
if ismember("BIOMES_cat", string(Meta.Properties.VariableNames)) && ~iscategorical(Meta.BIOMES_cat)
    Meta.BIOMES_cat = categorical(string(Meta.BIOMES_cat));
end
if ismember("KOPPEN_cat", string(Meta.Properties.VariableNames)) && ~iscategorical(Meta.KOPPEN_cat)
    Meta.KOPPEN_cat = categorical(string(Meta.KOPPEN_cat));
end

J = innerjoin(CityEquity, Meta, 'Keys','City');

% Outcomes (parallel definitions)
% ΔLST = Bottom − Top (°C/decade) used in your panel
if ismember('LME_DeltaSlope_LST', J.Properties.VariableNames)
    J.DeltaLST = J.LME_DeltaSlope_LST;
else
    J.DeltaLST = J.Slope_BotLST - J.Slope_TopLST;
end

% ΔNDVI = Top − Bottom (NDVI/decade) to mirror LME convention in typology code
if ismember('LME_DeltaSlope_NDVI', J.Properties.VariableNames)
    J.DeltaNDVI = J.LME_DeltaSlope_NDVI;
else
    J.DeltaNDVI = J.Slope_TopNDVI - J.Slope_BotNDVI;
end

% Valid rows only
ok = isfinite(J.DeltaLST) & isfinite(J.DeltaNDVI);
J  = J(ok,:);

%% ---------- Helpers (no groupsummary anywhere) ----------
function T = group_stats_no_gs(g, y)
    if ~iscategorical(g), g = categorical(string(g)); end
    [cats,~,ic] = unique(g,'stable');
    N  = accumarray(ic, 1);
    mu = accumarray(ic, y, [], @(x) mean(x,'omitnan'));
    sd = accumarray(ic, y, [], @(x) std(x, 'omitnan'));
    sem = sd ./ sqrt(max(N,1));
    T = table(cats, N, mu, sd, sem, 'VariableNames', {'Group','N','Mean','Std','SEM'});
end

function OUT = oneway_with_tukey(tbl, yname, gname)
    y = tbl.(yname);
    g = tbl.(gname);
    if ~iscategorical(g), g = categorical(string(g)); end

    % Drop singleton groups
    [cats,~,ic] = unique(g,'stable'); N = accumarray(ic,1);
    keepCats = cats(N>=2);
    rows = ismember(g, keepCats);
    y = y(rows); g = removecats(g(rows));

    OUT.groups = table(); OUT.tukey = table(); OUT.letters = table();
    OUT.g = g; OUT.y = y;
    if numel(categories(g)) < 2 || numel(y) < 3, return; end

    % Group stats (manual)
    OUT.groups = group_stats_no_gs(g, y);

    % Tukey letters (if ANOVA works)
    try
        [~,~,stats] = anova1(y, g, 'off');
        tk = multcompare(stats, 'ctype','tukey-kramer', 'display','off');
        T = array2table(tk, 'VariableNames', {'g1','g2','LowerCI','Diff','UpperCI','pAdj'});
        levs = categories(g);
        T.g1 = categorical(levs(T.g1));
        T.g2 = categorical(levs(T.g2));
        OUT.tukey = T;
        OUT.letters = compact_letters(T, categories(g), 0.05);
    catch
        % leave tukey/letters empty
    end
end

function CL = compact_letters(T, levelNames, alpha)
    L = string(levelNames(:));
    K = numel(L);
    sig = false(K,K);
    for r = 1:height(T)
        i = find(L == string(T.g1(r)), 1);
        j = find(L == string(T.g2(r)), 1);
        if ~isempty(i) && ~isempty(j) && T.pAdj(r) < alpha
            sig(i,j) = true; sig(j,i) = true;
        end
    end
    letters = strings(K,1);
    used = strings(0,1);
    for i = 1:K
        placed = false;
        for li = 1:numel(used)
            letter = used(li);
            ok = true;
            for j = 1:K
                if letters(j)==letter && sig(i,j)
                    ok = false; break;
                end
            end
            if ok
                letters(i) = letter; placed = true; break;
            end
        end
        if ~placed
            newL = char('A' + numel(used));
            used(end+1) = string(newL);
            letters(i)  = string(newL);
        end
    end
    CL = table(categorical(L), letters, 'VariableNames', {'Group','Letter'});
end

function plot_one(ax, tbl, yname, gname, titleStr, ylab)
    S = oneway_with_tukey(tbl, yname, gname);
    if isempty(S.groups) || height(S.groups)==0
        cla(ax); title(ax, [titleStr ' (insufficient groups)']); return
    end

    % Order groups by mean
    G = S.groups;
    [~,ord] = sort(G.Mean, 'descend');
    levs = string(G.Group(ord));
    g = reordercats(categorical(S.g, levs), levs);
    y = S.y;

    axes(ax); cla(ax); hold on; box on; grid on;

    % Box layer (fallback to boxplot if boxchart is unavailable)
    if exist('boxchart','file') == 2
        boxchart(g, y, 'BoxFaceColor',[0.8 0.85 0.95], 'WhiskerLineColor',[0.3 0.3 0.3]);
    else
        xnum = double(g);
        boxplot(y, xnum, 'Colors','k', 'Symbol','k.', 'OutlierSize',3);
        set(findobj(gca,'Tag','Box'),'LineWidth',1.2);
        set(gca,'XTick',1:numel(levs),'XTickLabel',cellstr(levs));
    end

    % Jittered points
    ug = categories(g);
    for ii = 1:numel(ug)
        idx = g==ug{ii};
        x0  = ii + 0.10*randn(nnz(idx),1);
        scatter(x0, y(idx), 18, 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0, ...
            'MarkerFaceColor',[0.2 0.4 0.7]);
    end

    % Mean ± SEM bars
    for ii = 1:height(G)
        gi = find(strcmp(ug, string(G.Group(ii))));
        if isempty(gi), continue; end
        m = G.Mean(ii); se = G.SEM(ii);
        plot([gi-0.18 gi+0.18], [m m], 'k-', 'LineWidth',1.5);
        if isfinite(se) && se>0
            line([gi gi],[m-se m+se], 'Color','k','LineWidth',1);
            line(gi+[-0.08 0.08], [m-se m-se], 'Color','k','LineWidth',1);
            line(gi+[-0.08 0.08], [m+se m+se], 'Color','k','LineWidth',1);
        end
    end

    % Tukey letters above boxes (if available)
    if ~isempty(S.letters)
        Ltab = S.letters;
        ylims = ylim; yr = ylims(2)-ylims(1);
        pad = 0.03*yr;
        for ii = 1:numel(ug)
            hit = find(string(Ltab.Group)==string(ug{ii}), 1);
            if ~isempty(hit)
                txt = char(Ltab.Letter(hit));
                ymax = max(y(g==ug{ii}));
                text(ii, ymax + pad, txt, 'HorizontalAlignment','center', 'FontWeight','bold', 'Color',[0.15 0.15 0.15]);
            end
        end
    end

    ylabel(ylab); title(titleStr);
    set(ax,'TickDir','out','Layer','top');
end

%% ---------- Build figures ----------
% Figure 1: ΔLST by BIOME and Köppen
f1 = figure('Color','w','Position',[80 80 1200 500]);
t1 = tiledlayout(f1,1,2,'TileSpacing','compact','Padding','loose');

ax1 = nexttile(t1,1);
if ismember('BIOMES_cat', J.Properties.VariableNames)
    plot_one(ax1, J, 'DeltaLST', 'BIOMES_cat', '\DeltaLST by Biome', '\DeltaLST = slope_{Top}-slope_{Bottom} (°C/decade)');
else
    title(ax1,'BIOME not available'); axis(ax1,'off');
end

ax2 = nexttile(t1,2);
if ismember('KOPPEN_cat', J.Properties.VariableNames)
    plot_one(ax2, J, 'DeltaLST', 'KOPPEN_cat', '\DeltaLST by Köppen', '\DeltaLST = slope_{Top}-slope_{Bottom} (°C/decade)');
else
    title(ax2,'Köppen not available'); axis(ax2,'off');
end

exportgraphics(f1, 'Fig_Climate_DeltaLST_Boxplots.png', 'Resolution', 300);
try, exportgraphics(f1, 'Fig_Climate_DeltaLST_Boxplots.pdf'); end

% Figure 2: ΔNDVI by BIOME and Köppen (parity with ΔLST)
f2 = figure('Color','w','Position',[80 80 1200 500]);
t2 = tiledlayout(f2,1,2,'TileSpacing','compact','Padding','loose');

ax3 = nexttile(t2,1);
if ismember('BIOMES_cat', J.Properties.VariableNames)
    plot_one(ax3, J, 'DeltaNDVI', 'BIOMES_cat', '\DeltaNDVI by Biome', '\DeltaNDVI = slope_{Top}-slope_{Bottom} (NDVI/decade)');
else
    title(ax3,'BIOME not available'); axis(ax3,'off');
end

ax4 = nexttile(t2,2);
if ismember('KOPPEN_cat', J.Properties.VariableNames)
    plot_one(ax4, J, 'DeltaNDVI', 'KOPPEN_cat', '\DeltaNDVI by Köppen', '\DeltaNDVI = slope_{Top}-slope_{Bottom} (NDVI/decade)');
else
    title(ax4,'Köppen not available'); axis(ax4,'off');
end

exportgraphics(f2, 'Fig_Climate_DeltaNDVI_Boxplots.png', 'Resolution', 300);
try, exportgraphics(f2, 'Fig_Climate_DeltaNDVI_Boxplots.pdf'); end

fprintf('Saved:\n  Fig_Climate_DeltaLST_Boxplots.{png,pdf}\n  Fig_Climate_DeltaNDVI_Boxplots.{png,pdf}\n');
