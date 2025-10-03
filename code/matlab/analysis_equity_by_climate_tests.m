%% equity_by_climate_from_CityEquity.m  (updated to report NDVI Balance)
% Tests whether equity changes differ by climate (Biome / Köppen)
% Requires: CityEquity (table) and A_HYPO_MEANS (table with BIOMES_cat / KOPPEN_cat)

%% ---------- Guards & merge ----------
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
if any(lower(keep)=="biomes_cat"), Meta.BIOMES_cat = categorical(string(Meta.BIOMES_cat)); end
if any(lower(keep)=="koppen_cat"), Meta.KOPPEN_cat = categorical(string(Meta.KOPPEN_cat)); end

J = innerjoin(CityEquity, Meta, 'Keys','City');

%% ---------- Outcomes (match typology figure) ----------
% ΔLST (Bottom–Top; positive = equity closing for LST)
if ismember('LME_DeltaSlope_LST', J.Properties.VariableNames)
    J.DeltaLST = J.LME_DeltaSlope_LST;
else
    J.DeltaLST = J.Slope_BotLST - J.Slope_TopLST;
end

% NDVI Balance Index: + means "bottom-up" (bottom improving more in magnitude)
% If precomputed as B_NDVI, use that; otherwise derive from per-tail slopes.
if ismember('B_NDVI', J.Properties.VariableNames)
    J.NDVI_Balance = J.B_NDVI;
else
    % robust to missing columns; set to NaN if not present
    if all(ismember({'Slope_BotNDVI','Slope_TopNDVI'}, J.Properties.VariableNames))
        J.NDVI_Balance = abs(J.Slope_BotNDVI) - abs(J.Slope_TopNDVI);
    else
        J.NDVI_Balance = NaN(height(J),1);
    end
end

% ΔNDVI (Top–Bottom; negative = equity closing for NDVI)
if ismember('LME_DeltaSlope_NDVI', J.Properties.VariableNames)
    J.DeltaNDVI = J.LME_DeltaSlope_NDVI;
elseif all(ismember({'Slope_TopNDVI','Slope_BotNDVI'}, J.Properties.VariableNames))
    J.DeltaNDVI = J.Slope_TopNDVI - J.Slope_BotNDVI;
end

%% ---------- Helper: one-way tests with singleton drop ----------
function OUT = run_oneway(tbl, yname, gname)
    % --- pull and coerce
    y = tbl.(yname);
    g = tbl.(gname);
    if ~iscategorical(g), g = categorical(string(g)); end

    % --- drop non-finite y, then singleton groups
    ok = isfinite(y) & ~ismissing(g);
    y = y(ok); g = g(ok);
    [cats,~,ic] = unique(g,'stable');
    Nper = accumarray(ic,1);
    keepCats = cats(Nper >= 2);
    rows = ismember(g, keepCats);
    y = y(rows); g = removecats(g(rows));

    % --- scaffold
    OUT = struct('yname', yname, ...
                 'gname', gname, ...
                 'levels', { {categories(g)} }, ...
                 'p_anova', NaN, ...
                 'eta2', NaN, ...
                 'p_kw', NaN, ...
                 'tukey', table(), ...
                 'groups', table());

    if numel(categories(g)) < 2 || numel(y) < 3
        return
    end

    % --- ANOVA (η² via anovan), fallback to anova1
    try
        [pA, AT] = anovan(y, {g}, 'display','off');   % one factor
        SSg = AT{2,2};  SSe = AT{3,2};
        OUT.eta2    = SSg / (SSg + SSe);
        OUT.p_anova = pA;
    catch
        OUT.p_anova = anova1(y, g, 'off');
        OUT.eta2    = NaN;
    end

    % --- Kruskal–Wallis
    OUT.p_kw = kruskalwallis(y, g, 'off');

    % --- Group stats
    [G, levs] = findgroups(g);
    N  = splitapply(@numel, y, G);
    MU = splitapply(@(x) mean(x,'omitnan'), y, G);
    SD = splitapply(@(x)  std(x,'omitnan'), y, G);
    SEM = SD ./ sqrt(N);
    OUT.groups = table(levs, N, MU, SD, SEM, ...
        'VariableNames', {'Group','N','Mean','Std','SEM'});

    % --- Tukey HSD if ANOVA significant
    if isfinite(OUT.p_anova) && OUT.p_anova < 0.05
        try
            [~,~,stats] = anova1(y, g, 'off');
            tk = multcompare(stats, 'ctype','tukey-kramer', 'display','off');
            T = array2table(tk, 'VariableNames', ...
                {'g1','g2','LowerCI','Diff','UpperCI','pAdj'});
            L = categories(g);
            T.g1 = categorical(L(T.g1));
            T.g2 = categorical(L(T.g2));
            OUT.tukey = T;
        catch
        end
    end
end

function print_block(titleStr, OUT)
    fprintf('\n=== %s ===\n', titleStr);
    if ~isempty(OUT.groups) && height(OUT.groups)>0
        disp(OUT.groups);
    else
        fprintf('(No groups to compare — after filtering, <2 levels or <3 rows)\n');
    end
    if ~isnan(OUT.p_anova)
        fprintf('ANOVA p = %.4g, η² = %.3f | KW p = %.4g\n', OUT.p_anova, OUT.eta2, OUT.p_kw);
    else
        fprintf('KW p = %.4g (ANOVA unavailable)\n', OUT.p_kw);
    end
    if ~isempty(OUT.tukey) && height(OUT.tukey)>0
        try
            nSig = sum(OUT.tukey.pAdj < 0.05);
            fprintf('Tukey (α=0.05): %d significant pair(s).\n', nSig);
        catch
        end
    end
end

%% ---------- Run tests ----------
Results = struct();

if ismember('BIOMES_cat', J.Properties.VariableNames)
    Results.BIOME = struct();
    Results.BIOME.DeltaLST      = run_oneway(J, 'DeltaLST',     'BIOMES_cat');
    Results.BIOME.NDVI_Balance  = run_oneway(J, 'NDVI_Balance', 'BIOMES_cat');
    if ismember('DeltaNDVI', J.Properties.VariableNames)
        Results.BIOME.DeltaNDVI = run_oneway(J, 'DeltaNDVI', 'BIOMES_cat');
    end
    print_block('Tests by BIOME (ΔLST)',        Results.BIOME.DeltaLST);
    print_block('Tests by BIOME (NDVI Balance)',Results.BIOME.NDVI_Balance);
    if isfield(Results.BIOME,'DeltaNDVI')
        print_block('Tests by BIOME (ΔNDVI)',   Results.BIOME.DeltaNDVI);
    end
end

if ismember('KOPPEN_cat', J.Properties.VariableNames)
    Results.KOPPEN = struct();
    Results.KOPPEN.DeltaLST      = run_oneway(J, 'DeltaLST',     'KOPPEN_cat');
    Results.KOPPEN.NDVI_Balance  = run_oneway(J, 'NDVI_Balance', 'KOPPEN_cat');
    if ismember('DeltaNDVI', J.Properties.VariableNames)
        Results.KOPPEN.DeltaNDVI = run_oneway(J, 'DeltaNDVI', 'KOPPEN_cat');
    end
    print_block('Tests by KÖPPEN (ΔLST)',        Results.KOPPEN.DeltaLST);
    print_block('Tests by KÖPPEN (NDVI Balance)',Results.KOPPEN.NDVI_Balance);
    if isfield(Results.KOPPEN,'DeltaNDVI')
        print_block('Tests by KÖPPEN (ΔNDVI)',   Results.KOPPEN.DeltaNDVI);
    end
end

assignin('base','EquityClimateResults', Results);
fprintf('\nSaved -> EquityClimateResults (group means, ANOVA/KW, Tukey if ANOVA sig.).\n');
