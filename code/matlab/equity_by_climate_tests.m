%% equity_by_climate_from_CityEquity.m  (fixed)
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
if ismember('LME_DeltaSlope_LST', J.Properties.VariableNames)
    J.DeltaLST = J.LME_DeltaSlope_LST;
else
    J.DeltaLST = J.Slope_BotLST - J.Slope_TopLST;
end
if ismember('B_NDVI', J.Properties.VariableNames)
    J.NDVI_Balance = J.B_NDVI;
else
    J.NDVI_Balance = abs(J.Slope_BotNDVI) - abs(J.Slope_TopNDVI);
end
if ismember('LME_DeltaSlope_NDVI', J.Properties.VariableNames)
    J.DeltaNDVI = J.LME_DeltaSlope_NDVI;
end

ok = isfinite(J.DeltaLST) & isfinite(J.NDVI_Balance);
J  = J(ok,:);

%% ---------- Helper: one-way tests with singleton drop ----------
function OUT = run_oneway(tbl, yname, gname)
    % --- pull and coerce
    y = tbl.(yname);
    g = tbl.(gname);
    if ~iscategorical(g), g = categorical(string(g)); end

    % --- drop singleton groups (e.g., your single Köppen class)
    [cats,~,ic] = unique(g,'stable');
    Nper = accumarray(ic,1);
    keepCats = cats(Nper >= 2);
    rows = ismember(g, keepCats);
    y = y(rows);
    g = removecats(g(rows));

    % --- scalar struct; wrap array fields in a cell
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

    % --- ANOVA (η² via anovan for robust SS), fallback to anova1
    try
        [pA, AT] = anovan(y, {g}, 'display','off');   % one factor
        % AT rows: {Source, SS, df, MS, F, Prob>F}
        SSg = AT{2,2};  SSe = AT{3,2};
        OUT.eta2    = SSg / (SSg + SSe);
        OUT.p_anova = pA;
    catch
        OUT.p_anova = anova1(y, g, 'off');
        OUT.eta2    = NaN;
    end

    % --- Kruskal–Wallis
    OUT.p_kw = kruskalwallis(y, g, 'off');

    % --- Group stats (mean, std, N, SEM) via splitapply (no groupsummary)
    [G, levs] = findgroups(g);
    N  = splitapply(@numel, y, G);
    MU = splitapply(@(x) mean(x,'omitnan'), y, G);
    SD = splitapply(@(x)  std(x,'omitnan'), y, G);
    SEM = SD ./ sqrt(N);

    GT = table(levs, N, MU, SD, SEM, ...
        'VariableNames', {'Group','N','Mean','Std','SEM'});
    OUT.groups = GT;

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
            % ignore Tukey errors silently
        end
    end
end



%% ---------- Run tests ----------
Results = struct();

if ismember('BIOMES_cat', J.Properties.VariableNames)
    fprintf('\n=== Tests by BIOME (singletons removed) ===\n');
    Results.BIOME.DeltaLST      = run_oneway(J, 'DeltaLST',     'BIOMES_cat');
    Results.BIOME.NDVI_Balance  = run_oneway(J, 'NDVI_Balance', 'BIOMES_cat');
    if ismember('DeltaNDVI', J.Properties.VariableNames)
        Results.BIOME.DeltaNDVI = run_oneway(J, 'DeltaNDVI', 'BIOMES_cat');
    end
    disp(Results.BIOME.DeltaLST.groups);
    fprintf('ANOVA p (ΔLST~Biome)=%.4g, η²=%.3f | KW p=%.4g\n', ...
        Results.BIOME.DeltaLST.p_anova, Results.BIOME.DeltaLST.eta2, Results.BIOME.DeltaLST.p_kw);
end

if ismember('KOPPEN_cat', J.Properties.VariableNames)
    fprintf('\n=== Tests by KÖPPEN (singletons removed) ===\n');
    Results.KOPPEN.DeltaLST      = run_oneway(J, 'DeltaLST',     'KOPPEN_cat');
    Results.KOPPEN.NDVI_Balance  = run_oneway(J, 'NDVI_Balance', 'KOPPEN_cat');
    if ismember('DeltaNDVI', J.Properties.VariableNames)
        Results.KOPPEN.DeltaNDVI = run_oneway(J, 'DeltaNDVI', 'KOPPEN_cat');
    end
    disp(Results.KOPPEN.DeltaLST.groups);
    fprintf('ANOVA p (ΔLST~Köppen)=%.4g, η²=%.3f | KW p=%.4g\n', ...
        Results.KOPPEN.DeltaLST.p_anova, Results.KOPPEN.DeltaLST.eta2, Results.KOPPEN.DeltaLST.p_kw);
end

assignin('base','EquityClimateResults', Results);
fprintf('\nSaved -> EquityClimateResults (group means, ANOVA/KW, Tukey for significant ANOVAs).\n');
