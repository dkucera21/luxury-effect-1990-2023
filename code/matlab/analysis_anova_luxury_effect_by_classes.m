%% anova_luxury_effect_by_classes.m  (SCRIPT)
% -------------------------------------------------------------------------
% Test whether city-level luxury-effect slopes differ by Köppen climate and
% by Biome. Robust to multiple column name variants.
%
% Requires in workspace: A_HYPO_MEANS (table)
% Acceptable column name aliases:
%   NDVI luxury slope:  "LUX_NDVI_per10k","NDVI_LUX_per10k",
%                       "NDVI_LUX_TREND","lux_ndvi_per10k",
%                       "lux_ndvi_slope","NDVI_LUX_SLOPE"
%   LST  luxury slope:  "LUX_LST_per10k","LST_LUX_per10k",
%                       "LST_LUX_TREND","lux_lst_per10k",
%                       "lux_lst_slope","LST_LUX_SLOPE"
%   Köppen category:    "KoppenClimate_cat","KOPPEN_cat","koppen","koppen_cat"
%   Biome category:     "BIOMES_cat","Biome","biome","biomes_cat"
% -------------------------------------------------------------------------

%% ------------------ Load & normalize ------------------
assert(exist('A_HYPO_MEANS','var')==1 && istable(A_HYPO_MEANS), ...
    'A_HYPO_MEANS not found. Build it first.');

T = A_HYPO_MEANS;
v = string(T.Properties.VariableNames);

% Helper: first matching column name from a list of aliases (case-insensitive)
find_col = @(aliases) v( find(ismember(lower(v), lower(string(aliases))), 1, 'first') );

% Aliases
luxN_aliases = ["LUX_NDVI_per10k","NDVI_LUX_per10k","NDVI_LUX_TREND","lux_ndvi_per10k","lux_ndvi_slope","NDVI_LUX_SLOPE"];
luxL_aliases = ["LUX_LST_per10k", "LST_LUX_per10k", "LST_LUX_TREND", "lux_lst_per10k", "lux_lst_slope", "LST_LUX_SLOPE"];
kop_aliases  = ["KoppenClimate_cat","KOPPEN_cat","koppen","koppen_cat"];
bio_aliases  = ["BIOMES_cat","Biome","biome","biomes_cat"];

% Resolve names present in the table
cLuxN = find_col(luxN_aliases);
cLuxL = find_col(luxL_aliases);
cKop  = find_col(kop_aliases);
cBio  = find_col(bio_aliases);

% Helpful assert with diagnostics
if isempty(cLuxN) || isempty(cLuxL) || isempty(cKop) || isempty(cBio)
    msg = compose("Columns present: %s", strjoin(v, ', '));
    need = [
        "Need one NDVI luxury column among: " + strjoin(luxN_aliases, ", ")
        "Need one LST  luxury column among: " + strjoin(luxL_aliases, ", ")
        "Need one Köppen column among:      " + strjoin(kop_aliases,  ", ")
        "Need one Biome column among:       " + strjoin(bio_aliases,  ", ")
    ];
    error("Required columns missing.\n%s\n\n%s", strjoin(need,newline), msg);
end

% Pull series
ndvi   = T.(cLuxN);
lst    = T.(cLuxL);
kopAll = T.(cKop);
bioAll = T.(cBio);

% Ensure categoricals
if ~iscategorical(kopAll), kopAll = categorical(string(kopAll)); end
if ~iscategorical(bioAll), bioAll = categorical(string(bioAll)); end

% Label strings for printing (use the actual names found)
labelLuxN = char(cLuxN);
labelLuxL = char(cLuxL);

% Tiny wrapper for consistent outputs
run_oneway = @(y,g) struct('anova', local_anova1(y,g), 'kw', local_kw(y,g));

%% ------------------ NDVI slope ------------------
ok0 = isfinite(ndvi) & ~isundefined(kopAll);
res_ndvi_koppen_all = run_oneway(ndvi(ok0), removecats(kopAll(ok0)));

% Keep Köppen classes with n >= 2 for one-way & two-way
kop0     = removecats(kopAll(ok0));
cats     = categories(kop0); counts = countcats(kop0);
keepCats = cats(counts >= 2);
keepK    = ismember(kop0, categorical(keepCats));
ndviK    = ndvi(ok0); ndviK = ndviK(keepK);
kopK     = removecats(kop0(keepK));

% Biome one-way (NDVI)
okB = isfinite(ndvi) & ~isundefined(bioAll);
res_ndvi_biome = run_oneway(ndvi(okB), removecats(bioAll(okB)));

% Two-way (Köppen × Biome) on filtered rows
bio0 = removecats(bioAll(ok0));     % subset by ok0
bioK = removecats(bio0(keepK));     % then subset by keepK
if numel(categories(kopK)) >= 2 && numel(categories(bioK)) >= 2
    p_ndvi_tw = anovan(ndviK, {kopK, bioK}, 'model','interaction', ...
        'varnames',{'Koppen','Biome'}, 'display','off');
else
    p_ndvi_tw = NaN;
end

%% ------------------ LST slope ------------------
ok0L = isfinite(lst) & ~isundefined(kopAll);
res_lst_koppen_all = run_oneway(lst(ok0L), removecats(kopAll(ok0L)));

kop0L     = removecats(kopAll(ok0L));
catsL     = categories(kop0L); countsL = countcats(kop0L);
keepCatsL = catsL(countsL >= 2);
keepKL    = ismember(kop0L, categorical(keepCatsL));
lstK      = lst(ok0L); lstK = lstK(keepKL);
kopKL     = removecats(kop0L(keepKL));

% Biome one-way (LST)
okBL = isfinite(lst) & ~isundefined(bioAll);
res_lst_biome = run_oneway(lst(okBL), removecats(bioAll(okBL)));

bio0L = removecats(bioAll(ok0L));
bioKL = removecats(bio0L(keepKL));
if numel(categories(kopKL)) >= 2 && numel(categories(bioKL)) >= 2
    p_lst_tw = anovan(lstK, {kopKL, bioKL}, 'model','interaction', ...
        'varnames',{'Koppen','Biome'}, 'display','off');
else
    p_lst_tw = NaN;
end

%% ------------------ Print compact summaries ------------------
fprintf('\n=== %s ~ Köppen (one-way; ORIGINAL) ===\n', labelLuxN);      local_print(res_ndvi_koppen_all);
fprintf('=== %s ~ Köppen (one-way; FILTERED n>=2) ===\n', labelLuxN); local_print(res_oneway_safe(res_ndvi_koppen_all, ndviK, kopK));
fprintf('=== %s ~ Biome  (one-way) ===\n',                           labelLuxN); local_print(res_ndvi_biome);
fprintf('NDVI slope two-way ANOVA (Köppen, Biome, Interaction): %s\n', mat2str(p_ndvi_tw,3));

fprintf('\n=== %s ~ Köppen (one-way; ORIGINAL) ===\n', labelLuxL);      local_print(res_lst_koppen_all);
fprintf('=== %s ~ Köppen (one-way; FILTERED n>=2) ===\n', labelLuxL); local_print(res_oneway_safe(res_lst_koppen_all, lstK, kopKL));
fprintf('=== %s ~ Biome  (one-way) ===\n',                           labelLuxL); local_print(res_lst_biome);
fprintf('LST slope two-way ANOVA (Köppen, Biome, Interaction): %s\n',  mat2str(p_lst_tw,3));

dropped = setdiff(cats, keepCats);
if ~isempty(dropped)
    fprintf('\nDropped Köppen classes (n=1): %s\n', strjoin(dropped, ', '));
end

%% ------------------ Tukey (on slopes) ------------------
Tukey_NDVI_Biome  = tukey_table(ndvi, bioAll);
Tukey_NDVI_Koppen = tukey_table(ndvi, kopAll);
Tukey_LST_Biome   = tukey_table(lst,  bioAll);
Tukey_LST_Koppen  = tukey_table(lst,  kopAll);

fprintf('\n=== Tukey %s ~ Biome (significant pairs) ===\n', labelLuxN);
disp(Tukey_NDVI_Biome(Tukey_NDVI_Biome.pAdj<0.05, :))
fprintf('\n=== Tukey %s ~ Köppen (significant pairs) ===\n', labelLuxN);
disp(Tukey_NDVI_Koppen(Tukey_NDVI_Koppen.pAdj<0.05, :))
fprintf('\n=== Tukey %s ~ Biome (significant pairs) ===\n', labelLuxL);
disp(Tukey_LST_Biome(Tukey_LST_Biome.pAdj<0.05, :))
fprintf('\n=== Tukey %s ~ Köppen (significant pairs) ===\n', labelLuxL);
disp(Tukey_LST_Koppen(Tukey_LST_Koppen.pAdj<0.05, :))

% Optional: group means for context
grpMean = @(y,g) groupsummary(table(y,g), 'g', 'mean','y');
fprintf('\nGroup means (%s by Biome; singletons only for context):\n', labelLuxN); disp(grpMean(ndvi, bioAll))
fprintf('\nGroup means (%s by Biome; singletons only for context):\n', labelLuxL); disp(grpMean(lst,  bioAll))
fprintf('\nGroup means (%s by Köppen; singletons only for context):\n', labelLuxN); disp(grpMean(ndvi, kopAll))
fprintf('\nGroup means (%s by Köppen; singletons only for context):\n', labelLuxL); disp(grpMean(lst,  kopAll))

%% ===================== Local helpers =====================
function out = res_oneway_safe(~, y, g)
    if nargin<2 || isempty(y) || isempty(g) || numel(unique(g))<2
        out = struct('anova',struct('p',NaN,'eta2',NaN,'tbl',[],'tukey',[],'k',0,'n',0), ...
                     'kw',   struct('p',NaN,'tbl',[],'stats',[]));
    else
        out = struct('anova', local_anova1(y,g), 'kw', local_kw(y,g));
    end
end

function local_print(res)
    % Pretty-print a one-way result struct from run_oneway()
    if ~isstruct(res) || ~isfield(res,'anova') || ~isfield(res,'kw')
        fprintf('[no result]\n'); return
    end
    a = res.anova; k = res.kw;
    if isempty(a) || ~isfield(a,'p')
        fprintf('[no ANOVA result]\n');
    else
        pA  = a.p;   if isempty(pA) || ~isfinite(pA), pAstr = 'NaN'; else, pAstr = sprintf('%.3g', pA); end
        e2  = a.eta2; if isempty(e2) || ~isfinite(e2), e2str = 'NaN'; else, e2str = sprintf('%.3f', e2); end
        fprintf('  one-way ANOVA: p = %s,  eta^2 = %s,  k=%d, n=%d\n', pAstr, e2str, a.k, a.n);
    end
    if isempty(k) || ~isfield(k,'p')
        fprintf('  Kruskal–Wallis: [no result]\n');
    else
        pK = k.p;  if isempty(pK) || ~isfinite(pK), pKstr = 'NaN'; else, pKstr = sprintf('%.3g', pK); end
        fprintf('  Kruskal–Wallis: p = %s\n', pKstr);
    end
end

function out = local_anova1(y, g)
    ok = isfinite(y) & ~isundefined(g);
    y  = y(ok); g = removecats(g(ok));
    if numel(categories(g)) < 2, error('Need ≥2 groups for ANOVA.'); end
    [p, tbl, stats] = anova1(y, g, 'off');
    SS_effect = tbl{2,2}; SS_error = tbl{3,2};
    eta2 = SS_effect / (SS_effect + SS_error);
    tukey = multcompare(stats, 'CType','tukey-kramer', 'Display','off');
    out = struct('p',p, 'eta2',eta2, 'tbl',{tbl}, 'tukey',tukey, ...
                 'k', numel(categories(g)), 'n', numel(y));
end

function out = local_kw(y, g)
    ok = isfinite(y) & ~isundefined(g);
    y  = y(ok); g = removecats(g(ok));
    if numel(categories(g)) < 2
        out = struct('p', NaN, 'tbl', [], 'stats', []); return
    end
    [p, tbl, stats] = kruskalwallis(y, g, 'off');
    out = struct('p',p, 'tbl',{tbl}, 'stats',stats);
end

function out = tukey_table(y, g, alpha)
    if nargin<3, alpha = 0.05; end
    ok = isfinite(y) & ~isundefined(g);
    y  = y(ok); g = removecats(g(ok));

    % --- Drop groups with n<2 (avoid singletons in post-hoc) ---
    cats = categories(g); nByCat = countcats(g);
    keepCats = cats(nByCat >= 2);
    keep = ismember(g, categorical(keepCats));
    y = y(keep); g = removecats(g(keep));
    if numel(categories(g)) < 2, out = table(); return; end

    [~, ~, stats] = anova1(y, g, 'off');
    C = multcompare(stats, 'CType','tukey-kramer', 'Alpha',alpha, 'Display','off');  % [i j lo diff hi p]
    names = local_names_from_stats(stats);
    kMeans = size(stats.means,1);
    if numel(names) ~= kMeans, names = categories(g); end
    i = C(:,1); j = C(:,2);
    keepIJ = i>=1 & i<=numel(names) & j>=1 & j<=numel(names);
    C = C(keepIJ,:); i = C(:,1); j = C(:,2);
    out = table( string(names(i)), string(names(j)), C(:,4), C(:,3), C(:,5), C(:,6), ...
        'VariableNames', {'Group1','Group2','Diff','CI_L','CI_U','pAdj'} );
    out = sortrows(out, {'pAdj','Diff'});
end

function names = local_names_from_stats(stats)
    names = strings(0,1);
    if isfield(stats,'gnames')
        if iscell(stats.gnames), names = string(stats.gnames(:));
        elseif ischar(stats.gnames), names = string(cellstr(stats.gnames));
        else, try, names = string(stats.gnames(:)); catch, names = strings(0,1); end
        end
    end
end
