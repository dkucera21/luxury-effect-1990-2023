%% === Update predictor tables with Climate/NHGIS + LME luxury-effect slopes, then Z-score ===
% Prereqs in memory: CENSUS_TABLES_rebuilt (struct), A_HYPO_MEANS, A_HYPO_MEDIANS, A_HYPO_MEANS_TEMPORAL
% Also CityListMaster (or MasterCities/CityNames as fallback).

baseDir = 'E:\LUXURY_NATL_FINAL';   % adjust if needed

%% -------------------- Settings for mixed-effects --------------------
yearsUse  = [1990 2000 2010 2020 2023];
thresh    = 50;      % PCT_OVERLAP cutoff
minTracts = 5;       % minimum valid tracts per city-year
useWinsor = false;   % optional
winzPct   = [1 99];

assert(exist('CENSUS_TABLES_rebuilt','var')==1 && isstruct(CENSUS_TABLES_rebuilt), ...
    'CENSUS_TABLES_rebuilt not found or not a struct.');

% Resolve master city list
if exist('CityListMaster','var') && ~isempty(CityListMaster)
    CityListMaster = string(CityListMaster(:));
elseif exist('MasterCities','var') && istable(MasterCities) && ismember('City', MasterCities.Properties.VariableNames)
    CityListMaster = string(MasterCities.City);
elseif exist('CityNames','var') && ~isempty(CityNames)
    CityListMaster = string(CityNames(:));
else
    error('No valid city list found (CityListMaster, MasterCities.City, or CityNames).');
end
CityListMaster = unique(strtrim(CityListMaster));

killPats = ["_SSM_", "PercMale"];   % columns containing these substrings will be dropped


%% -------------------- Ensure core A_HYPO_* are in memory --------------------
needVar = @(nm) assert(evalin('base',sprintf('exist(''%s'',''var'')==1',nm)), '%s must exist.', nm);
needVar('A_HYPO_MEANS');           A_HYPO_MEANS.City = string(A_HYPO_MEANS.City);
needVar('A_HYPO_MEDIANS');         A_HYPO_MEDIANS.City = string(A_HYPO_MEDIANS.City);
needVar('A_HYPO_MEANS_TEMPORAL');  A_HYPO_MEANS_TEMPORAL.City = string(A_HYPO_MEANS_TEMPORAL.City);
if ~ismember('Year',A_HYPO_MEANS_TEMPORAL.Properties.VariableNames)
    error('A_HYPO_MEANS_TEMPORAL must contain a Year column.');
end
A_HYPO_MEANS_TEMPORAL.Year = double(A_HYPO_MEANS_TEMPORAL.Year);

%% -------------------- Read external CSVs (Climate + NHGIS) --------------------
rd  = @(name) local_read_csv_if_exists(baseDir, name);
jC  = @(A,B) local_outerjoin_city(A,B);
jCY = @(A,B) local_outerjoin_cityyear(A,B);

C_MEANS   = rd('A_CLIMATE_HYPO_MEANS.csv');
C_MEDIANS = rd('A_CLIMATE_HYPO_MEDIANS.csv');
C_TEMP    = rd('A_CLIMATE_HYPO_MEANS_TEMPORAL.csv');
C_TR      = rd('A_CLIMATE_HYPO_TRENDS.csv');
C_TRP     = rd('A_CLIMATE_HYPO_TRENDS_pvals.csv');

N_MEANS   = rd('TABLE_NHGIS_MEANS.csv');
N_MEDIANS = rd('TABLE_NHGIS_MEDIANS.csv');
N_TEMP    = rd('TABLE_NHGIS_MEANS_TEMPORAL.csv');   % “MEANS_TEMPORAL”
N_TR      = rd('TABLE_NHGIS_MEAN_TREND.csv');
N_TRP     = rd('TABLE_NHGIS_MEAN_TREND_pvals.csv');

% Join City-level
if ~isempty(C_MEANS),   A_HYPO_MEANS   = jC(A_HYPO_MEANS,   C_MEANS);   end
if ~isempty(N_MEANS),   A_HYPO_MEANS   = jC(A_HYPO_MEANS,   N_MEANS);   end
if ~isempty(C_MEDIANS), A_HYPO_MEDIANS = jC(A_HYPO_MEDIANS, C_MEDIANS); end
if ~isempty(N_MEDIANS), A_HYPO_MEDIANS = jC(A_HYPO_MEDIANS, N_MEDIANS); end

% Join (City,Year)
if ~isempty(C_TEMP), A_HYPO_MEANS_TEMPORAL = jCY(A_HYPO_MEANS_TEMPORAL, C_TEMP); end
if ~isempty(N_TEMP), A_HYPO_MEANS_TEMPORAL = jCY(A_HYPO_MEANS_TEMPORAL, N_TEMP); end

assignin('base','A_HYPO_MEANS',A_HYPO_MEANS);
assignin('base','A_HYPO_MEDIANS',A_HYPO_MEDIANS);
assignin('base','A_HYPO_MEANS_TEMPORAL',A_HYPO_MEANS_TEMPORAL);

%% -------------------- Build/refresh trends from temporal if needed --------------------
if evalin('base','exist(''A_HYPO_TRENDS'',''var'')~=1 || exist(''A_HYPO_TRENDS_pvals'',''var'')~=1')
    [TRENDS, PVALS] = local_build_trends_from_temporal(A_HYPO_MEANS_TEMPORAL, 1990, 2023);
else
    TRENDS = evalin('base','A_HYPO_TRENDS');       TRENDS.City = string(TRENDS.City);
    PVALS  = evalin('base','A_HYPO_TRENDS_pvals'); PVALS.City  = string(PVALS.City);
end

% Merge external trend CSVs (if provided)
if ~isempty(C_TR),   TRENDS = jC(TRENDS, C_TR); end
if ~isempty(C_TRP),  PVALS  = jC(PVALS,  C_TRP); end
if ~isempty(N_TR),   TRENDS = jC(TRENDS, N_TR); end
if ~isempty(N_TRP),  PVALS  = jC(PVALS,  N_TRP); end

assignin('base','A_HYPO_TRENDS',       TRENDS);
assignin('base','A_HYPO_TRENDS_pvals', PVALS);

% ---- Drop disallowed columns from trends before prefix/append ----
A_HYPO_TRENDS = local_drop_vars_by_pattern(A_HYPO_TRENDS, killPats, "City");
assignin('base','A_HYPO_TRENDS',A_HYPO_TRENDS);


%% -------------------- Add mixed-effects luxury-effect slopes (per $10k) --------------------
% Uses tract-level data in CENSUS_TABLES_rebuilt; coverage=PCT_OVERLAP; income=RAW_INCOME_CPI.
% Adds: LME_LE_MEAN_NDVI_per10k, LME_LE_MEAN_LST_per10k (right after City).

% Ensure A_HYPO_TRENDS exists and has City
if ~(exist('A_HYPO_TRENDS','var')==1 && istable(A_HYPO_TRENDS)), A_HYPO_TRENDS = table(CityListMaster,'VariableNames',{'City'}); end
A_HYPO_TRENDS.City = string(A_HYPO_TRENDS.City);
A_HYPO_TRENDS = local_outerjoin_city(table(CityListMaster,'VariableNames',{'City'}), A_HYPO_TRENDS);

% Vars needed by the persistent function used inside fit_metric_percity
pctNames    = {'PCT_OVERLAP'};   %#ok<NASGU>  % harmonized spec
incRawCands = {'RAW_INCOME_CPI'};    %#ok<NASGU>
tblNames    = string(fieldnames(CENSUS_TABLES_rebuilt));  %#ok<NASGU>
numify      = @(x) str2double(regexprep(string(x),'[,\$%]','')); %#ok<NASGU>

% Run LME for MEAN_NDVI and MEAN_LST
pc_ndvi = fit_metric_percity("MEAN_NDVI");
pc_lst  = fit_metric_percity("MEAN_LST");

pc_ndvi.Properties.VariableNames{'slope_per10k'} = 'LME_LE_MEAN_NDVI_per10k';
pc_lst .Properties.VariableNames{'slope_per10k'} = 'LME_LE_MEAN_LST_per10k';

% Drop existing if present and join
for nm = ["LME_LE_MEAN_NDVI_per10k","LME_LE_MEAN_LST_per10k"]
    if ismember(nm, string(A_HYPO_TRENDS.Properties.VariableNames))
        A_HYPO_TRENDS(:, nm) = [];
    end
end
A_HYPO_TRENDS = local_outerjoin_city(A_HYPO_TRENDS, pc_ndvi);
A_HYPO_TRENDS = local_outerjoin_city(A_HYPO_TRENDS, pc_lst);
A_HYPO_TRENDS = movevars(A_HYPO_TRENDS, {'LME_LE_MEAN_NDVI_per10k','LME_LE_MEAN_LST_per10k'}, 'After', 'City');

assignin('base','A_HYPO_TRENDS',A_HYPO_TRENDS);
fprintf('Mixed-effects slopes added to A_HYPO_TRENDS.\n');

%% -------------------- Drop all-NaN trends cols & prefix trends, then append MEAN_* --------------------
[A_HYPO_TRENDS, killed] = local_drop_all_nan_cols(A_HYPO_TRENDS, "City");
assignin('base','A_HYPO_TRENDS',A_HYPO_TRENDS);
if exist('A_HYPO_TRENDS_pvals','var')==1 && istable(A_HYPO_TRENDS_pvals)
    A_HYPO_TRENDS_pvals(:, intersect(killed, string(A_HYPO_TRENDS_pvals.Properties.VariableNames))) = [];
    assignin('base','A_HYPO_TRENDS_pvals',A_HYPO_TRENDS_pvals);
end

% Prefix trends columns with TRENDS_ (numeric, non-constant), then append MEAN_*
[A_HYPO_TRENDS] = local_prefix_trends_and_append_means(A_HYPO_TRENDS, A_HYPO_MEANS);
assignin('base','A_HYPO_TRENDS', A_HYPO_TRENDS);

% ---- Final pass: drop any remaining matches (e.g., MEAN_* or TRENDS_* carrying patterns) ----
A_HYPO_TRENDS = local_drop_vars_by_pattern(A_HYPO_TRENDS, killPats, "City");
assignin('base','A_HYPO_TRENDS', A_HYPO_TRENDS);


%% -------------------- Build Z-scored versions with mean imputation --------------------
A_HYPO_MEANS_Z   = local_zscore_table(A_HYPO_MEANS,   keyVars="City");
A_HYPO_MEDIANS_Z = local_zscore_table(A_HYPO_MEDIANS, keyVars="City");
A_HYPO_TRENDS_Z  = local_zscore_table(A_HYPO_TRENDS,  keyVars="City");

assignin('base','A_HYPO_MEANS_Z',   A_HYPO_MEANS_Z);
assignin('base','A_HYPO_MEDIANS_Z', A_HYPO_MEDIANS_Z);
assignin('base','A_HYPO_TRENDS_Z',  A_HYPO_TRENDS_Z);

% ---- Drop disallowed columns from base city/city-year tables ----
A_HYPO_MEANS           = local_drop_vars_by_pattern(A_HYPO_MEANS,           killPats, "City");
A_HYPO_MEDIANS         = local_drop_vars_by_pattern(A_HYPO_MEDIANS,         killPats, "City");
A_HYPO_MEANS_TEMPORAL  = local_drop_vars_by_pattern(A_HYPO_MEANS_TEMPORAL,  killPats, ["City","Year"]);

assignin('base','A_HYPO_MEANS',A_HYPO_MEANS);
assignin('base','A_HYPO_MEDIANS',A_HYPO_MEDIANS);
assignin('base','A_HYPO_MEANS_TEMPORAL',A_HYPO_MEANS_TEMPORAL);


fprintf('Z-scored tables created: A_HYPO_MEANS_Z, A_HYPO_MEDIANS_Z, A_HYPO_TRENDS_Z.\n');

%% =============================== local fns ===============================
function T = local_read_csv_if_exists(baseDir, fname)
    fp = fullfile(baseDir, fname);
    if exist(fp,'file')~=2, T = []; return; end
    T = readtable(fp);
    if ismember('City', T.Properties.VariableNames), T.City = string(T.City); end
    if ismember('Year', T.Properties.VariableNames), T.Year = double(T.Year); end
end

function A = local_outerjoin_city(A,B)
    if isempty(B), return; end
    A.City = string(A.City); B.City = string(B.City);
    ov = intersect(string(A.Properties.VariableNames), string(B.Properties.VariableNames));
    ov = setdiff(ov, "City");
    if ~isempty(ov), A(:, ov) = []; end
    A = outerjoin(A, B, 'Keys','City', 'MergeKeys',true, 'Type','left');
end

function A = local_outerjoin_cityyear(A,B)
    if isempty(B), return; end
    A.City = string(A.City); B.City = string(B.City);
    A.Year = double(A.Year); B.Year = double(B.Year);
    ov = intersect(string(A.Properties.VariableNames), string(B.Properties.VariableNames));
    ov = setdiff(ov, ["City","Year"]);
    if ~isempty(ov), A(:, ov) = []; end
    A = outerjoin(A, B, 'Keys',{'City','Year'}, 'MergeKeys',true, 'Type','left');
end

function [SLOPE_TAB, PVAL_TAB] = local_build_trends_from_temporal(T, y0, y1)
    T = T(T.Year>=y0 & T.Year<=y1,:);
    T.City = string(T.City); T.Year = double(T.Year);
    vars = string(T.Properties.VariableNames);
    vars = setdiff(vars, ["City","Year"], 'stable');
    keep = false(size(vars));
    for i=1:numel(vars)
        c = T.(vars{i});
        keep(i) = (isnumeric(c)||islogical(c)||isduration(c)) && (nnz(isfinite(c))>0);
    end
    vars = vars(keep);
    cities = unique(T.City);
    nC = numel(cities);

    S = array2table(nan(nC,numel(vars)), 'VariableNames', matlab.lang.makeValidName(vars));
    P = array2table(nan(nC,numel(vars)), 'VariableNames', matlab.lang.makeValidName(vars));

    for ci=1:nC
        TT = T(T.City==cities(ci),:);
        if height(TT)<5 || numel(unique(TT.Year))<5, continue; end
        t_dn = datenum(datetime(TT.Year,1,1));
        for vi=1:numel(vars)
            y = TT.(vars{vi}); if isduration(y), y=seconds(y); end; y = double(y);
            ok = isfinite(t_dn) & isfinite(y);
            if nnz(ok)<5 || std(y(ok))==0, continue; end
            mdl = fitlm(t_dn(ok), y(ok));
            S{ci,vi} = mdl.Coefficients.Estimate(2) * 365.25; % annual slope
            P{ci,vi} = mdl.Coefficients.pValue(2);
        end
    end
    SLOPE_TAB = table(cities, 'VariableNames',{'City'}); SLOPE_TAB = [SLOPE_TAB S];
    PVAL_TAB  = table(cities, 'VariableNames',{'City'}); PVAL_TAB  = [PVAL_TAB  P];
end

function [T, killed] = local_drop_all_nan_cols(T, keyName)
    T.City = string(T.City);
    vars = string(T.Properties.VariableNames);
    cand = setdiff(vars, keyName, 'stable');
    isAllNaNCol = false(size(cand));
    for i = 1:numel(cand)
        col = T.(cand{i});
        if isnumeric(col) || islogical(col)
            isAllNaNCol(i) = all(~isfinite(col));
        elseif isduration(col)
            isAllNaNCol(i) = all(~isfinite(seconds(col)));
        elseif isdatetime(col)
            isAllNaNCol(i) = all(~isfinite(datenum(col)));
        else
            isAllNaNCol(i) = false;
        end
    end
    killed = cand(isAllNaNCol);
    if ~isempty(killed), T(:, killed) = []; end
end

function T = local_prefix_trends_and_append_means(Ttr, Tmn)
    % Prefix numeric, non-constant trend columns with TRENDS_, then append MEAN_* from A_HYPO_MEANS
    Ttr.City = string(Ttr.City);
    Tmn.City = string(Tmn.City);

    isNumericLike = @(x) isnumeric(x) || islogical(x) || isduration(x) || isdatetime(x);
    toNum = @(x) (isduration(x) .* seconds(x)) + (isdatetime(x) .* datenum(x)) + (~(isduration(x)||isdatetime(x)) .* double(x));
    isConstant = @(col) (~isNumericLike(col)) || (all(~isfinite(toNum(col)))) || ((max(toNum(col),[],'omitnan') - min(toNum(col),[],'omitnan'))==0);

    vTr = string(Ttr.Properties.VariableNames);
    candTr = setdiff(vTr, "City", 'stable');
    keepRename = false(size(candTr));
    for i = 1:numel(candTr)
        col = Ttr.(candTr{i});
        keepRename(i) = isNumericLike(col) && ~isConstant(col);
    end
    renameList = candTr(keepRename);
    if ~isempty(renameList)
        newNames = "TRENDS_" + renameList;
        Ttr = renamevars(Ttr, renameList, newNames);
    end

    vMn = string(Tmn.Properties.VariableNames);
    candMn = setdiff(vMn, "City", 'stable');
    keepMn = false(size(candMn));
    for i = 1:numel(candMn)
        col = Tmn.(candMn{i});
        keepMn(i) = isNumericLike(col) && ~isConstant(col);
    end
    useMn = candMn(keepMn);

    if isempty(useMn)
        T = Ttr; return;
    end
    Tmn_sub = Tmn(:, ["City", useMn]);
    Tmn_prefixed = Tmn_sub;
    Tmn_prefixed.Properties.VariableNames(2:end) = cellstr("MEAN_" + useMn);

    dup = intersect(string(Tmn_prefixed.Properties.VariableNames), string(Ttr.Properties.VariableNames));
    dup = setdiff(dup, "City");
    if ~isempty(dup), Ttr(:, dup) = []; end

    T = outerjoin(Ttr, Tmn_prefixed, 'Keys','City', 'MergeKeys',true, 'Type','left');
    T = movevars(T, 'City', 'Before', 1);
end

% ======================= ROBUST LME FUNCTION (unchanged logic) =======================
function perCity = fit_metric_percity(metricName)
    % Access outer-scope variables
    persistent YEARS CITYLIST TBLNAMES CENSUS pctNames incRawCands numify thresh minTracts useWinsor winzPct
    if isempty(YEARS)
        YEARS       = evalin('base','yearsUse');
        CITYLIST    = evalin('base','CityListMaster');
        TBLNAMES    = evalin('base','tblNames');
        CENSUS      = evalin('base','CENSUS_TABLES_rebuilt');
        pctNames    = evalin('base','pctNames');
        incRawCands = evalin('base','incRawCands');
        numify      = evalin('base','numify');
        thresh      = evalin('base','thresh');
        minTracts   = evalin('base','minTracts');
        useWinsor   = evalin('base','useWinsor');
        winzPct     = evalin('base','winzPct');
    end

    CityCol = {}; YearCol = []; IncCol = []; YCol = [];

    % -------- Build pooled tall table across years & cities ----------
    for ci = 1:numel(CITYLIST)
        city = CITYLIST(ci);
        for yr = YEARS
            tname = sprintf('T_%d_%s', yr, city);
            if ~ismember(tname, TBLNAMES), continue; end
            T = CENSUS.(tname); if ~istable(T) || height(T)==0, continue; end

            V = string(T.Properties.VariableNames);
            vPct = first_match_ci(V, pctNames);           if vPct == "", continue; end
            vY   = first_match_ci(V, {char(metricName)}); if vY   == "", continue; end
            vInc = first_match_ci(V, incRawCands);        if vInc == "", continue; end

            pct = numify(T.(vPct));
            mxc = max(pct(isfinite(pct)));
            if isfinite(mxc) && mxc <= 1.5, pct = 100*pct; end
            pct(pct > 100) = 100;

            x10 = numify(T.(vInc)) ./ 10000;   % per $10k
            yy  = numify(T.(vY));

            ok = isfinite(pct) & pct>=thresh & isfinite(x10) & isfinite(yy);
            if nnz(ok) < minTracts, continue; end

            x = winsor_if(x10(ok), useWinsor, winzPct);
            y = winsor_if(yy(ok),  useWinsor, winzPct);
            if std(x,'omitnan')==0, continue; end

            n = numel(x);
            CityCol(end+1:end+n,1) = repmat({char(city)}, n, 1); %#ok<AGROW>
            YearCol(end+1:end+n,1) = repmat(yr,   n, 1);        %#ok<AGROW>
            IncCol(end+1:end+n,1)  = x;                         %#ok<AGROW>
            YCol(end+1:end+n,1)    = y;                         %#ok<AGROW>
        end
    end

    % If nothing valid, return NaNs
    if isempty(YearCol)
        perCity = table(string(CITYLIST), NaN(numel(CITYLIST),1), 'VariableNames', {'City','slope_per10k'});
        return
    end

    Tall = table(categorical(CityCol), YearCol, IncCol, YCol, ...
                 'VariableNames', {'City','year','inc10k','y'});
    Tall.City = categorical(string(Tall.City), CITYLIST);
    Tall = Tall(~ismissing(Tall.City), :);

    if numel(unique(Tall.year)) < 2
        perCity = table(string(CITYLIST), NaN(numel(CITYLIST),1), 'VariableNames', {'City','slope_per10k'});
        return
    end

    Tall.yearC = Tall.year - mean(Tall.year,'omitnan');

    % -------- Fit LME (try random slope for inc10k; fallback to RI) ----------
    hasSlope = true;
    try
        mdl = fitlme(Tall, 'y ~ inc10k*yearC + (1 + inc10k | City)');
    catch
        mdl = fitlme(Tall, 'y ~ inc10k*yearC + (1 | City)');
        hasSlope = false;
    end

    % Fixed-effect income slope
    feNames = string(mdl.CoefficientNames);
    jInc = find(feNames=="inc10k" | contains(lower(feNames),"inc10k"), 1);
    b_inc = NaN; if ~isempty(jInc), b_inc = mdl.Coefficients.Estimate(jInc); end

    % Try to extract random-slope BLUPs; if absent or empty, mark hasSlope=false
    re_inc_map = containers.Map('KeyType','char','ValueType','double');
    if hasSlope
        REraw = randomEffects(mdl);
        REtab = local_to_table(REraw);
        REtab = local_harmonize_re_table(REtab);

        nameCol = string(REtab.Name);
        want = contains(lower(nameCol), "inc10k");  % tolerant across releases
        REinc = REtab(want, :);

        if isempty(REinc)
            hasSlope = false;
        else
            lev = string(REinc.Level);
            est = double(REinc.Estimate);
            for r = 1:numel(est)
                key = lower(strtrim(lev(r)));
                if ~isempty(key)
                    re_inc_map(key) = est(r);
                end
            end
        end
    end

    % -------- Compose per-city slope --------
    perCity = table(string(CITYLIST), NaN(numel(CITYLIST),1), 'VariableNames', {'City','slope_per10k'});

    if hasSlope && isfinite(b_inc)
        % FE + RE for cities that have a BLUP; else fallback OLS for that city
        for i = 1:numel(CITYLIST)
            c = string(CITYLIST(i));
            key = lower(strtrim(c));
            if isKey(re_inc_map, key)
                perCity.slope_per10k(i) = b_inc + re_inc_map(key);
            else
                perCity.slope_per10k(i) = local_city_ols(Tall, c); % fallback
            end
        end
    else
        % No random slope survived → per-city OLS fallback for everyone
        for i = 1:numel(CITYLIST)
            perCity.slope_per10k(i) = local_city_ols(Tall, string(CITYLIST(i)));
        end
    end
end

function name = first_match_ci(V, cands)
    V = string(V); VL = lower(V);
    cands = string(cands);
    for i = 1:numel(cands)
        j = find(VL == lower(cands(i)), 1, 'first');
        if ~isempty(j), name = V(j); return; end
    end
    name = "";
end

function out = winsor_if(x, doIt, lim)
    if ~doIt, out = x; return; end
    if numel(x) < 50, lim = [5 95]; end
    lo = prctile(x, lim(1)); hi = prctile(x, lim(2));
    out = min(max(x, lo), hi);
end

% ----------- ROBUST conversion helpers for randomEffects output -----------
function T = local_to_table(RE)
    % Safely coerce randomEffects output to a table (or empty table)
    if isempty(RE)
        T = table(); return;
    end
    if istable(RE)
        T = RE; return;
    end
    if isa(RE,'dataset') %#ok<DISASA>
        T = dataset2table(RE); return;
    end
    if isstruct(RE)
        try
            T = struct2table(RE); return;
        catch
            T = table(); return;
        end
    end
    % Unknown type → empty table
    T = table();
end

function Tbl = local_harmonize_re_table(T)
    % Return a table(Name, Level, Estimate) robust across MATLAB releases
    if isempty(T) || width(T)==0 || height(T)==0
        Tbl = table(string.empty(0,1), string.empty(0,1), double.empty(0,1), ...
                    'VariableNames', {'Name','Level','Estimate'});
        return
    end

    v = string(T.Properties.VariableNames);

    % --- Name column ---
    idxName = find(ismember(lower(v), lower(["Name","Effect","EffName","REName"])), 1);
    if isempty(idxName)
        % first text-like column
        isTextCol = varfun(@(x) iscellstr(x) || isstring(x) || iscategorical(x), T, 'OutputFormat','uniform');
        idxName = find(isTextCol,1);
        if isempty(idxName), idxName = 1; end
    end
    varName_Name = char(v(idxName));  % char for maximum compatibility
    Name = string(T.(varName_Name));

    % --- Level column ---
    idxLev = find(ismember(lower(v), lower(["Level","Group"])), 1);
    if isempty(idxLev)
        isTextCol = varfun(@(x) iscellstr(x) || isstring(x) || iscategorical(x), T, 'OutputFormat','uniform');
        isTextCol(idxName) = false;
        idxLev = find(isTextCol,1);
        if isempty(idxLev), idxLev = idxName; end
    end
    varName_Lev = char(v(idxLev));
    Level = string(T.(varName_Lev));

    % --- Estimate column ---
    idxEst = find(ismember(lower(v), lower(["Estimate","RE","Value","Est"])), 1);
    if isempty(idxEst)
        isNumCol = varfun(@(x) isnumeric(x), T, 'OutputFormat','uniform');
        idxEst = find(isNumCol,1);
        if isempty(idxEst)
            Estimate = double.empty(0,1);
        else
            Estimate = double(T.(char(v(idxEst))));
        end
    else
        Estimate = double(T.(char(v(idxEst))));
    end

    Tbl = table(Name, Level, Estimate);
end

function b = local_city_ols(Tall, cityName)
    m = strcmpi(string(Tall.City), string(cityName));
    if nnz(m) < 3 || std(Tall.inc10k(m),'omitnan')==0
        b = NaN; return
    end
    X = [ones(nnz(m),1), Tall.inc10k(m)];
    y = Tall.y(m);
    beta = X \ y;
    b = beta(2);
end

function T = local_drop_vars_by_pattern(T, patterns, keyVars)
    % Remove columns whose names match any of the given patterns (case-insensitive).
    % keyVars are preserved even if they match a pattern.
    if isempty(T), return; end
    v = string(T.Properties.VariableNames);
    keyVars = string(keyVars(:));
    work = setdiff(v, keyVars, 'stable');
    kill = false(size(work));
    for i = 1:numel(work)
        namei = work(i);
        for p = 1:numel(patterns)
            if ~isempty(regexpi(namei, patterns(p), 'once'))
                kill(i) = true; break;
            end
        end
    end
    T(:, work(kill)) = [];
end


function TZ = local_zscore_table(T, varargin)
    % Mean-impute missing values, then z-score each numeric/logical/duration/datetime column.
    % Skip: key columns (keyVars), categorical, constant columns.
    p = inputParser; addParameter(p,'keyVars',"City"); parse(p,varargin{:});
    keyVars = string(p.Results.keyVars);
    TZ = T;  % copy
    vars = string(T.Properties.VariableNames);
    work = setdiff(vars, keyVars, 'stable');

    isNumericLike = @(x) isnumeric(x) || islogical(x) || isduration(x) || isdatetime(x);

    for i = 1:numel(work)
        vn = work(i);
        col = T.(vn);

        % skip non-numeric-like and categoricals
        if iscategorical(col) || ~isNumericLike(col), continue; end

        % coerce to numeric vector
        if isduration(col), x = seconds(col);
        elseif isdatetime(col), x = datenum(col);
        else, x = double(col);
        end

        % constant? skip
        xFinite = x(isfinite(x));
        if isempty(xFinite) || (max(xFinite) - min(xFinite) == 0)
            continue
        end

        % mean impute
        mu = mean(xFinite,'omitnan');
        x(~isfinite(x)) = mu;

        % standardize
        sd = std(x, 'omitnan');
        if sd>0, xz = (x - mu) ./ sd; else, xz = zeros(size(x)); end

        % write back
        TZ.(vn) = xz;
    end
end
