%% build_predictor_tables_from_census_PLUS_LUX.m
% Uses ONLY columns that exist in the harmonized tables:
% PCT_OVERLAP, CENSUS_TRACT_AREA,
% HOUSING_DENSITY, POP_DENSITY, PERCENT_OCCUPIED_HOUSING, PERCENT_OWNED_HOUSING,
% PERCENT_HISPANIC, PERCENT_WHITE, PERCENT_BLACK, PERCENT_ASIAN,
% PERCENT_BACHELOR, PERCENT_GRADUATE, PERCENT_HS,
% RAW_INCOME_CPI, MEAN_NDVI, MEAN_LST, MEDIAN_NDVI, MEDIAN_LST.

%% -------------------------- Requirements --------------------------
assert(exist('CENSUS_TABLES_rebuilt','var')==1 && isstruct(CENSUS_TABLES_rebuilt), ...
    'CENSUS_TABLES_rebuilt (struct) must be in the workspace.');
assert(exist('CityListMaster','var')==1, ...
    'CityListMaster must be in the workspace.');

CityListMaster = unique(string(strtrim(CityListMaster(:))));
nCities = numel(CityListMaster);
fprintf('\n[START] Building predictor tables for %d cities ...\n', nCities);

%% -------------------------- Config --------------------------
CoverageThreshold = 50;  % PCT_OVERLAP cutoff (0–100)
MinTractsPerCY    = 5;
YearPattern       = '^T_(\d{4})n?_(.+)$';

% Hard-locked variable list (present in harmonized) — NO Shape_* fields
ALLVARS = ["PCT_OVERLAP","CENSUS_TRACT_AREA", ...
           "HOUSING_DENSITY","POP_DENSITY","PERCENT_OCCUPIED_HOUSING","PERCENT_OWNED_HOUSING", ...
           "PERCENT_HISPANIC","PERCENT_WHITE","PERCENT_BLACK","PERCENT_ASIAN", ...
           "PERCENT_BACHELOR","PERCENT_GRADUATE","PERCENT_HS","RAW_INCOME",...
           "RAW_INCOME_CPI","MEAN_NDVI","MEAN_LST","MEDIAN_NDVI","MEDIAN_LST"];

NUMVARS_TO_AGG = setdiff(ALLVARS, "PCT_OVERLAP");   % aggregate everything except coverage

%% -------------------------- Build (City,Year) shell --------------------------
fns = string(fieldnames(CENSUS_TABLES_rebuilt));
tok = regexp(fns, YearPattern, 'tokens', 'once');
has = ~cellfun(@isempty, tok);
YearsAll  = cellfun(@(t) str2double(t{1}), tok(has));
CitiesAll = string(cellfun(@(t) t{2}, tok(has), 'UniformOutput', false));
keep = ismember(CitiesAll, CityListMaster);
CY   = table(CitiesAll(keep), YearsAll(keep), 'VariableNames', {'City','Year'});
CY   = unique(sortrows(CY, {'City','Year'}), 'rows');
CY.City = string(CY.City); CY.Year = double(CY.Year);

%% -------------------------- Ensure shells exist --------------------------
A_HYPO_MEANS            = pf_ensure_city_shell('A_HYPO_MEANS',           CityListMaster);
A_HYPO_MEDIANS          = pf_ensure_city_shell('A_HYPO_MEDIANS',         CityListMaster);
A_HYPO_MEANS_TEMPORAL   = pf_ensure_temporal_shell('A_HYPO_MEANS_TEMPORAL',   CY);
A_HYPO_MEDIANS_TEMPORAL = pf_ensure_temporal_shell('A_HYPO_MEDIANS_TEMPORAL', CY);

%% -------------------------- Aggregate per city-year --------------------------
rowsMEAN = []; rowsMEDN = [];

for i = 1:numel(fns)
    tname = fns(i);
    tk = regexp(tname, YearPattern, 'tokens', 'once');
    if isempty(tk), continue; end
    yr   = str2double(tk{1});
    city = string(tk{2});
    if ~ismember(city, CityListMaster), continue; end

    T = CENSUS_TABLES_rebuilt.(tname);
    if ~istable(T) || height(T)==0, continue; end

    % Coverage (only column that exists)
    if ~ismember("PCT_OVERLAP", string(T.Properties.VariableNames)), continue; end
    cov = double(T.PCT_OVERLAP);

    recMean = struct('City',city,'Year',yr);
    recMed  = struct('City',city,'Year',yr);
    for vn = NUMVARS_TO_AGG(:).'
        vchar = char(vn);
        if ismember(vn, string(T.Properties.VariableNames))
            val = double(T.(vchar));
            val(val == -222222222) = NaN;
            ok = isfinite(cov) & cov >= CoverageThreshold & isfinite(val);
            if nnz(ok) >= MinTractsPerCY
                recMean.(vchar) = mean(val(ok), 'omitnan');
                recMed.(vchar)  = median(val(ok), 'omitnan');
            else
                recMean.(vchar) = NaN; recMed.(vchar) = NaN;
            end
        else
            recMean.(vchar) = NaN; recMed.(vchar) = NaN;
        end
    end

    rowsMEAN = [rowsMEAN; recMean]; %#ok<AGROW>
    rowsMEDN = [rowsMEDN; recMed];  %#ok<AGROW>
end

MEAN_T = pf_struct2temporal(rowsMEAN);
MEDN_T = pf_struct2temporal(rowsMEDN);

% Keep exact city set
if ~isempty(MEAN_T), MEAN_T = MEAN_T(ismember(MEAN_T.City, CityListMaster), :); end
if ~isempty(MEDN_T), MEDN_T = MEDN_T(ismember(MEDN_T.City, CityListMaster), :); end

%% -------- incZ_among_cities (cross-city Z of RAW_INCOME_CPI per year) --------
MEAN_T.incZ_among_cities = NaN(height(MEAN_T),1);
MEDN_T.incZ_among_cities = NaN(height(MEDN_T),1);

if ismember('RAW_INCOME_CPI', string(MEAN_T.Properties.VariableNames))
    for yy = unique(MEAN_T.Year).'
        m = MEAN_T.Year==yy & isfinite(MEAN_T.RAW_INCOME_CPI);
        mu = mean(MEAN_T.RAW_INCOME_CPI(m),'omitnan');
        sd = std( MEAN_T.RAW_INCOME_CPI(m),'omitnan');
        if isfinite(sd) && sd>0
            MEAN_T.incZ_among_cities(m) = (MEAN_T.RAW_INCOME_CPI(m) - mu) ./ sd;
        end
    end
end
if ismember('RAW_INCOME_CPI', string(MEDN_T.Properties.VariableNames))
    for yy = unique(MEDN_T.Year).'
        m = MEDN_T.Year==yy & isfinite(MEDN_T.RAW_INCOME_CPI);
        mu = mean(MEDN_T.RAW_INCOME_CPI(m),'omitnan');
        sd = std( MEDN_T.RAW_INCOME_CPI(m),'omitnan');
        if isfinite(sd) && sd>0
            MEDN_T.incZ_among_cities(m) = (MEDN_T.RAW_INCOME_CPI(m) - mu) ./ sd;
        end
    end
end

%% -------------------------- Luxury-effect slopes --------------------------
% slope of outcome ~ (RAW_INCOME_CPI / 10,000) within city-year
rowsLuxMean = []; rowsLuxMed  = [];
for i = 1:numel(fns)
    tname = fns(i);
    tk = regexp(tname, YearPattern, 'tokens', 'once');
    if isempty(tk), continue; end
    yr   = str2double(tk{1});
    city = string(tk{2});
    if ~ismember(city, CityListMaster), continue; end

    T = CENSUS_TABLES_rebuilt.(tname);
    if ~istable(T) || height(T)==0, continue; end
    V = string(T.Properties.VariableNames);

    if ~all(ismember(["PCT_OVERLAP","RAW_INCOME_CPI"], V)), continue; end
    cov = double(T.PCT_OVERLAP);
    inc = double(T.RAW_INCOME_CPI);
    inc(~isfinite(inc)) = NaN;
    x10k = inc ./ 10000;

    s_mean_ndvi_10k = pf_slope_if_ok(x10k, pf_numcol(T,V,'MEAN_NDVI'),   cov, CoverageThreshold, MinTractsPerCY);
    s_mean_lst_10k  = pf_slope_if_ok(x10k, pf_numcol(T,V,'MEAN_LST'),    cov, CoverageThreshold, MinTractsPerCY);
    s_med_ndvi_10k  = pf_slope_if_ok(x10k, pf_numcol(T,V,'MEDIAN_NDVI'), cov, CoverageThreshold, MinTractsPerCY);
    s_med_lst_10k   = pf_slope_if_ok(x10k, pf_numcol(T,V,'MEDIAN_LST'),  cov, CoverageThreshold, MinTractsPerCY);

    rowsLuxMean = [rowsLuxMean; struct('City',city,'Year',yr, ...
        'LUX_NDVI_per10k', s_mean_ndvi_10k, 'LUX_LST_per10k',  s_mean_lst_10k)]; %#ok<AGROW>
    rowsLuxMed  = [rowsLuxMed;  struct('City',city,'Year',yr, ...
        'LUX_NDVI_per10k', s_med_ndvi_10k,  'LUX_LST_per10k',  s_med_lst_10k)];  %#ok<AGROW>
end
LUX_MEAN_T = pf_struct2temporal(rowsLuxMean);
LUX_MEDN_T = pf_struct2temporal(rowsLuxMed);

%% -------------------------- Vegetative cooling slopes ----------------------
% slope of LST on NDVI
rowsVCMean = []; rowsVCMed  = [];
for i = 1:numel(fns)
    tname = fns(i);
    tk = regexp(tname, YearPattern, 'tokens', 'once');
    if isempty(tk), continue; end
    yr   = str2double(tk{1});
    city = string(tk{2});
    if ~ismember(city, CityListMaster), continue; end

    T = CENSUS_TABLES_rebuilt.(tname);
    if ~istable(T) || height(T)==0, continue; end
    V = string(T.Properties.VariableNames);
    if ~ismember("PCT_OVERLAP", V), continue; end

    cov    = double(T.PCT_OVERLAP);
    ndmean = pf_numcol(T,V,'MEAN_NDVI');
    lsmean = pf_numcol(T,V,'MEAN_LST');
    ndmed  = pf_numcol(T,V,'MEDIAN_NDVI');
    lsmed  = pf_numcol(T,V,'MEDIAN_LST');

    b_mean = pf_pair_slope_if_ok(ndmean, lsmean, cov, CoverageThreshold, MinTractsPerCY);
    b_med  = pf_pair_slope_if_ok(ndmed,  lsmed,  cov, CoverageThreshold, MinTractsPerCY);

    rowsVCMean = [rowsVCMean; struct('City',city,'Year',yr,'VEGCOOL_MEAN_LST_perNDVI',   b_mean)]; %#ok<AGROW>
    rowsVCMed  = [rowsVCMed;  struct('City',city,'Year',yr,'VEGCOOL_MEDIAN_LST_perNDVI', b_med)];  %#ok<AGROW>
end
VC_MEAN_T = pf_struct2temporal(rowsVCMean);
VC_MEDN_T = pf_struct2temporal(rowsVCMed);

%% -------------------------- Merge temporal tables --------------------------
drop_MEAN_names = setdiff(string(MEAN_T.Properties.VariableNames),   ["City","Year"]);
drop_MEDN_names = setdiff(string(MEDN_T.Properties.VariableNames),   ["City","Year"]);

A_HYPO_MEANS_TEMPORAL   = pf_dropIfExists(A_HYPO_MEANS_TEMPORAL,   drop_MEAN_names);
A_HYPO_MEDIANS_TEMPORAL = pf_dropIfExists(A_HYPO_MEDIANS_TEMPORAL, drop_MEDN_names);

A_HYPO_MEANS_TEMPORAL   = outerjoin(A_HYPO_MEANS_TEMPORAL,   MEAN_T,     'Keys',{'City','Year'}, 'MergeKeys',true, 'Type','left');
A_HYPO_MEDIANS_TEMPORAL = outerjoin(A_HYPO_MEDIANS_TEMPORAL, MEDN_T,     'Keys',{'City','Year'}, 'MergeKeys',true, 'Type','left');

% add LUX
A_HYPO_MEANS_TEMPORAL   = pf_dropIfExists(A_HYPO_MEANS_TEMPORAL,   setdiff(string(LUX_MEAN_T.Properties.VariableNames), ["City","Year"]));
A_HYPO_MEDIANS_TEMPORAL = pf_dropIfExists(A_HYPO_MEDIANS_TEMPORAL, setdiff(string(LUX_MEDN_T.Properties.VariableNames), ["City","Year"]));
A_HYPO_MEANS_TEMPORAL   = outerjoin(A_HYPO_MEANS_TEMPORAL,   LUX_MEAN_T, 'Keys',{'City','Year'}, 'MergeKeys',true, 'Type','left');
A_HYPO_MEDIANS_TEMPORAL = outerjoin(A_HYPO_MEDIANS_TEMPORAL, LUX_MEDN_T, 'Keys',{'City','Year'}, 'MergeKeys',true, 'Type','left');

% add VEGCOOL
A_HYPO_MEANS_TEMPORAL   = pf_dropIfExists(A_HYPO_MEANS_TEMPORAL,   setdiff(string(VC_MEAN_T.Properties.VariableNames), ["City","Year"]));
A_HYPO_MEDIANS_TEMPORAL = pf_dropIfExists(A_HYPO_MEDIANS_TEMPORAL, setdiff(string(VC_MEDN_T.Properties.VariableNames), ["City","Year"]));
A_HYPO_MEANS_TEMPORAL   = outerjoin(A_HYPO_MEANS_TEMPORAL,   VC_MEAN_T, 'Keys',{'City','Year'}, 'MergeKeys',true, 'Type','left');
A_HYPO_MEDIANS_TEMPORAL = outerjoin(A_HYPO_MEDIANS_TEMPORAL, VC_MEDN_T, 'Keys',{'City','Year'}, 'MergeKeys',true, 'Type','left');

% Keep exact city set
A_HYPO_MEANS_TEMPORAL   = A_HYPO_MEANS_TEMPORAL(  ismember(string(A_HYPO_MEANS_TEMPORAL.City),   CityListMaster), :);
A_HYPO_MEDIANS_TEMPORAL = A_HYPO_MEDIANS_TEMPORAL(ismember(string(A_HYPO_MEDIANS_TEMPORAL.City), CityListMaster), :);

% RAW_INCOME_CPI to end if present
A_HYPO_MEANS_TEMPORAL   = pf_move_to_end(A_HYPO_MEANS_TEMPORAL,   "RAW_INCOME_CPI");
A_HYPO_MEDIANS_TEMPORAL = pf_move_to_end(A_HYPO_MEDIANS_TEMPORAL, "RAW_INCOME_CPI");

% Also move RAW_INCOME_CPI to end in MEAN_T/MEDN_T shells
MEAN_T = pf_move_to_end(MEAN_T, "RAW_INCOME_CPI");
MEDN_T = pf_move_to_end(MEDN_T, "RAW_INCOME_CPI");

%% -------------------------- Collapse to non-temporal -----------------------
vars_MEAN_all = setdiff(string(A_HYPO_MEANS_TEMPORAL.Properties.VariableNames),   ["City","Year"]);
vars_MEDN_all = setdiff(string(A_HYPO_MEDIANS_TEMPORAL.Properties.VariableNames), ["City","Year"]);

aggMEAN   = pf_collapseAcrossYears(A_HYPO_MEANS_TEMPORAL,   vars_MEAN_all, @mean);
aggMEDIAN = pf_collapseAcrossYears(A_HYPO_MEDIANS_TEMPORAL, vars_MEDN_all, @median);

% 2023 snapshot of incZ_among_cities (if present)
if any(A_HYPO_MEANS_TEMPORAL.Year==2023) && ismember('incZ_among_cities', string(A_HYPO_MEANS_TEMPORAL.Properties.VariableNames))
    Z23m = unique(A_HYPO_MEANS_TEMPORAL(A_HYPO_MEANS_TEMPORAL.Year==2023, {'City','incZ_among_cities'}), 'rows');
    Z23m.Properties.VariableNames{'incZ_among_cities'} = 'incZ_among_cities_2023';
    aggMEAN = outerjoin(aggMEAN, Z23m, 'Keys','City', 'MergeKeys',true, 'Type','left');
end
if any(A_HYPO_MEDIANS_TEMPORAL.Year==2023) && ismember('incZ_among_cities', string(A_HYPO_MEDIANS_TEMPORAL.Properties.VariableNames))
    Z23d = unique(A_HYPO_MEDIANS_TEMPORAL(A_HYPO_MEDIANS_TEMPORAL.Year==2023, {'City','incZ_among_cities'}), 'rows');
    Z23d.Properties.VariableNames{'incZ_among_cities'} = 'incZ_among_cities_2023';
    aggMEDIAN = outerjoin(aggMEDIAN, Z23d, 'Keys','City', 'MergeKeys',true, 'Type','left');
end

% Collapse LUX and VEGCOOL by city
luxVars = ["LUX_NDVI_per10k","LUX_LST_per10k"];
luxMEAN_byCity = pf_collapseAcrossYears(A_HYPO_MEANS_TEMPORAL,   luxVars, @mean);
luxMEDN_byCity = pf_collapseAcrossYears(A_HYPO_MEDIANS_TEMPORAL, luxVars, @median);

vcMEAN_byCity  = pf_collapseAcrossYears(VC_MEAN_T, "VEGCOOL_MEAN_LST_perNDVI",   @mean);
vcMEDN_byCity  = pf_collapseAcrossYears(VC_MEDN_T, "VEGCOOL_MEDIAN_LST_perNDVI", @median);

% Merge into collapsed tables
aggMEAN   = pf_dropIfExists(aggMEAN,   [luxVars, "VEGCOOL_MEAN_LST_perNDVI"]);
aggMEDIAN = pf_dropIfExists(aggMEDIAN, [luxVars, "VEGCOOL_MEDIAN_LST_perNDVI"]);

aggMEAN   = outerjoin(aggMEAN,   luxMEAN_byCity, 'Keys','City', 'MergeKeys',true, 'Type','left');
aggMEAN   = outerjoin(aggMEAN,   vcMEAN_byCity,  'Keys','City', 'MergeKeys',true, 'Type','left');

aggMEDIAN = outerjoin(aggMEDIAN, luxMEDN_byCity, 'Keys','City', 'MergeKeys',true, 'Type','left');
aggMEDIAN = outerjoin(aggMEDIAN, vcMEDN_byCity,  'Keys','City', 'MergeKeys',true, 'Type','left');

% Merge onto non-temporal shells
A_HYPO_MEANS   = pf_dropIfExists(A_HYPO_MEANS,   string(aggMEAN.Properties.VariableNames));
A_HYPO_MEDIANS = pf_dropIfExists(A_HYPO_MEDIANS, string(aggMEDIAN.Properties.VariableNames));
A_HYPO_MEANS   = outerjoin(A_HYPO_MEANS,   aggMEAN,   'Keys','City', 'MergeKeys',true, 'Type','left');
A_HYPO_MEDIANS = outerjoin(A_HYPO_MEDIANS, aggMEDIAN, 'Keys','City', 'MergeKeys',true, 'Type','left');

% RAW_INCOME_CPI to end
A_HYPO_MEANS   = pf_move_to_end(A_HYPO_MEANS,   "RAW_INCOME_CPI");
A_HYPO_MEDIANS = pf_move_to_end(A_HYPO_MEDIANS, "RAW_INCOME_CPI");

%% ---------- Attach City metadata (if available) ----------
metaPath = fullfile(pwd, 'City_Preliminary_Vars.csv');
if exist(metaPath,'file') == 2
    meta = readtable(metaPath);

    % Expect at least: City, BIOMES_cat, KOPPEN_cat, Latitude, Longitude
    v = string(meta.Properties.VariableNames);
    if ~ismember("City", v), error('City_Preliminary_Vars.csv must contain a City column.'); end
    meta.City = string(meta.City);

    % Keep only known columns (ignore extras)
    keep = intersect(["City","BIOMES_cat","KOPPEN_cat","Latitude","Longitude"], v, 'stable');
    meta = meta(:, keep);

    % Ensure types
    if ismember("BIOMES_cat", keep) && ~iscategorical(meta.BIOMES_cat), meta.BIOMES_cat = categorical(meta.BIOMES_cat); end
    if ismember("KOPPEN_cat", keep) && ~iscategorical(meta.KOPPEN_cat), meta.KOPPEN_cat = categorical(meta.KOPPEN_cat); end
    if ismember("Latitude",  keep),  meta.Latitude  = double(meta.Latitude);  end
    if ismember("Longitude", keep),  meta.Longitude = double(meta.Longitude); end

    addMeta = @(A) local_leftjoin_city(A, meta);

    if exist('A_HYPO_MEANS','var')==1,             A_HYPO_MEANS             = addMeta(A_HYPO_MEANS);             end
    if exist('A_HYPO_MEDIANS','var')==1,           A_HYPO_MEDIANS           = addMeta(A_HYPO_MEDIANS);           end
    if exist('A_HYPO_MEANS_TEMPORAL','var')==1,    A_HYPO_MEANS_TEMPORAL    = addMeta(A_HYPO_MEANS_TEMPORAL);    end
    if exist('A_HYPO_MEDIANS_TEMPORAL','var')==1,  A_HYPO_MEDIANS_TEMPORAL  = addMeta(A_HYPO_MEDIANS_TEMPORAL);  end

    fprintf('Merged BIOMES_cat, KOPPEN_cat, Latitude, Longitude from %s\n', metaPath);
else
    warning('City_Preliminary_Vars.csv not found in %s — skipping metadata merge.', pwd);
end

%% -------------------------- Save to base & report --------------------------
assignin('base','A_HYPO_MEANS',            A_HYPO_MEANS);
assignin('base','A_HYPO_MEDIANS',          A_HYPO_MEDIANS);
assignin('base','A_HYPO_MEANS_TEMPORAL',   A_HYPO_MEANS_TEMPORAL);
assignin('base','A_HYPO_MEDIANS_TEMPORAL', A_HYPO_MEDIANS_TEMPORAL);

fprintf('\n--- Predictor tables built (%d cities) ---\n', numel(CityListMaster));
fprintf('A_HYPO_MEANS            : %d rows, %d vars\n', height(A_HYPO_MEANS),            width(A_HYPO_MEANS));
fprintf('A_HYPO_MEDIANS          : %d rows, %d vars\n', height(A_HYPO_MEDIANS),          width(A_HYPO_MEDIANS));
fprintf('A_HYPO_MEANS_TEMPORAL   : %d rows, %d vars\n', height(A_HYPO_MEANS_TEMPORAL),   width(A_HYPO_MEANS_TEMPORAL));
fprintf('A_HYPO_MEDIANS_TEMPORAL : %d rows, %d vars\n', height(A_HYPO_MEDIANS_TEMPORAL), width(A_HYPO_MEDIANS_TEMPORAL));

%% ====================== Local helpers (pf_*) ======================
function T = pf_ensure_city_shell(varName, CityList)
    Seed = table(CityList, 'VariableNames', {'City'});
    Seed.City = string(Seed.City);
    if evalin('base', sprintf('exist(''%s'',''var'')==1', varName))
        T0 = evalin('base', varName);
        if istable(T0) && ismember('City', T0.Properties.VariableNames)
            T0.City = string(T0.City);
            T = outerjoin(Seed, T0, 'Keys','City', 'MergeKeys',true, 'Type','left');
        else
            T = Seed;
        end
    else
        T = Seed;
    end
    T = movevars(T, 'City', 'Before', 1);
end

function T = pf_ensure_temporal_shell(varName, CY)
    Seed = CY;
    Seed.City = string(Seed.City);
    Seed.Year = double(Seed.Year);
    if evalin('base', sprintf('exist(''%s'',''var'')==1', varName))
        T0 = evalin('base', varName);
        if istable(T0) && all(ismember({'City','Year'}, T0.Properties.VariableNames))
            T0.City = string(T0.City);
            T0.Year = double(T0.Year);
            T = outerjoin(Seed, T0, 'Keys', {'City','Year'}, 'MergeKeys',true, 'Type','left');
        else
            T = Seed;
        end
    else
        T = Seed;
    end
    T = movevars(T, {'City','Year'}, 'Before', 1);
end

function agg = pf_collapseAcrossYears(TT, varList, reducerFcn)
% Collapse a City–Year table to one row per City using mean/median, etc.
% - TT must contain columns: City, Year, and vars in varList (if present).
% - reducerFcn: e.g., @mean or @median (applied with 'omitnan').

    if isempty(TT) || ~istable(TT)
        agg = table(string.empty(0,1), 'VariableNames', {'City'});
        return
    end

    % Ensure types
    if ~isstring(TT.City), TT.City = string(TT.City); end
    if ~isnumeric(TT.Year), TT.Year = double(TT.Year); end

    % One row per City
    G   = findgroups(TT.City);
    CU  = splitapply(@(c) c(1), TT.City, G);
    agg = table(CU, 'VariableNames', {'City'});

    % Normalize var list
    varList = string(varList(:));

    % Reducer that omits NaNs
    fnum = @(x) reducerFcn(x, 'omitnan');

    for k = 1:numel(varList)
        vn = varList(k);
        if ~ismember(vn, string(TT.Properties.VariableNames))
            agg.(vn) = NaN(height(agg),1);
            continue
        end

        col = TT.(vn);
        % Accept numeric/logical only (your harmonized tables are numeric)
        if isnumeric(col) || islogical(col)
            agg.(vn) = splitapply(fnum, double(col), G);
        else
            % Non-numeric → fill with NaN
            agg.(vn) = NaN(height(agg),1);
        end
    end
end

function A = pf_move_to_end(A, varName)
    if isempty(A) || ~istable(A), return; end
    v = string(A.Properties.VariableNames);
    hit = find(strcmpi(v, string(varName)), 1, 'first');
    if ~isempty(hit) && hit < numel(v)
        A = movevars(A, v(hit), 'After', v(end));
    end
end

function y = pf_numcol(T, V, name)
    if ismember(name, V)
        raw = T.(char(name));
        y = double(raw);
        y(y == -222222222) = NaN;
    else
        y = NaN(height(T),1);
    end
end

function b1 = pf_slope_if_ok(x10k, y, cov, covThr, nMin)
    ok = isfinite(x10k) & isfinite(y) & isfinite(cov) & cov >= covThr;
    if nnz(ok) < nMin || var(x10k(ok),'omitnan')==0, b1 = NaN; return; end
    x = x10k(ok); yy = y(ok);
    xb = mean(x); yb = mean(yy);
    Sxx = sum((x - xb).^2);
    if Sxx == 0, b1 = NaN; else, b1 = sum((x - xb).*(yy - yb)) / Sxx; end
end

function b1 = pf_pair_slope_if_ok(x_ndvi, y_lst, cov, covThr, nMin)
    ok = isfinite(x_ndvi) & isfinite(y_lst) & isfinite(cov) & cov >= covThr;
    if nnz(ok) < nMin || var(x_ndvi(ok),'omitnan')==0, b1 = NaN; return; end
    x = x_ndvi(ok); yy = y_lst(ok);
    xb = mean(x); yb = mean(yy);
    Sxx = sum((x - xb).^2);
    if Sxx == 0, b1 = NaN; else, b1 = sum((x - xb).*(yy - yb)) / Sxx; end
end

function TT = pf_struct2temporal(rows)
    if isempty(rows)
        TT = table(string.empty(0,1), zeros(0,1), 'VariableNames', {'City','Year'});
    else
        TT = struct2table(rows);
        if ~ismember('City', TT.Properties.VariableNames), TT.City = string(TT.City); end
        TT.City = string(TT.City);
        if ~ismember('Year', TT.Properties.VariableNames), TT.Year = double(TT.Year); end
        TT.Year = double(TT.Year);
    end
end

function A = pf_dropIfExists(A, names)
    if isempty(A) || isempty(names), return; end
    vA = string(A.Properties.VariableNames);
    names = string(names);
    names = setdiff(names, ["City","Year"]);
    kill = ismember(vA, names);
    A(:, kill) = [];
end

function A = local_leftjoin_city(A, meta)
    if isempty(A) || ~istable(A), return; end
    if ~ismember('City', A.Properties.VariableNames), return; end
    A.City = string(A.City);
    meta.City = string(meta.City);
    % Drop pre-existing metadata columns to avoid suffixes
    drop = intersect({'BIOMES_cat','KOPPEN_cat','Latitude','Longitude'}, string(A.Properties.VariableNames));
    if ~isempty(drop), A(:, drop) = []; end
    A = outerjoin(A, meta, 'Keys','City', 'MergeKeys',true, 'Type','left');
end
