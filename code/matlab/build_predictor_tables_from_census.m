%% build_predictor_tables_from_census.m
% Builds predictor tables from CENSUS_TABLES_new (in base workspace)
% Filters cities STRICTLY: a city is included only if *every* city-year
% present for that city has >= minTracts tracts with coverage >= threshold.
% Also excludes specified non-CONUS cities.

%% Pull inputs from base workspace
try
    CENSUS_TABLES_new = evalin('base','CENSUS_TABLES_new');
catch
    error('CENSUS_TABLES_new not found in base workspace. Load it first, then run this script.');
end

% If these predictor tables already exist in base, pull them; otherwise create placeholders
A_HYPO_MEANS            = fetchOrEmpty('A_HYPO_MEANS');
A_HYPO_MEDIANS          = fetchOrEmpty('A_HYPO_MEDIANS');
A_HYPO_MEANS_TEMPORAL   = fetchOrEmpty('A_HYPO_MEANS_TEMPORAL');
A_HYPO_MEDIANS_TEMPORAL = fetchOrEmpty('A_HYPO_MEDIANS_TEMPORAL');

%% ---- configuration ----
threshold   = 50;    % % overlap cutoff (0–100) per tract
minTracts   = 5;     % minimum valid tracts per city-year (by coverage only)
pctNames    = {'PCT_OVERLAP','PercentAreaNDVI','PCT_COVER','PCT_COV'}; % coverage candidates
excludeList = string({'Anchorage','Fairbanks','HiloHI','Honolulu','Ponce','SanJuan'});

commonVars = { ...
  'POP_DENSITY','HOUSING_DENSITY','PERCENT_UNDER_5', ...
  'PERCENT_OCCUPIED_HOUSING','PERCENT_OWNED_HOUSING', ...
  'PERCENT_HISPANIC','PERCENT_WHITE','PERCENT_BLACK', ...
  'PERCENT_NATIVE','PERCENT_ASIAN','PERCENT_BACHELOR', ...
  'PERCENT_GRADUATE','PERCENT_HS','incZ'};

meanOnlyVars   = {'MEAN_NDVI','MEAN_LST'};
medianOnlyVars = {'MEDIAN_NDVI','MEDIAN_LST'};

varsForMeansTemporal   = [commonVars, meanOnlyVars];
varsForMediansTemporal = [commonVars, medianOnlyVars];

%% ---- pass 0: determine strict inclusion per city (every year must pass) ----
tblNames = fieldnames(CENSUS_TABLES_new);
passRows = []; % struct array: City, Year, Pass (logical)

for i = 1:numel(tblNames)
    tname = tblNames{i};
    tok = regexp(tname, '^T_(\d{4})n?_(.+)$', 'tokens', 'once');
    if isempty(tok), continue; end
    yr   = str2double(tok{1});
    city = string(tok{2});

    T = CENSUS_TABLES_new.(tname);
    if ~istable(T) || height(T)==0, continue; end

    pctColName = findFirstVar(T.Properties.VariableNames, pctNames);
    if isempty(pctColName)
        pass = false; % no coverage column → treat as failing this year
    else
        pct = getNumCol(T, pctColName);
        pass = nnz(pct >= threshold & isfinite(pct)) >= minTracts;
    end
    passRows = [passRows; struct('City',city,'Year',yr,'Pass',logical(pass))]; %#ok<AGROW>
end

if isempty(passRows)
    error('Could not parse any city-year tables from CENSUS_TABLES_new (expected names like T_YYYY_CITY).');
end

PassTab = struct2table(passRows);
% Strict rule: include city only if ALL its years present pass
G   = findgroups(PassTab.City);
ok  = splitapply(@all, PassTab.Pass, G);
CU  = splitapply(@(c) c(1), PassTab.City, G);
includeCities = CU(ok);

% Remove excluded cities explicitly
includeCities = setdiff(includeCities, excludeList);

% If nothing passes, fail fast
if isempty(includeCities)
    error('No cities meet the strict coverage rule across all city-years (after exclusions).');
end

%% ---- pass 1: per-city-year summaries (we still compute, then filter to includeCities) ----
rowsMean   = []; % struct array for temporal means
rowsMedian = []; % struct array for temporal medians

for i = 1:numel(tblNames)
    tname = tblNames{i};
    tok = regexp(tname, '^T_(\d{4})n?_(.+)$', 'tokens', 'once');
    if isempty(tok), continue; end
    yr   = str2double(tok{1});
    city = string(tok{2});

    T = CENSUS_TABLES_new.(tname);
    if ~istable(T) || height(T)==0, continue; end

    pctColName = findFirstVar(T.Properties.VariableNames, pctNames);
    if isempty(pctColName), continue; end
    pct = getNumCol(T, pctColName);

    % --- MEAN set (per-city-year) ---
    recM = struct('City', city, 'Year', yr);
    for v = 1:numel(varsForMeansTemporal)
        vn  = varsForMeansTemporal{v};
        col = getNumCol(T, vn);
        if isempty(col)
            recM.(vn) = NaN; continue;
        end
        valid = pct >= threshold & isfinite(col);
        recM.(vn) = iif(nnz(valid) >= minTracts, mean(col(valid), 'omitnan'), NaN);
    end
    rowsMean = [rowsMean; recM]; %#ok<AGROW>

    % --- MEDIAN set (per-city-year) ---
    recMd = struct('City', city, 'Year', yr);
    for v = 1:numel(varsForMediansTemporal)
        vn  = varsForMediansTemporal{v};
        col = getNumCol(T, vn);
        if isempty(col)
            recMd.(vn) = NaN; continue;
        end
        valid = pct >= threshold & isfinite(col);
        recMd.(vn) = iif(nnz(valid) >= minTracts, median(col(valid), 'omitnan'), NaN);
    end
    rowsMedian = [rowsMedian; recMd]; %#ok<AGROW>
end

% Convert to tables and filter to included cities
meanTemporal   = toTemporalTable(rowsMean, includeCities);
medianTemporal = toTemporalTable(rowsMedian, includeCities);

%% ---- pass 2: merge into *_TEMPORAL (left-join on City,Year) ----
A_HYPO_MEANS_TEMPORAL = ensureTable(A_HYPO_MEANS_TEMPORAL);
A_HYPO_MEDIANS_TEMPORAL = ensureTable(A_HYPO_MEDIANS_TEMPORAL);

A_HYPO_MEANS_TEMPORAL.City = string(A_HYPO_MEANS_TEMPORAL.City);
A_HYPO_MEANS_TEMPORAL.Year = double(A_HYPO_MEANS_TEMPORAL.Year);
A_HYPO_MEANS_TEMPORAL = dropIfExists(A_HYPO_MEANS_TEMPORAL, varsForMeansTemporal);
A_HYPO_MEANS_TEMPORAL = outerjoin(A_HYPO_MEANS_TEMPORAL, meanTemporal, ...
    'Keys', {'City','Year'}, 'MergeKeys', true, 'Type', 'left');
A_HYPO_MEANS_TEMPORAL = A_HYPO_MEANS_TEMPORAL(ismember(A_HYPO_MEANS_TEMPORAL.City, includeCities), :);

A_HYPO_MEDIANS_TEMPORAL.City = string(A_HYPO_MEDIANS_TEMPORAL.City);
A_HYPO_MEDIANS_TEMPORAL.Year = double(A_HYPO_MEDIANS_TEMPORAL.Year);
A_HYPO_MEDIANS_TEMPORAL = dropIfExists(A_HYPO_MEDIANS_TEMPORAL, varsForMediansTemporal);
A_HYPO_MEDIANS_TEMPORAL = outerjoin(A_HYPO_MEDIANS_TEMPORAL, medianTemporal, ...
    'Keys', {'City','Year'}, 'MergeKeys', true, 'Type', 'left');
A_HYPO_MEDIANS_TEMPORAL = A_HYPO_MEDIANS_TEMPORAL(ismember(A_HYPO_MEDIANS_TEMPORAL.City, includeCities), :);

%% ---- pass 3: collapse across years and merge into non-temporal ----
aggMean   = collapseAcrossYears(meanTemporal,   varsForMeansTemporal,   @mean);
aggMedian = collapseAcrossYears(medianTemporal, varsForMediansTemporal, @median);

A_HYPO_MEANS   = filterAndMergeNonTemporal(A_HYPO_MEANS,   aggMean,   includeCities, varsForMeansTemporal);
A_HYPO_MEDIANS = filterAndMergeNonTemporal(A_HYPO_MEDIANS, aggMedian, includeCities, varsForMediansTemporal);

%% ---- summary ----
fprintf('\nSTRICT filter: city kept only if EVERY city-year has coverage ≥ %g%% AND ≥ %d tracts.\n', threshold, minTracts);
fprintf('Excluded (hard list): %s\n', strjoin(excludeList, ', '));
fprintf('Included cities: %d\n', numel(includeCities));
fprintf('  • A_HYPO_MEANS_TEMPORAL: %d rows, +%d cols\n', ...
    height(A_HYPO_MEANS_TEMPORAL), numel(varsForMeansTemporal));
fprintf('  • A_HYPO_MEDIANS_TEMPORAL: %d rows, +%d cols\n', ...
    height(A_HYPO_MEDIANS_TEMPORAL), numel(varsForMediansTemporal));
fprintf('  • A_HYPO_MEANS: %d rows, +%d cols (city-avg across years)\n', ...
    height(A_HYPO_MEANS), numel(varsForMeansTemporal));
fprintf('  • A_HYPO_MEDIANS: %d rows, +%d cols (city-median across years)\n', ...
    height(A_HYPO_MEDIANS), numel(varsForMediansTemporal));

%% ---- push results back to base ----
assignin('base','A_HYPO_MEANS',A_HYPO_MEANS);
assignin('base','A_HYPO_MEDIANS',A_HYPO_MEDIANS);
assignin('base','A_HYPO_MEANS_TEMPORAL',A_HYPO_MEANS_TEMPORAL);
assignin('base','A_HYPO_MEDIANS_TEMPORAL',A_HYPO_MEDIANS_TEMPORAL);
assignin('base','CityNames',includeCities);

clearvars -except CENSUS_TABLES_new A_HYPO_MEANS A_HYPO_MEANS_TEMPORAL A_HYPO_MEDIANS A_HYPO_MEDIANS_TEMPORAL

%% ====== local helpers ======
function T = fetchOrEmpty(varName)
    try
        T = evalin('base', varName);
        if ~istable(T), T = table(string.empty(0,1), zeros(0,1), 'VariableNames', {'City','Year'}); end
    catch
        T = table(string.empty(0,1), zeros(0,1), 'VariableNames', {'City','Year'});
    end
end

function name = findFirstVar(vnames, candidates)
    name = '';
    for k = 1:numel(candidates)
        j = find(strcmpi(vnames, candidates{k}), 1, 'first');
        if ~isempty(j), name = vnames{j}; return; end
    end
end

function col = getNumCol(T, name)
    vnames = T.Properties.VariableNames;
    j = find(strcmpi(vnames, name), 1, 'first');
    if isempty(j)
        col = [];
        return;
    end
    raw = T.(vnames{j});
    if isnumeric(raw)
        col = double(raw);
    else
        col = str2double(regexprep(string(raw), '[,\$%]', ''));
    end
    % Treat sentinel -222222222 as missing
    col(col == -222222222) = NaN;
end

function out = iif(cond, a, b), if cond, out=a; else, out=b; end, end

function TT = toTemporalTable(rows, includeCities)
    if isempty(rows)
        TT = table(string.empty(0,1), zeros(0,1), 'VariableNames', {'City','Year'});
    else
        TT = struct2table(rows);
        TT.City = string(TT.City);
        TT.Year = double(TT.Year);
        if ~isempty(includeCities)
            TT = TT(ismember(TT.City, includeCities), :);
        end
    end
end

function TT = ensureTable(TT)
    if isempty(TT)
        TT = table(string.empty(0,1), zeros(0,1), 'VariableNames', {'City','Year'});
    end
end

function agg = collapseAcrossYears(TT, varList, reducerFcn)
    if isempty(TT) || height(TT)==0
        agg = table(string.empty(0,1), 'VariableNames', {'City'});
        return
    end
    G = findgroups(TT.City);
    CU = splitapply(@(c) c(1), TT.City, G);
    agg = table(CU, 'VariableNames', {'City'});
    for v = 1:numel(varList)
        vn = varList{v};
        if ismember(vn, TT.Properties.VariableNames)
            if isequal(func2str(reducerFcn),'median')
                agg.(vn) = splitapply(@(x) median(x,'omitnan'), TT.(vn), G);
            else
                agg.(vn) = splitapply(@(x) mean(x,'omitnan'), TT.(vn), G);
            end
        else
            agg.(vn) = NaN(height(agg),1);
        end
    end
end

function A = filterAndMergeNonTemporal(A, agg, includeCities, varList)
    A = ensureCityOnly(A);
    if ~isempty(includeCities)
        A = A(ismember(string(A.City), includeCities), :);
        agg = agg(ismember(string(agg.City), includeCities), :);
    end
    A = dropIfExists(A, varList);
    A = outerjoin(A, agg, 'Keys', 'City', 'MergeKeys', true, 'Type', 'left');
end

function A = ensureCityOnly(A)
    if isempty(A)
        A = table(string.empty(0,1), 'VariableNames', {'City'});
        return
    end
    if ~ismember('City', A.Properties.VariableNames)
        error('Non-temporal table is missing "City" column.');
    end
    A.City = string(A.City);
end

function A = dropIfExists(A, names)
    if isempty(A) || isempty(names), return; end
    vA = A.Properties.VariableNames;
    kill = false(1, numel(vA));
    for i = 1:numel(vA)
        if any(strcmpi(vA{i}, names))
            kill(i) = true;
        end
    end
    A(:, kill) = [];
end
