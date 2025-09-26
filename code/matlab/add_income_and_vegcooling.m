%% add_income_and_vegcooling.m
% Adds standardized income (incZ) per city-year table and computes vegetative
% cooling (VEG_COOLING) slopes from tract NDVI/LST, producing:
%   - A_HYPO_MEDIANS_TEMPORAL (per city-year)
%   - A_HYPO_MEANS_TEMPORAL   (per city-year)
%   - A_HYPO_MEDIANS          (collapsed city-level, median of annual slopes)
%   - A_HYPO_MEANS            (collapsed city-level, mean   of annual slopes)
%
% Requirements:
%   - Variable `CENSUS_TABLES_new` in workspace (struct of tables T_YYYY_CITY)
%   - Statistics and Machine Learning Toolbox (for fitlm; falls back to OLS)
%
% Outputs:
%   - Updates variables in workspace (no assignin/clear)
%   - Optionally saves MAT (toggle saveOut at bottom)

if ~exist('CENSUS_TABLES_new','var') || ~isstruct(CENSUS_TABLES_new)
    error('CENSUS_TABLES_new not found. Run your build script first.');
end


%% =========================
%  1) Vegetative cooling slopes from tract NDVI/LST
% =========================
fprintf('\n=== Computing vegetative cooling slopes (VEG_COOLING) ===\n');

% Settings (adjust to taste)
threshold  = 50;                       % % coverage cutoff (0–100)
minTracts  = 5;                        % minimum valid tracts per city-year
pctNames   = {'PCT_OVERLAP','PercentAreaNDVI'};  % first matching column is used
useRobust  = false;                    % set true to request robust regression

% Helper to clean numeric
cleanNum = @(x) str2double(regexprep(string(x),'[,\$%]',''));

rowsMed  = [];  % struct array for MEDIAN_* path
rowsMean = [];  % struct array for MEAN_*   path

tblNames = fieldnames(CENSUS_TABLES_new);
for i = 1:numel(tblNames)
    tname = tblNames{i};

    % Expect T_YEAR_CITY or T_YEARn_CITY
    tok = regexp(tname, '^T_(\d{4})n?_(.+)$', 'tokens', 'once');
    if isempty(tok), continue; end
    yr   = str2double(tok{1});
    city = string(tok{2});

    T = CENSUS_TABLES_new.(tname);
    V = T.Properties.VariableNames;

    % Pick coverage column
    vPct = '';
    for p = 1:numel(pctNames)
        j = find(strcmpi(V, pctNames{p}), 1, 'first');
        if ~isempty(j), vPct = V{j}; break; end
    end
    if isempty(vPct), continue; end
    pct = cleanNum(T.(vPct));

    % -------- MEDIANS path --------
    jX = find(strcmpi(V,'MEDIAN_NDVI'), 1, 'first');
    jY = find(strcmpi(V,'MEDIAN_LST'),  1, 'first');
    if ~isempty(jX) && ~isempty(jY)
        x = cleanNum(T.(V{jX}));
        y = cleanNum(T.(V{jY}));
        valid = pct >= threshold & isfinite(x) & isfinite(y);
        if nnz(valid) >= minTracts && var(x(valid),'omitnan') > 0
            slope = slope_y_on_x(x(valid), y(valid), useRobust);
        else
            slope = NaN;
        end
        rowsMed = [rowsMed; struct('City',city,'Year',yr,'VEG_COOLING',slope)]; %#ok<AGROW>
    end

    % -------- MEANS path --------
    jX = find(strcmpi(V,'MEAN_NDVI'), 1, 'first');
    jY = find(strcmpi(V,'MEAN_LST'),  1, 'first');
    if ~isempty(jX) && ~isempty(jY)
        x = cleanNum(T.(V{jX}));
        y = cleanNum(T.(V{jY}));
        valid = pct >= threshold & isfinite(x) & isfinite(y);
        if nnz(valid) >= minTracts && var(x(valid),'omitnan') > 0
            slope = slope_y_on_x(x(valid), y(valid), useRobust);
        else
            slope = NaN;
        end
        rowsMean = [rowsMean; struct('City',city,'Year',yr,'VEG_COOLING',slope)]; %#ok<AGROW>
    end
end

% Convert to tables
MedTemporal  = struct2table_safe(rowsMed);
MeanTemporal = struct2table_safe(rowsMean);

% -------- Merge with *TEMPORAL predictor tables --------
% MEDIANS
if exist('A_HYPO_MEDIANS_TEMPORAL','var') && istable(A_HYPO_MEDIANS_TEMPORAL)
    A_HYPO_MEDIANS_TEMPORAL.City = string(A_HYPO_MEDIANS_TEMPORAL.City);
    A_HYPO_MEDIANS_TEMPORAL.Year = double(A_HYPO_MEDIANS_TEMPORAL.Year);
    if ismember('VEG_COOLING', A_HYPO_MEDIANS_TEMPORAL.Properties.VariableNames)
        A_HYPO_MEDIANS_TEMPORAL.VEG_COOLING = [];
    end
    A_HYPO_MEDIANS_TEMPORAL = outerjoin(A_HYPO_MEDIANS_TEMPORAL, MedTemporal, ...
        'Keys',{'City','Year'}, 'MergeKeys',true, 'Type','left');
else
    warning('A_HYPO_MEDIANS_TEMPORAL not found; creating from census slopes (MEDIANS).');
    A_HYPO_MEDIANS_TEMPORAL = MedTemporal;
end

% MEANS
if exist('A_HYPO_MEANS_TEMPORAL','var') && istable(A_HYPO_MEANS_TEMPORAL)
    A_HYPO_MEANS_TEMPORAL.City = string(A_HYPO_MEANS_TEMPORAL.City);
    A_HYPO_MEANS_TEMPORAL.Year = double(A_HYPO_MEANS_TEMPORAL.Year);
    if ismember('VEG_COOLING', A_HYPO_MEANS_TEMPORAL.Properties.VariableNames)
        A_HYPO_MEANS_TEMPORAL.VEG_COOLING = [];
    end
    A_HYPO_MEANS_TEMPORAL = outerjoin(A_HYPO_MEANS_TEMPORAL, MeanTemporal, ...
        'Keys',{'City','Year'}, 'MergeKeys',true, 'Type','left');
else
    warning('A_HYPO_MEANS_TEMPORAL not found; creating from census slopes (MEANS).');
    A_HYPO_MEANS_TEMPORAL = MeanTemporal;
end

% -------- Collapse per city and merge with non-temporal predictor tables --------
% MEDIANS → median across years
aggMed = collapse_city(MedTemporal, @median);
if exist('A_HYPO_MEDIANS','var') && istable(A_HYPO_MEDIANS)
    A_HYPO_MEDIANS.City = string(A_HYPO_MEDIANS.City);
    if ismember('VEG_COOLING', A_HYPO_MEDIANS.Properties.VariableNames)
        A_HYPO_MEDIANS.VEG_COOLING = [];
    end
    A_HYPO_MEDIANS = outerjoin(A_HYPO_MEDIANS, aggMed, 'Keys','City', 'MergeKeys',true, 'Type','left');
else
    warning('A_HYPO_MEDIANS not found; creating from MEDIANS temporal slopes.');
    A_HYPO_MEDIANS = aggMed;
end

% MEANS → mean across years
aggMean = collapse_city(MeanTemporal, @mean);
if exist('A_HYPO_MEANS','var') && istable(A_HYPO_MEANS)
    A_HYPO_MEANS.City = string(A_HYPO_MEANS.City);
    if ismember('VEG_COOLING', A_HYPO_MEANS.Properties.VariableNames)
        A_HYPO_MEANS.VEG_COOLING = [];
    end
    A_HYPO_MEANS = outerjoin(A_HYPO_MEANS, aggMean, 'Keys','City', 'MergeKeys',true, 'Type','left');
else
    warning('A_HYPO_MEANS not found; creating from MEANS temporal slopes.');
    A_HYPO_MEANS = aggMean;
end

fprintf('\nVEG_COOLING added to:\n');
fprintf('  • A_HYPO_MEDIANS_TEMPORAL (per city-year, MEDIAN_* path)\n');
fprintf('  • A_HYPO_MEANS_TEMPORAL   (per city-year, MEAN_* path)\n');
fprintf('  • A_HYPO_MEDIANS          (city-level median of annual slopes)\n');
fprintf('  • A_HYPO_MEANS            (city-level mean of annual slopes)\n');

%% =========================
%  Local helper functions
% =========================
function T = struct2table_safe(rows)
    if isempty(rows)
        T = table(string.empty(0,1), double.empty(0,1), double.empty(0,1), ...
            'VariableNames', {'City','Year','VEG_COOLING'});
    else
        T = struct2table(rows);
        if ~ismember('City', T.Properties.VariableNames)
            T.City = string(T.City);
        else
            T.City = string(T.City);
        end
        T.Year = double(T.Year);
    end
end

function agg = collapse_city(TemporalTbl, aggfun)
    if isempty(TemporalTbl)
        agg = table(string.empty(0,1),'VariableNames',{'City'});
        agg.VEG_COOLING = [];
        return
    end
    G  = findgroups(string(TemporalTbl.City));
    CU = splitapply(@(c) c(1), string(TemporalTbl.City), G);
    agg = table(CU, 'VariableNames', {'City'});
    agg.VEG_COOLING = splitapply(@(z) aggfun(z,'omitnan'), TemporalTbl.VEG_COOLING, G);
end

function b1 = slope_y_on_x(x, y, useRob)
% Returns slope of y ~ x. Uses fitlm if available; otherwise OLS.
    x = x(:); y = y(:);
    good = isfinite(x) & isfinite(y);
    x = x(good); y = y(good);
    if numel(x) < 2 || var(x)==0
        b1 = NaN; return
    end
    try
        if license('test','Statistics_Toolbox') && exist('fitlm','file') == 2
            if useRob
                mdl = fitlm(table(x,y),'y ~ x','RobustOpts','on');
            else
                mdl = fitlm(table(x,y),'y ~ x');
            end
            b1 = mdl.Coefficients.Estimate(2);
            return
        end
    catch
        % fall through to OLS
    end
    % Closed-form OLS fallback
    xbar = mean(x); ybar = mean(y);
    Sxx = sum((x - xbar).^2);
    if Sxx == 0, b1 = NaN; else, b1 = sum((x - xbar).*(y - ybar)) / Sxx; end
end