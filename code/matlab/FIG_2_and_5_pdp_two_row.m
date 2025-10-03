%% pdp_two_row.m
% Repository: luxury-effect-analysis
% Purpose: Two-row PDPs with reproducible CV + stability selection + elastic net
%          Supports either MEANS (A_HYPO_MEANS_Z) or TRENDS (A_HYPO_TRENDS_Z)
%          Uses NON-Z tables for dependent variables to preserve original units.

%% =========================
%% USER INPUT (edit these)
%% =========================
% Data source: 'MEANS' (A_HYPO_MEANS_Z / DVs from A_HYPO_MEANS)
%              'TRENDS' (A_HYPO_TRENDS_Z / DVs from CityTrends_* tables; DV=Slope)
whichTable = 'MEANS';   % 'MEANS' or 'TRENDS'

% Dependent variables (labels only; for TRENDS, values come from CityTrends_* tables)
% - For MEANS (non-z DV table A_HYPO_MEANS):   LUX_NDVI_per10k, LUX_LST_per10k
% - For TRENDS (CityTrends_* tables):          defaults to MEDIAN_* slopes (fallbacks built in)
respNDVI = '';  % leave '' to use default for chosen whichTable
respLST  = '';  % leave '' to use default for chosen whichTable

% Master city set (must exist in workspace)
% CityListMaster = ... (must be present)

% Predictor-domain rule:
% - Drop predictors from DV's own domain (drop *NDVI* when DV is NDVI; drop *LST* when DV is LST)
% - Exception: allow NDVI predictors for LST DV (classic VC path)
allowNDVIonLST = true;

%% Pathway-specific predictor excludes
% Names: exact matches.  Patterns: regular expressions (case-insensitive).
% Defaults reproduce your earlier logic (esp. NDVI anti-leak).
ndviExcludeNames    = [ "TRENDS_LME_LE_MEAN_LST_per10k","MEAN_NHGIS_PercentBLACKNHP","MEAN_NHGIS_PercentBLACK"];  % reduce redundant selected predictors
ndviExcludePatterns = [ "^LUX_NDVI($|_)", "^dLE_.*NDVI.*" ];  % anti-leak for NDVI DV

lstExcludeNames     = [ ];
lstExcludePatterns  = [ ];

%% Reproducibility & model knobs
rngSeed          = 42;           % global RNG seed
alphaEN          = 0.8;          % elastic-net mixing (0=ridge,1=lasso)
lambdaGrid       = logspace(-4,0,100);  % fixed lambda path
KfoldTarget      = 5;            % target folds (auto 2..5 with min n rules)
R_stab_LST       = 300;          % stability selection repeats for LST DV
R_stab_NDVI      = 300;          % stability selection repeats for NDVI DV
subFrac_LST      = 0.80;         % subsample frac per stability repeat (LST)
subFrac_NDVI     = 0.80;         % subsample frac per stability repeat (NDVI)
tauStable_NDVI   = 0.60;         % stability frequency threshold for NDVI DV
topCap_stepwise  = 6;            % max number of stepwise additions per DV

useParallelFits  = true;         % use parpool for CV fits when available
plotWidthPx      = 1400;         % figure width
plotHeightPx     = 800;          % figure height

%% Universal predictor exclusions (NEVER allowed as predictors)
universalExclude = [ ...
    "TRENDS_RAW_INCOME", ...
    "LUX_LST_per10k", ...
    "MEAN_LUX_LST_per10k", ...
    "RAW_INCOME", ...
    "mean_SSM_2023_pm_median", ...
    "mean_SSM_2023_am_median", ...
    "Longitude" , ...
    "era5_soil_temperature_level_2", ...
    "era5_soil_temperature_level_1", ...
    "NHGIS_PercMale", ...
	"era5_skin_temperature",... %correlate with landsat LST
    "TRENDS_Longitude", ...
    "MEAN_Longitude",...
    "TRENDS_CENSUS_TRACT_AREA", ...
	"TRENDS_LUX_LST_per10k",...
	"MEAN_CENSUS_TRACT_AREA",...
	'MEAN_NHGIS_PercOWNOCCmortgage',... % Not a reasonable predictor
	"TRENDS_LUX_LST_per10k",...
	"TRENDS_LME_TREND_dLE_MEAN_LST_per10k_perYR",...
	"TRENDS_LME_TREND_dLE_MEAN_NDVI_per10k_perYR",...
];

%% Data-source-specific soft exclusions
keepSpei = ["spei_SPEI_12_month","spei_SPEI_06_month","spei_SPEI_48_month" ];   % only allow these SPEI
dropIfContains = ["gldas"];        % drop any var containing 'gldas' (case-insensitive)

%% =========================
%% Setup & table selection
%% =========================
assert(exist('CityListMaster','var')==1 && ~isempty(CityListMaster), ...
    'CityListMaster not found or empty.');
citySet = unique(strtrim(string(CityListMaster(:))));

switch upper(strtrim(whichTable))
    case 'MEANS'
        assert(exist('A_HYPO_MEANS_Z','var')==1 && istable(A_HYPO_MEANS_Z), ...
            'A_HYPO_MEANS_Z not found or not a table.');
        assert(exist('A_HYPO_MEANS','var')==1 && istable(A_HYPO_MEANS), ...
            'A_HYPO_MEANS (non-z) not found or not a table.');
        tblZ   = A_HYPO_MEANS_Z;
        tblRAW = A_HYPO_MEANS;
        if isempty(respNDVI), respNDVI = 'LUX_NDVI_per10k'; end
        if isempty(respLST),  respLST  = 'LUX_LST_per10k';  end
        tableTag = 'A_HYPO_MEANS';

    case 'TRENDS'
        assert(exist('A_HYPO_TRENDS_Z','var')==1 && istable(A_HYPO_TRENDS_Z), ...
            'A_HYPO_TRENDS_Z not found or not a table.');
        tblZ   = A_HYPO_TRENDS_Z;

        % ---- Build tblRAW from CityTrends_* tables (use Slope) ----
        % Preference: MEDIAN_*; fallback to MEAN_*; fallback to shims CityTrends_NDVI/LST
        haveMed = exist('CityTrends_MEDIAN_NDVI','var')==1 && istable(CityTrends_MEDIAN_NDVI) && ...
                  exist('CityTrends_MEDIAN_LST','var')==1  && istable(CityTrends_MEDIAN_LST);
        haveMean= exist('CityTrends_MEAN_NDVI','var')==1   && istable(CityTrends_MEAN_NDVI)   && ...
                  exist('CityTrends_MEAN_LST','var')==1    && istable(CityTrends_MEAN_LST);
        haveShim= exist('CityTrends_NDVI','var')==1        && istable(CityTrends_NDVI)        && ...
                  exist('CityTrends_LST','var')==1         && istable(CityTrends_LST);

        if haveMed
            Tn = CityTrends_MEDIAN_NDVI; Tl = CityTrends_MEDIAN_LST;
            if isempty(respNDVI), respNDVI = 'LE_trend_NDVI'; end
            if isempty(respLST),  respLST  = 'LE_trend_LST';  end
        elseif haveMean
            Tn = CityTrends_MEAN_NDVI;   Tl = CityTrends_MEAN_LST;
            if isempty(respNDVI), respNDVI = 'LE_trend_NDVI'; end
            if isempty(respLST),  respLST  = 'LE_trend_LST';  end
        elseif haveShim
            Tn = CityTrends_NDVI;        Tl = CityTrends_LST;
            if isempty(respNDVI), respNDVI = 'LE_trend_NDVI'; end
            if isempty(respLST),  respLST  = 'LE_trend_LST';  end
        else
            error(['CityTrends_* tables not found. Expected MEDIAN or MEAN versions, ' ...
                   'or CityTrends_NDVI/LST shims.']);
        end

        % Harmonize and keep only City, Slope
        Tn = Tn(:, intersect(["City","Slope"], string(Tn.Properties.VariableNames)));
        Tl = Tl(:, intersect(["City","Slope"], string(Tl.Properties.VariableNames)));
        Tn.City = string(Tn.City); Tl.City = string(Tl.City);

        % Join to make a tiny raw DV table with named columns for NDVI and LST slopes
        J = outerjoin(Tn, Tl, 'Keys','City', 'MergeKeys', true, 'Type','left', ...
                      'LeftVariables', {'City','Slope'}, ...
                      'RightVariables', {'Slope'});
        % Rename Slope columns to user labels
        J.Properties.VariableNames{'Slope_Tn'} = respNDVI; % older MATLAB names auto like Slope_Tn/Slope_Tl
        if ~ismember(respNDVI, string(J.Properties.VariableNames))
            % Handle cases where name became Slope_left
            vn = string(J.Properties.VariableNames);
            jN = find(contains(vn,'Slope') & contains(vn, {'Tn','left'}, 'IgnoreCase',true), 1);
            if ~isempty(jN), J.Properties.VariableNames{jN} = respNDVI; end
        end
        J.Properties.VariableNames{'Slope_Tl'} = respLST;
        if ~ismember(respLST, string(J.Properties.VariableNames))
            vn = string(J.Properties.VariableNames);
            jL = find(contains(vn,'Slope') & contains(vn, {'Tl','right'}, 'IgnoreCase',true), 1);
            if ~isempty(jL), J.Properties.VariableNames{jL} = respLST; end
        end

        % Final DV table aligned later in build_xy
        tblRAW = table(J.City, J.(respNDVI), J.(respLST), 'VariableNames', {'City', respNDVI, respLST});

        tableTag = 'A_HYPO_TRENDS';

    otherwise
        error('whichTable must be ''MEANS'' or ''TRENDS''.');
end

% Keep valid-city rows in Z (predictors)
assert(ismember('City', tblZ.Properties.VariableNames), 'Predictor table needs a City column.');
tblZ = tblZ(ismember(strtrim(string(tblZ.City)), citySet), :);

% Harmonize NON-Z DV table City
assert(ismember('City', tblRAW.Properties.VariableNames), 'NON-Z table needs a City column.');
tblRAW.City = string(tblRAW.City);

% RNG + parallel
rng(rngSeed, 'twister');
if useParallelFits && isempty(gcp('nocreate')), try parpool; catch, end, end
optsPar = statset('UseParallel', ~isempty(gcp('nocreate')));
optsSS  = statset('UseParallel', false);   % keep stability loop deterministic
getKfold = @(m) max(2, min(5, floor(m/2))); %#ok<NASGU>

%% =========================
%% Build exclusion masks (top-level)
%% =========================
allVars0 = string(tblZ.Properties.VariableNames);

% SPEI rule: keep only keepSpei, drop other spei*
isSpei   = contains(lower(allVars0), 'spei');
speiAll  = allVars0(isSpei);
speiKeep = speiAll(ismember(lower(speiAll), lower(keepSpei)));
speiDrop = setdiff(speiAll, speiKeep, 'stable');

% GLDAS drop
gldasDrop = allVars0( contains(lower(allVars0), lower(dropIfContains)) );

% Universal list (intersect with present variables)
univDrop = intersect(universalExclude, allVars0, 'stable');

baseExclusions = unique([ "City","BIOMES_cat","KOPPEN_cat", speiDrop, gldasDrop, univDrop ], 'stable');

%% EXTRA LEAKAGE GUARD FOR TRENDS: drop any luxury-effect–derived predictors
if strcmpi(whichTable,'TRENDS')
    allVars0 = string(tblZ.Properties.VariableNames);

    % Known exact names you already flagged + a few common variants
    leLeakNames = [
        "TRENDS_LME_TREND_dLE_MEAN_LST_per10k_perYR"
        "TRENDS_LME_TREND_dLE_MEAN_NDVI_per10k_perYR"
        "TRENDS_LME_LE_MEAN_LST_per10k"
        "TRENDS_LME_LE_MEAN_NDVI_per10k"
        "LUX_NDVI_per10k"
        "LUX_LST_per10k"
        "TRENDS_LUX_NDVI_per10k"
        "TRENDS_LUX_LST_per10k"
    ];

    % Regex patterns to catch *any* luxury-effect derived feature names
    % (covers dLE_, LE_, LME_LE_, LUX_* and variants for NDVI/LST/MEAN/MEDIAN, trend, perYR, etc.)
    leLeakRegex = [
        "(^|_)dLE_"
        "(^|_)LE_"
        "(^|_)LME_LE_"
        "(^|_)LUX_(NDVI|LST)"
    ];

    leakMask = ismember(allVars0, leLeakNames);
    for p = leLeakRegex(:)'
        leakMask = leakMask | ~cellfun('isempty', regexpi(allVars0, char(p), 'once'));
    end

    leakDrop = allVars0(leakMask);
    if ~isempty(leakDrop)
        baseExclusions = unique([baseExclusions(:); leakDrop(:)], 'stable');
        fprintf('[Guard] Dropped %d LE-derived predictors in TRENDS mode:\n  %s\n', ...
            numel(leakDrop), strjoin(cellstr(leakDrop), ', '));
    end
end


%% =========================
%% Helper: Build predictors & y for a given response
%% =========================
function [tblPred, yRaw, predsAll] = build_xy( ...
        tblZloc, tblRAWloc, resp, baseEx, allowNDVIonLST, roleStr, ...
        ndviExcludeNames, ndviExcludePatterns, lstExcludeNames, lstExcludePatterns)

    % y from NON-Z table, aligned by City to tblZ rows
    if ~ismember(resp, string(tblRAWloc.Properties.VariableNames))
        error('Dependent variable %s not found in NON-Z table.', resp);
    end
    Z   = tblZloc(:, 'City');
    RAW = tblRAWloc(:, {'City', char(resp)});
    J   = outerjoin(Z, RAW, 'Keys','City', 'MergeKeys', true, 'Type','left');
    yRaw = double(J.(char(resp)));

    % ---- Assemble exclusion sets ----
    allVars = string(tblZloc.Properties.VariableNames);

    % Domain-based (“drop DV’s own domain”) — keep NDVI for LST if allowed
    respLower = lower(resp);
    dropDomain = strings(0,1);
    if contains(respLower, 'lst')
        dropDomain = allVars(contains(lower(allVars),'lst'));  % drop LST-domain preds
        if allowNDVIonLST
            dropDomain = setdiff(dropDomain, allVars(contains(lower(allVars),'ndvi')), 'stable');
        end
    elseif contains(respLower, 'ndvi')
        dropDomain = allVars(contains(lower(allVars),'ndvi')); % drop NDVI-domain preds
    end

    % Pathway-specific excludes (names + regex patterns)
    role = upper(string(roleStr));   % 'NDVI' or 'LST'
    dropNames = strings(0,1);
    dropRegex = false(size(allVars));

    switch role
        case "NDVI"
            dropNames = intersect(string(ndviExcludeNames), allVars, 'stable');
            for p = ndviExcludePatterns(:)'
                dropRegex = dropRegex | ~cellfun('isempty', regexpi(allVars, char(p), 'once'));
            end
        case "LST"
            dropNames = intersect(string(lstExcludeNames), allVars, 'stable');
            for p = lstExcludePatterns(:)'
                dropRegex = dropRegex | ~cellfun('isempty', regexpi(allVars, char(p), 'once'));
            end
        otherwise
            error('roleStr must be ''NDVI'' or ''LST''.');
    end
    dropByPattern = allVars(dropRegex);

    % ---- Final candidate set ----
    preds0 = setdiff(allVars, unique([ ...
                    baseEx(:); ...
                    dropDomain(:); ...
                    dropNames(:); ...
                    dropByPattern(:); ...
                    string(resp) ]), 'stable');

    % Keep numeric, non-degenerate
    keep = true(size(preds0));
    for i = 1:numel(preds0)
        v = tblZloc.(preds0{i});
        if ~isnumeric(v), keep(i)=false; continue; end
        s = std(v,'omitnan'); if ~isfinite(s) || s==0, keep(i)=false; end
    end
    predsAll = preds0(keep);

    % Keep rows with finite y
    rowKeep = isfinite(yRaw);
    tblPred = tblZloc(rowKeep, :);
    yRaw    = yRaw(rowKeep);

    % Drop any predictors that became all-NaN after row filter
    if ~isempty(predsAll)
        X = tblPred{:, predsAll};
        colFinite = false(size(predsAll));
        for j = 1:numel(predsAll), colFinite(j) = any(isfinite(X(:,j))); end
        predsAll = predsAll(colFinite);
    end
end

%% =========================
%% LST DV (bottom row) — NDVI predictors allowed
%% =========================
[tbl_LSTsrc, y_LST, preds_LST] = build_xy( ...
    tblZ, tblRAW, respLST, baseExclusions, true, 'LST', ...
    ndviExcludeNames, ndviExcludePatterns, lstExcludeNames, lstExcludePatterns);

if isempty(preds_LST), error('No usable predictors for %s after exclusions.', respLST); end
X1 = tbl_LSTsrc{:, preds_LST}; n1 = height(tbl_LSTsrc);

% Stability selection for LST DV
countsSel1 = zeros(1, numel(preds_LST)); runsEff1 = 0;
for r = 1:R_stab_LST
    rng(1000+r, 'twister');
    nSub = max(5, round(subFrac_LST*n1));
    if nSub >= n1
        idxSub = (1:n1)';                  % no subsample if tiny n
    else
        idxSub = randsample(n1, nSub, false);
    end

    if numel(idxSub) < 3, continue; end
    Kfold_eff = max(2, min(5, floor(numel(idxSub)/2)));
    try
        cvSS = cvpartition(numel(idxSub), 'KFold', Kfold_eff);
        [B,F] = lasso(X1(idxSub,:), y_LST(idxSub), ...
                      'CV', cvSS, 'Lambda', lambdaGrid, ...
                      'Alpha', alphaEN, 'Standardize', false, ...
                      'Options', optsSS);
        idx = F.Index1SE;
        if ~isempty(idx)
            countsSel1 = countsSel1 + (B(:,idx) ~= 0)'; runsEff1 = runsEff1 + 1;
        end
    catch
        continue
    end
end
if runsEff1 > 0
    freq1 = countsSel1 / runsEff1; [~,ord1] = sort(freq1,'descend');
    cand1 = preds_LST(ord1(freq1(ord1)>0));
else
    cand1 = string.empty(0,1);
end
if isempty(cand1)
    rvec = corr(X1, y_LST, 'rows','pairwise'); [~,ordU] = sort(abs(rvec),'descend');
    nTop1 = min(topCap_stepwise, sum(isfinite(rvec))); cand1 = preds_LST(ordU(1:nTop1));
else
    nTop1 = min(topCap_stepwise, numel(cand1)); cand1 = cand1(1:nTop1);
end

tbl2 = tbl_LSTsrc; tbl2.y = y_LST;
mdlStep1 = stepwiselm(tbl2,'y ~ 1', 'Upper', sprintf('y ~ %s', strjoin(cellstr(cand1),' + ')), ...
    'PEnter',0,'PRemove',-Inf,'NSteps',nTop1,'Criterion','rsquared','Verbose',0);
sel_LST = setdiff(mdlStep1.PredictorNames, '(Intercept)', 'stable');

% Final model fit (always use y)
tblFit_LST = tbl_LSTsrc; tblFit_LST.y = y_LST;
mdl_LST = fitlm(tblFit_LST, sprintf('y ~ %s', strjoin(cellstr(sel_LST),' + ')));

R2f1    = mdl_LST.Rsquared.Ordinary;
dR2_1   = nan(numel(sel_LST),1);
for k=1:numel(sel_LST)
    keepk = sel_LST; keepk(k)=[];
    mdlk  = fitlm(tblFit_LST, sprintf('y ~ %s', strjoin(cellstr(keepk),' + ')));
    dR2_1(k) = R2f1 - mdlk.Rsquared.Ordinary;
end
imp_LST = (sum(dR2_1)>0) .* (dR2_1 / max(sum(dR2_1), eps));

%% ==========================
%% NDVI DV (top row) — drop NDVI-domain predictors
%% ==========================
[tbl_NDVIsrc, y_NDVI, preds_NDVI] = build_xy( ...
    tblZ, tblRAW, respNDVI, baseExclusions, false, 'NDVI', ...
    ndviExcludeNames, ndviExcludePatterns, lstExcludeNames, lstExcludePatterns);
if isempty(preds_NDVI), error('No usable predictors for %s after exclusions.', respNDVI); end
X2 = tbl_NDVIsrc{:, preds_NDVI}; n2 = height(tbl_NDVIsrc);

% One full-data CV fit (optional reference)
if n2 >= 3
    cvp_full = cvpartition(n2, 'KFold', max(2, min(5, floor(n2/2))));
    [B2, FI2] = lasso(X2, y_NDVI, 'Alpha',0.7, 'CV',cvp_full, 'Standardize',false, ...
        'NumLambda',150, 'LambdaRatio',1e-4, 'Options',optsPar);
    coef2 = B2(:, FI2.IndexMinMSE); selLasso2_full = preds_NDVI(coef2 ~= 0);
else
    selLasso2_full = string.empty(0,1);
end

% Stability selection for NDVI DV
nBoot = R_stab_NDVI; subRate = subFrac_NDVI;
cand2  = strings(0,1); runsEff2 = 0; selCounts2 = containers.Map('KeyType','char','ValueType','double');
for b = 1:nBoot
    rng(2000 + b, 'twister');
    nSub = max(5, round(subRate*n2));
    if nSub >= n2
        idx = (1:n2)';                     % no subsample if tiny n
    else
        idx = randsample(n2, nSub, false);
    end

    if numel(idx) < 3, continue; end
    try
        cvp_b = cvpartition(numel(idx), 'KFold', max(2, min(5, floor(numel(idx)/2))));
        [Bb, FIb] = lasso(X2(idx,:), y_NDVI(idx), 'Alpha',0.7, 'CV',cvp_b, 'Standardize',false, ...
            'NumLambda',120, 'LambdaRatio',1e-4, 'Options',optsPar);
        coefb = Bb(:, FIb.IndexMinMSE); selb = preds_NDVI(coefb ~= 0);
        runsEff2 = runsEff2 + 1;
        if ~isempty(selb)
            cand2 = [cand2; string(selb(:))]; %#ok<AGROW>
            for s = 1:numel(selb)
                key = char(selb{s});
                if isKey(selCounts2, key), selCounts2(key) = selCounts2(key) + 1;
                else, selCounts2(key) = 1; end
            end
        end
    catch
        continue
    end
end

if runsEff2 == 0 || isempty(cand2)
    stabTbl2 = table(string.empty(0,1), [], [], 'VariableNames', {'Predictor','Freq','Frac'});
else
    u = unique(cand2,'stable'); freq = zeros(numel(u),1);
    for i = 1:numel(u), key = char(u(i)); if isKey(selCounts2, key), freq(i) = selCounts2(key); end, end
    frac = freq / max(1, runsEff2);
    stabTbl2 = table(u, freq, frac, 'VariableNames', {'Predictor','Freq','Frac'});
    stabTbl2 = sortrows(stabTbl2, {'Frac','Predictor'}, {'descend','ascend'});
end
disp('[Stability] NDVI DV top candidates by frequency:'); disp(stabTbl2);

stableSet = string.empty(0,1);
if ~isempty(stabTbl2), stableSet = stabTbl2.Predictor( stabTbl2.Frac >= tauStable_NDVI ); end
if isempty(stableSet) && ~isempty(stabTbl2)
    K = min(12, height(stabTbl2)); stableSet = stabTbl2.Predictor(1:K);
end
stableSet = unique([stableSet; string(selLasso2_full(:))], 'stable');

if isempty(stableSet)
    r = corr(X2, y_NDVI, 'rows','pairwise'); [~,ord] = sort(abs(r),'descend');
    K = min(topCap_stepwise, sum(isfinite(r))); stableSet = string(preds_NDVI(ord(1:K)));
    fprintf('ENet selected none; fallback to top-%d univariate predictors.\n', K);
end

nTop2 = min(topCap_stepwise, numel(stableSet));
tbl2 = tbl_NDVIsrc; tbl2.y = y_NDVI;
upperForm2 = sprintf('y ~ %s', strjoin(cellstr(stableSet), ' + '));
mdlStep2 = stepwiselm(tbl2, 'y ~ 1', 'Upper', upperForm2, ...
    'PEnter',0, 'PRemove',-Inf, 'NSteps',nTop2, 'Criterion','rsquared', 'Verbose',0);
sel_NDVI = setdiff(mdlStep2.PredictorNames, '(Intercept)', 'stable');

% Final model fit (always use y)
tblFit_NDVI = tbl_NDVIsrc; tblFit_NDVI.y = y_NDVI;
mdl_NDVI = fitlm(tblFit_NDVI, sprintf('y ~ %s', strjoin(cellstr(sel_NDVI),' + ')));
R2f2     = mdl_NDVI.Rsquared.Ordinary;

dR2_2 = nan(numel(sel_NDVI),1);
for k = 1:numel(sel_NDVI)
    keepk = sel_NDVI; keepk(k) = [];
    mdlk  = fitlm(tblFit_NDVI, sprintf('y ~ %s', strjoin(cellstr(keepk),' + ')));
    dR2_2(k) = R2f2 - mdlk.Rsquared.Ordinary;
end
imp_NDVI = (sum(dR2_2)>0) .* (dR2_2 / max(sum(dR2_2), eps));

%% =========================
%% Visualization (2-row PDP layout)
%% =========================
nCols = 1 + max([numel(sel_LST), numel(sel_NDVI)]);
cmap  = parula(256);
figure('Color','w','Units','pixels','Position',[100 100 plotWidthPx plotHeightPx]);
tiledlayout(2, nCols, 'TileSpacing','compact','Padding','compact');

% ----- Row 1: NDVI DV -----
ax1 = nexttile(1);
v = y_NDVI; xLo = prctile(v,5); xHi = prctile(v,95);
histogram(ax1, v, 'FaceColor',[.7 .7 .7],'EdgeColor','none','Normalization','count');
xlim(ax1,[xLo xHi]); title(ax1, sprintf('%s distribution', strrep(respNDVI,'_','\_')), 'FontSize',9);
ylabel(ax1,'# Cities','FontSize',9); xlabel(ax1,strrep(respNDVI,'_','\_'),'FontSize',9);

mu1 = varfun(@(z) mean(z,'omitnan'), tbl_NDVIsrc(:, sel_NDVI)); mu1.Properties.VariableNames = sel_NDVI;
scores1 = imp_NDVI ./ max(imp_NDVI + (max(imp_NDVI)==0));
for j = 1:numel(sel_NDVI)
    ax = nexttile(1 + j); hold(ax,'on');
    vpred = tbl_NDVIsrc.(sel_NDVI{j});
    a = prctile(vpred,5); b = prctile(vpred,95);
    xGrid = linspace(a,b,120)'; 
    P = repmat(mu1{1,:}, numel(xGrid), 1); P = array2table(P,'VariableNames',sel_NDVI);
    P.(sel_NDVI{j}) = xGrid;
    [yHat, yCI] = predict(mdl_NDVI, P, 'Prediction','curve');
    cidx = min(max(1 + round(scores1(j)*255), 1), 256);
    plot(ax, xGrid, yHat, '-', 'LineWidth', 1.8, 'Color', cmap(cidx,:));
    fill(ax,[xGrid;flipud(xGrid)],[yCI(:,1);flipud(yCI(:,2))],cmap(cidx,:), ...
        'FaceAlpha',0.18,'EdgeColor','none');
    title(ax, sprintf('Imp=%.2f', scores1(j)), 'Interpreter','none','FontSize',9);
    xlabel(ax, sel_NDVI{j}, 'Interpreter','none','FontSize',9);
    if j==1, ylabel(ax, strrep(respNDVI,'_','\_'),'FontSize',9); end
    grid(ax,'off'); box(ax,'off'); set(ax,'FontSize',9);
    hold(ax,'off');
end

% ----- Row 2: LST DV -----
ax2 = nexttile(nCols + 1);
v = y_LST; xLo = prctile(v,5); xHi = prctile(v,95);
histogram(ax2, v, 'FaceColor',[.7 .7 .7],'EdgeColor','none','Normalization','count');
xlim(ax2,[xLo xHi]); title(ax2, sprintf('%s distribution', strrep(respLST,'_','\_')), 'FontSize',9);
ylabel(ax2,'# Cities','FontSize',9); xlabel(ax2,strrep(respLST,'_','\_'),'FontSize',9);

mu2 = varfun(@(z) mean(z,'omitnan'), tbl_LSTsrc(:, sel_LST)); mu2.Properties.VariableNames = sel_LST;
scores2 = imp_LST ./ max(imp_LST + (max(imp_LST)==0));
for j = 1:numel(sel_LST)
    ax = nexttile(nCols + 1 + j); hold(ax,'on');
    vpred = tbl_LSTsrc.(sel_LST{j});
    a = prctile(vpred,5); b = prctile(vpred,95);
    xGrid = linspace(a,b,120)';
    P = repmat(mu2{1,:}, numel(xGrid), 1); P = array2table(P,'VariableNames',sel_LST);
    P.(sel_LST{j}) = xGrid;
    [yHat, yCI] = predict(mdl_LST, P, 'Prediction','curve');
    cidx = min(max(1 + round(scores2(j)*255), 1), 256);
    plot(ax, xGrid, yHat, '-', 'LineWidth', 1.8, 'Color', cmap(cidx,:));
    fill(ax,[xGrid;flipud(xGrid)],[yCI(:,1);flipud(yCI(:,2))],cmap(cidx,:), ...
        'FaceAlpha',0.18,'EdgeColor','none');
    title(ax, sprintf('Imp=%.2f', scores2(j)), 'Interpreter','none','FontSize',9);
    xlabel(ax, sel_LST{j}, 'Interpreter','none','FontSize',9);
    if j==1, ylabel(ax, strrep(respLST,'_','\_'),'FontSize',9); end
    grid(ax,'off'); box(ax,'off'); set(ax,'FontSize',9);
    hold(ax,'off');
end

% Global cosmetics
axAll = findall(gcf,'Type','axes'); set(axAll,'XGrid','off','YGrid','off','ZGrid','off','XMinorGrid','off','YMinorGrid','off');
sgtitle(sprintf('Top predictors — %s (Z predictors, NON-Z DVs or CityTrends slopes) — %s', ...
        strrep(respNDVI,'_','\_'), tableTag), 'FontSize',12);

%% =========================
%% Fit summaries
%% =========================
mkRow = @(label, mdl, nPred) {label, mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted, mdl.RMSE, ...
                      mdl.NumObservations, nPred, strjoin(mdl.PredictorNames, ', ')};
FitSummary = cell2table( ...
   [mkRow(respLST,  mdl_LST,  numel(sel_LST)); ...
    mkRow(respNDVI, mdl_NDVI, numel(sel_NDVI))], ...
   'VariableNames', {'Response','R2','R2Adj','RMSE','N','NumPredictors','Predictors'});
disp(FitSummary);
fprintf('Adjusted R² — %s : %.3f (R²=%.3f, k=%d, n=%d)\n', ...
    respLST,  mdl_LST.Rsquared.Adjusted,  mdl_LST.Rsquared.Ordinary,  numel(sel_LST),  mdl_LST.NumObservations);
fprintf('Adjusted R² — %s: %.3f (R²=%.3f, k=%d, n=%d)\n', ...
    respNDVI, mdl_NDVI.Rsquared.Adjusted, mdl_NDVI.Rsquared.Ordinary, numel(sel_NDVI), mdl_NDVI.NumObservations);
