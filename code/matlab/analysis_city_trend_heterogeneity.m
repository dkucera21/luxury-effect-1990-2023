%% city_trend_heterogeneity_from_struct.m
% Tests whether city-specific trends (inc10k × yearC) vary across cities,
% running the full pipeline for ALL FOUR outcomes:
%   MEDIAN_NDVI, MEAN_NDVI, MEDIAN_LST, MEAN_LST
%
% Inputs expected in workspace:
%   - CENSUS_TABLES_rebuilt (struct of city-year tract tables)
%   - CityListMaster (string/cellstr of valid city names) or MasterCities/CityNames
%
% Outputs to base workspace:
%   - CityTrends_MEDIAN_NDVI, CityTrends_MEAN_NDVI,
%     CityTrends_MEDIAN_LST,  CityTrends_MEAN_LST   (each: City, Slope, SE, CI_L, CI_U, nYears)
%   - Back-compat shims:
%       CityTrends_NDVI = CityTrends_MEDIAN_NDVI
%       CityTrends_LST  = CityTrends_MEDIAN_LST
%
% Printed per outcome:
%   - Primary LMM (REML) inc10k×yearC fixed effect context (via var comps), BLUP spread, reliability
%   - (Guarded) ML LRT for adding random slope inc10k×yearC by City
%   - Meta-style heterogeneity of per-city LE trends (Q, I^2, tau^2) with robust guards


%% ---------------- USER SETTINGS ----------------
yearsUse     = [1990 2000 2010 2020 2023];
coveragePct  = 50;     % keep tracts with >= this cover (%)
minTracts    = 5;      % per city-year minimum valid tracts
useWinsor    = false;  % winsorize X and Y within city-year
winzPct      = [1 99]; % [5 95] auto if n<50
tryKR        = true;   % try Kenward–Roger DF (falls back silently)
% ------------------------------------------------

metrics = struct( ...
  'name', {'MEDIAN_NDVI','MEAN_NDVI','MEDIAN_LST','MEAN_LST'}, ...
  'pretty',{'Median NDVI','Mean NDVI','Median LST','Mean LST'} );

%% Guards and common lookups
assert(exist('CENSUS_TABLES_rebuilt','var')==1 && isstruct(CENSUS_TABLES_rebuilt), ...
    'CENSUS_TABLES_rebuilt is required and must be a struct.');

if exist('CityListMaster','var')==1 && ~isempty(CityListMaster)
    CityListMaster = string(CityListMaster(:));
elseif exist('MasterCities','var')==1 && istable(MasterCities) && ismember('City', MasterCities.Properties.VariableNames)
    CityListMaster = string(MasterCities.City);
elseif exist('CityNames','var')==1 && ~isempty(CityNames)
    CityListMaster = string(CityNames(:));
else
    error('No valid city list (CityListMaster / MasterCities.City / CityNames).');
end
CityListMaster = unique(strtrim(CityListMaster));
tblNames = string(fieldnames(CENSUS_TABLES_rebuilt));
numify   = @(x) str2double(regexprep(string(x),'[,\$%]',''));

pctNames    = {'PCT_OVERLAP','PCT_COVER','PCT_COV','PercentAreaNDVI'};
incRawCands = [{'RAW_INCOME_CPI'}, ...
    {'S1903_C03_001E','S1903_003_001E','S1903_C02_001E','DP3_C112','P080A001', ...
     'MEDIAN_INCOME','MED_INCOME','MED_INCOME_TOTAL','MEDIAN_HH_INCOME','MedianIncome','HCT036001'}];

try, if isempty(gcp('nocreate')), parpool; end, end

%% ---------- RUN LOOP OVER ALL FOUR OUTCOMES ----------
for m = 1:numel(metrics)
    metricChoice = metrics(m).name;
    fprintf('\n=== Heterogeneity of city-specific trends for %s ===\n', metricChoice);

    %% Build tall tract-level table (per outcome)
    nC = numel(CityListMaster);
    TallParts = cell(nC,1);

    parfor ci = 1:nC
        city = CityListMaster(ci);
        CityCol = {}; YearCol = []; Xcol = []; Ycol = [];

        for yr = yearsUse
            tname = sprintf('T_%d_%s', yr, city);
            if ~ismember(tname, tblNames), continue; end
            T = CENSUS_TABLES_rebuilt.(tname);
            if ~istable(T) || height(T)==0, continue; end

            V   = string(T.Properties.VariableNames);
            vPct = first_match(V, pctNames);        if vPct == "", continue; end
            vY   = first_match(V, {metricChoice});  if vY   == "", continue; end
            vInc = first_match(V, incRawCands);     if vInc == "", continue; end

            pct = numify(T.(vPct));
            if max(pct(isfinite(pct))) <= 1.5, pct = 100*pct; end
            x10 = numify(T.(vInc)) ./ 10000;   % per $10k
            yv  = numify(T.(vY));

            valid = isfinite(pct) & pct >= coveragePct & isfinite(x10) & isfinite(yv);
            if nnz(valid) < minTracts, continue; end

            x = x10(valid); y = yv(valid);
            if useWinsor
                lim = winzPct; if numel(x) < 50, lim = [5 95]; end
                x = local_winsor(x, lim); y = local_winsor(y, lim);
            end
            if std(x,'omitnan')==0, continue; end

            n = numel(x);
            CityCol(end+1:end+n,1) = repmat({char(city)}, n, 1); %#ok<AGROW>
            YearCol(end+1:end+n,1) = repmat(yr, n, 1);          %#ok<AGROW>
            Xcol(end+1:end+n,1)    = x;                          %#ok<AGROW>
            Ycol(end+1:end+n,1)    = y;                          %#ok<AGROW>
        end

        if ~isempty(Xcol)
            TallParts{ci} = table(categorical(CityCol), YearCol, Xcol, Ycol, ...
                'VariableNames', {'City','year','inc10k','y'});
        else
            TallParts{ci} = table('Size',[0 4], 'VariableTypes',{'categorical','double','double','double'}, ...
                'VariableNames', {'City','year','inc10k','y'});
        end
    end

    Tall = vertcat(TallParts{:});
    Tall.City = categorical(string(Tall.City), CityListMaster);
    Tall.City = removecats(Tall.City);
    Tall = Tall(isfinite(Tall.y) & isfinite(Tall.inc10k) & isfinite(Tall.year) & ~ismissing(Tall.City), :);
    if numel(unique(Tall.year)) < 2
        warning('%s: Too few unique years after filtering. Skipping.', metricChoice);
        continue
    end
    Tall.yearC = Tall.year - mean(Tall.year,'omitnan');

    %% Primary REML fit (for variance/BLUP reporting)
    dfm = 'Kenward-Roger'; if ~tryKR, dfm = 'Satterthwaite'; end
    formRE = 'y ~ inc10k*yearC + (1 + inc10k + yearC + inc10k:yearC | City)';
    mdlRE  = try_fitlme(Tall, formRE, dfm);  % REML
    varRS  = extract_rs_variance(mdlRE,'inc10k:yearC');

    % BLUP SD & reliability
    [blup_delta, blup_se] = get_city_rs_blups(mdlRE, 'inc10k:yearC');
    sdBLUP = std(blup_delta,'omitnan');
    rho = NaN;
    if ~isempty(blup_delta) && ~isempty(blup_se)
        varBLUP = sdBLUP.^2;
        meanSE2 = mean(blup_se.^2, 'omitnan');
        rho     = varBLUP / (varBLUP + meanSE2);
    end

    %% ML LRT for random slope (guarded)
    pLRT = NaN; LR = NaN; dDF = NaN;
    try
        form_full = 'y ~ inc10k*yearC + (1 + inc10k + yearC + inc10k:yearC | City)';
        form_red  = 'y ~ inc10k*yearC + (1 + inc10k + yearC | City)';
        mdl_full_ML = fitlme(Tall, form_full, 'FitMethod','ML');
        mdl_red_ML  = fitlme(Tall, form_red,  'FitMethod','ML');
        Cmp = compare(mdl_red_ML, mdl_full_ML);  % larger model second
        if istable(Cmp) && all(ismember({'LRStat','DeltaDF','pValue'}, Cmp.Properties.VariableNames))
            LR  = double(Cmp.LRStat(end));
            dDF = double(Cmp.DeltaDF(end));
            pLRT= double(Cmp.pValue(end));
        end
    catch
        % leave NaN
    end

    %% ===== City-year LE and meta-style heterogeneity =====
    % City-year LE: slope of y ~ inc10k within each City×Year
    G  = findgroups(Tall.City, Tall.year);
    LE = splitapply(@(x,y) local_slope(x,y), Tall.inc10k, Tall.y, G);
    C_first = splitapply(@(c) c(1), Tall.City, G);
    Y_first = splitapply(@(y) y(1),  Tall.year, G);
    TcityYear = table(C_first, Y_first, LE, 'VariableNames',{'City','year','LE'});
    TcityYear = TcityYear(isfinite(TcityYear.LE), :);
    TcityYear.yearC = TcityYear.year - mean(TcityYear.year);

    % Per-city OLS slope of LE ~ yearC (centered within city) and its SE
    uC = unique(string(TcityYear.City));
    k  = numel(uC);
    City   = strings(k,1);
    Slope  = nan(k,1);
    SE     = nan(k,1);
    CI_L   = nan(k,1);
    CI_U   = nan(k,1);
    nYears = nan(k,1);

    for i = 1:k
        Ti = TcityYear(string(TcityYear.City)==uC{i}, :);
        yrsU = unique(double(Ti.year));
        nYears(i) = numel(yrsU);
        City(i)   = string(uC{i});
        if numel(yrsU) < 3, continue; end
        yc = double(Ti.year) - mean(double(Ti.year),'omitnan');
        if std(yc)==0, continue; end
        yi = double(Ti.LE);

        X  = [ones(numel(yc),1) yc];
        XtX = X.'*X; if rcond(XtX) < 1e-12, continue; end
        b  = XtX \ (X.'*yi);
        r  = yi - X*b;
        dof= max(1, numel(yi)-2);
        s2 = sum(r.^2)/dof;
        Vb = s2 * inv(XtX);

        Slope(i) = b(2);
        SE(i)    = sqrt(max(eps, Vb(2,2)));
        tcrit    = tinv(0.975, dof);
        CI_L(i)  = Slope(i) - tcrit*SE(i);
        CI_U(i)  = Slope(i) + tcrit*SE(i);
    end

    % Keep valid rows
    kOK = isfinite(Slope) & isfinite(SE) & (SE > 0);
    CityTrends = table(City(kOK), Slope(kOK), SE(kOK), CI_L(kOK), CI_U(kOK), nYears(kOK), ...
        'VariableNames', {'City','Slope','SE','CI_L','CI_U','nYears'});

    % Assign per-outcome table to base with clear names
    outName = "CityTrends_" + string(metricChoice);
    assignin('base', outName, CityTrends);

    % Back-compat shims for MEDIAN_* only (once when relevant)
    if metricChoice == "MEDIAN_NDVI"
        assignin('base','CityTrends_NDVI', CityTrends);
    elseif metricChoice == "MEDIAN_LST"
        assignin('base','CityTrends_LST',  CityTrends);
    end

    % Meta heterogeneity (guarded)
    beta_c = CityTrends.Slope; se_c = CityTrends.SE;
    hasMeta = numel(beta_c) >= 3 && all(isfinite(se_c)) && any(se_c > 0);
    if hasMeta
        wFE     = 1 ./ (se_c.^2);
        beta_FE = sum(wFE .* beta_c) / sum(wFE);
        Q       = sum(wFE .* (beta_c - beta_FE).^2);
        dfQ     = numel(beta_c) - 1;
        pQ      = 1 - chi2cdf(max(Q,0), dfQ);
        I2      = max(0, (Q - dfQ) / max(Q, eps));
        cJ      = sum(wFE) - (sum(wFE.^2) / sum(wFE));
        tau2    = 0; if cJ > 0, tau2 = max(0, (Q - dfQ) / cJ); end
        wRE     = 1 ./ (se_c.^2 + tau2);
        beta_RE = sum(wRE .* beta_c) / sum(wRE);

        fprintf('\n[Meta-style heterogeneity of per-city LE trends]\n');
        fprintf('  Cochran Q = %.3f (df=%d), p = %.3g;  I^2 = %.1f%%%%;  tau^2 = %.6g\n', ...
            Q, dfQ, pQ, 100*I2, tau2);
        fprintf('  Mean trend: fixed-effects = %.6g, random-effects = %.6g (units per year)\n', ...
            beta_FE, beta_RE);
        fprintf('  Usable cities: %d (median years per city = %g)\n', numel(beta_c), median(CityTrends.nYears,'omitnan'));
        fprintf('  Median per-city slope = %.6g, IQR = [%.6g, %.6g]\n', ...
            median(beta_c,'omitnan'), quantile(beta_c,0.25), quantile(beta_c,0.75));
    else
        fprintf('\n[Meta-style heterogeneity of per-city LE trends]\n');
        fprintf('  Not enough cities with valid slope & SE (k=%d). Q/I^2/τ^2 not computed.\n', numel(beta_c));
    end

    %% Final report blocks (guarded where needed)
    fprintf('\n[Model comparison (ML)]\n');
    if isfinite(pLRT) && isfinite(LR) && isfinite(dDF)
        fprintf('  LRT for random city slope (inc10k×yearC): Δdf=%d, LR=%.3f, p=%g\n', dDF, LR, pLRT);
        if pLRT < 0.05
            fprintf('  ⇒ Evidence that trend differs across cities (heterogeneity present).\n');
        else
            fprintf('  ⇒ No strong evidence that trend differs across cities (homogeneous slopes plausible).\n');
        end
    else
        fprintf('  LRT unavailable in this run (model comparison failed or not supported).\n');
    end

    fprintf('\n[REML diagnostics]\n');
    fprintf('  Random-slope variance (inc10k×yearC): %.6g\n', varRS);
    if isfinite(sdBLUP)
        fprintf('  BLUP SD across cities (deltas around FE): %.6g\n', sdBLUP);
    end
    if isfinite(rho)
        fprintf('  Reliability rho of city trend BLUPs: %.3f (implies max explainable R^2 ≈ %.3f)\n', ...
            rho, max(0,min(1,rho)));
    end

    fprintf('\n— Done (%s) —\n', metrics(m).pretty);
end

%% -------- local helpers (script-scoped) --------
function v = first_match(V, cands)
    v = "";
    VL = lower(string(V));
    for i = 1:numel(cands)
        j = find(VL == lower(string(cands{i})), 1, 'first');
        if ~isempty(j), v = string(V(j)); return; end
    end
end

function xw = local_winsor(x, lim)
    lo = prctile(x, lim(1)); hi = prctile(x, lim(2));
    xw = min(max(x, lo), hi);
end

function mdl = try_fitlme(ds, form, dfmeth)
    try
        mdl = fitlme(ds, form, 'FitMethod','REML', 'DFMethod', dfmeth);
    catch
        mdl = fitlme(ds, form, 'FitMethod','REML');
    end
end

function varRS = extract_rs_variance(mdl, term)
    varRS = NaN;
    try
        CP = covarianceParameters(mdl);
        M = [];
        if iscell(CP) && ~isempty(CP) && isnumeric(CP{1})
            M = CP{1};
        elseif iscell(CP) && ~isempty(CP) && iscell(CP{1}) && ~isempty(CP{1}) && isnumeric(CP{1}{1})
            M = CP{1}{1};
        end
        if ~isempty(M)
            names = {'(Intercept)','inc10k','yearC','inc10k:yearC'};
            idx = find(strcmpi(names, term), 1);
            if ~isempty(idx) && idx <= size(M,1)
                varRS = M(idx,idx);
            end
        end
    catch
        % leave NaN
    end
end

function [delta, se] = get_city_rs_blups(mdl, term)
    delta = []; se = [];
    try
        [b, names, stats] = randomEffects(mdl);
        Name = []; SE = [];
        if istable(names) && any(strcmpi('Name', names.Properties.VariableNames))
            Name = string(names.Name);
        end
        if isstruct(stats) && isfield(stats,'SE')
            SE = stats.SE(:);
        elseif istable(stats) && any(strcmpi('SE', stats.Properties.VariableNames))
            SE = double(stats.SE);
        else
            SE = NaN(size(b(:)));
        end
        want = contains(lower(Name), lower(term));
        delta = b(want); se = SE(want);
        if ~isempty(delta), return; end
    catch
    end
    try
        RE_tbl = randomEffects(mdl);
        v = string(RE_tbl.Properties.VariableNames);
        colName  = find(ismember(lower(v), lower(["Name","Effect","EffName","REName"])),1);
        colEst   = find(ismember(lower(v), lower(["Estimate","RE","Value","Est"])),1);
        if isempty(colEst) || isempty(colName), return; end
        Name = string(RE_tbl.(v(colName)));
        Est  = double(RE_tbl.(v(colEst)));
        if any(strcmpi(v,'SE')), SE = double(RE_tbl.('SE')); else, SE = NaN(size(Est)); end
        want = contains(lower(Name), lower(term));
        delta = Est(want); se = SE(want);
    catch
    end
end

function b = local_slope(x, y)
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    if numel(x) < 3 || std(x)==0, b = NaN; return; end
    X = [ones(numel(x),1) x];
    beta = X\y;
    b = beta(2);
end
