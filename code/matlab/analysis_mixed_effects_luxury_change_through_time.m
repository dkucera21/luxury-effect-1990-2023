%% ========== Script: Mixed-Effects Luxury-Effect Change (per $10k per year) ==========
% Fits (primary, tract-level): y ~ (RAW_INCOME_CPI/10k) * centered(year)
% Random effects by City (rich structure with fallbacks).
% Prints: fixed effect (per $10k/yr), 95% CI, p (REML + Satterthwaite/KR*),
% variance components, Nakagawa R^2, practical translations, LOCO,
% two-stage confirmation, AND diagnostics to assess whether between-city
% trend heterogeneity is real or mostly noise:
%  - LRT for City random slope on inc10k:yearC (refit with ML for LRT)
%  - Random-slope variance (if available)
%  - BLUPs & SEs for inc10k:yearC by City
%  - Reliability rho and implied max explainable cross-city R^2
% *If DF methods aren't supported, falls back silently.

% ---------------------- USER SETTINGS ----------------------
yearsUse  = [1990 2000 2010 2020 2023];
thresh    = 50;      % PCT_OVERLAP cutoff (>=)
minTracts = 5;       % minimum valid tracts per city-year
useWinsor = false;   % light winsorization for tract-level x,y
winzPct   = [1 99];  % [5 95] auto for n<50
tryKR     = true;    % try Kenward–Roger df if available
% -----------------------------------------------------------

try, if isempty(gcp('nocreate')), parpool; end, end

metrics = struct( ...
  'y',   {'MEDIAN_NDVI','MEAN_NDVI','MEDIAN_LST','MEAN_LST'}, ...
  'tag', {'MEDIAN_NDVI','MEAN_NDVI','MEDIAN_LST','MEAN_LST'} );

%% ---------- inputs / guards ----------
assert(exist('CENSUS_TABLES_rebuilt','var')==1 && isstruct(CENSUS_TABLES_rebuilt), ...
    'CENSUS_TABLES_rebuilt not found or not a struct.');

% Master city list
if exist('CityListMaster','var') && ~isempty(CityListMaster)
    CityListMaster = string(CityListMaster(:));
elseif exist('MasterCities','var') && istable(MasterCities) && ismember('City', MasterCities.Properties.VariableNames)
    CityListMaster = string(MasterCities.City);
elseif exist('CityNames','var') && ~isempty(CityNames)
    warning('CityListMaster not found; falling back to CityNames (all cities).');
    CityListMaster = string(CityNames(:));
else
    error('No valid city list found (CityListMaster, MasterCities.City, or CityNames).');
end
CityListMaster = unique(strtrim(CityListMaster));
nC = numel(CityListMaster);

tblNames = string(fieldnames(CENSUS_TABLES_rebuilt));
numify   = @(x) str2double(regexprep(string(x),'[,\$%]',''));

% Coverage / income candidates
pctNames    = {'PCT_OVERLAP','PCT_COVER','PCT_COV','PercentAreaNDVI'};
incRawCands = [{'RAW_INCOME_CPI'}, ... % preferred (CPI-deflated)
    {'S1903_C03_001E','S1903_003_001E','S1903_C02_001E','DP3_C112','P080A001', ...
     'MEDIAN_INCOME','MED_INCOME','MED_INCOME_TOTAL','MEDIAN_HH_INCOME','MedianIncome','HCT036001'}];

fprintf('\n— Mixed-effects (global) change in luxury effect per $10k per YEAR —\n');

% Ensure trends tables exist / aligned
if ~exist('A_HYPO_TRENDS','var') || ~istable(A_HYPO_TRENDS)
    A_HYPO_TRENDS = table(CityListMaster, 'VariableNames', {'City'});
else
    A_HYPO_TRENDS.City = string(A_HYPO_TRENDS.City);
    A_HYPO_TRENDS = outerjoin(table(CityListMaster,'VariableNames',{'City'}), ...
                              A_HYPO_TRENDS, 'Keys','City','MergeKeys',true,'Type','left');
end
if ~exist('A_HYPO_TRENDS_pvals','var') || ~istable(A_HYPO_TRENDS_pvals)
    A_HYPO_TRENDS_pvals = table(CityListMaster, 'VariableNames', {'City'});
else
    A_HYPO_TRENDS_pvals.City = string(A_HYPO_TRENDS_pvals.City);
    A_HYPO_TRENDS_pvals = outerjoin(table(CityListMaster,'VariableNames',{'City'}), ...
                                    A_HYPO_TRENDS_pvals, 'Keys','City','MergeKeys',true,'Type','left');
end

% ---------- main loop over metrics ----------
for m = 1:numel(metrics)
    yName = metrics(m).y;
    tag   = metrics(m).tag;

    CityCol = {}; YearCol = []; IncCol = []; YCol = [];

    % ------- build tract-level table Tall (PARALLEL CITY LOOP) -------
    TallParts = cell(nC,1);
    parfor ci = 1:nC
        city = CityListMaster(ci);
        CityCol = {}; YearCol = []; IncCol = []; YCol = [];

        for yr = yearsUse
            tname = sprintf('T_%d_%s', yr, city);
            if ~ismember(tname, tblNames), continue; end

            T = CENSUS_TABLES_rebuilt.(tname);
            if ~istable(T) || height(T)==0, continue; end
            V = string(T.Properties.VariableNames);

            % coverage + metric + income
            vPct   = first_match(V, pctNames);      if vPct   == "", continue; end
            vY     = first_match(V, {yName});       if vY     == "", continue; end
            vInc   = first_match(V, incRawCands);   if vInc   == "", continue; end

            pct = numify(T.(vPct));
            pct(pct > 100) = 100;
            mxP = max(pct(isfinite(pct)));
            if isfinite(mxP) && mxP <= 1.5, pct = 100*pct; end

            x10 = numify(T.(vInc)) ./ 10000;        % per $10k
            yv  = numify(T.(vY));

            valid = isfinite(pct) & pct >= thresh & isfinite(x10) & isfinite(yv);
            if nnz(valid) < minTracts, continue; end

            x = x10(valid); y = yv(valid);
            if useWinsor
                lim = winzPct; if numel(x) < 50, lim = [5 95]; end
                x = local_winsor(x, lim); y = local_winsor(y, lim);
            end
            if std(x,'omitnan')==0, continue; end

            n = numel(x);
            CityCol(end+1:end+n,1) = repmat({char(city)}, n, 1); %#ok<AGROW>
            YearCol(end+1:end+n,1) = repmat(yr,   n, 1);        %#ok<AGROW>
            IncCol(end+1:end+n,1)  = x;                         %#ok<AGROW>
            YCol(end+1:end+n,1)    = y;                         %#ok<AGROW>
        end

        if ~isempty(IncCol)
            TallParts{ci} = table(categorical(CityCol), YearCol, IncCol, YCol, ...
                'VariableNames', {'City','year','inc10k','y'});
        else
            TallParts{ci} = table('Size',[0 4], 'VariableTypes',{'categorical','double','double','double'}, ...
                'VariableNames', {'City','year','inc10k','y'});
        end
    end

    Tall = vertcat(TallParts{:});

    % Clean & keep only valid rows
    Tall.City = categorical(string(Tall.City), CityListMaster);
    Tall.City = removecats(Tall.City);
    rOK = isfinite(Tall.y) & isfinite(Tall.inc10k) & isfinite(Tall.year) & ~ismissing(Tall.City);
    Tall = Tall(rOK, :);

    % Need ≥2 unique years
    uYears = unique(Tall.year);
    if numel(uYears) < 2
        warning('%s: only one year present after filtering (%s). Cannot estimate change over time.', ...
            tag, strjoin(string(uYears),','));
        write_nan_results(tag);
        continue
    end

    % Center year (global centering so interaction = change per calendar year)
    Tall.yearC = Tall.year - mean(Tall.year,'omitnan');

    % Need some rows and >1 city
    if height(Tall) < 50 || numel(categories(Tall.City)) < 2
        warning('%s: too few rows/cities after filtering (n=%d, nCities=%d).', ...
            tag, height(Tall), numel(categories(Tall.City)));
        write_nan_results(tag);
        continue
    end

    % ---------- Practical contrast: within-city income IQR ----------
    IQRs = varfun(@iqr, Tall, 'InputVariables','inc10k', ...
                  'GroupingVariables',{'City','year'});
    vn = IQRs.Properties.VariableNames;
    colI = find(contains(vn,'inc10k','IgnoreCase',true),1,'last');
    incIQR_mean = mean(IQRs.(vn{colI}), 'omitnan');

    % spans for reporting (years)
    span_1990_2023 = (2023 - 1990);
    span_2000_2023 = (2023 - 2000);

    % ---------- PRIMARY LMM (tract-level; REML) ----------
    % Rich RE structure, fallback ladder
    reTry = { ...
        '(1 + inc10k + yearC + inc10k:yearC | City)', ...
        '(1 + inc10k + yearC | City)', ...
        '(1 + inc10k | City)' };

    mdl   = [];
    usedForm = '';
    usedDF   = 'Satterthwaite'; if tryKR, usedDF = 'Kenward-Roger'; end

    for r = 1:numel(reTry)
        form = sprintf('y ~ inc10k*yearC + %s', reTry{r});
        try
            mdl = fitlme_df(Tall, form, usedDF);     % REML
            usedForm = form; break
        catch
            if r == 1 && tryKR
                try
                    mdl = fitlme_df(Tall, form, 'Satterthwaite');
                    usedForm = form; usedDF = 'Satterthwaite'; break
                catch, end
            end
        end
    end
    if isempty(mdl)
        usedForm = 'y ~ inc10k*yearC + (1 | City)';
        mdl = fitlme_df(Tall, usedForm, 'Satterthwaite');
        usedDF = 'Satterthwaite';
    end

    % Locate interaction
    cn = string(mdl.CoefficientNames);
    isInt = (cn=="inc10k:yearC") | (cn=="yearC:inc10k") | (contains(cn,"inc10k") & contains(cn,"yearC") & contains(cn,":"));
    k = find(isInt,1);
    if isempty(k)
        beta    = NaN; pval = NaN; ci_beta = [NaN NaN];
    else
        beta    = mdl.Coefficients.Estimate(k);           % (y per $10k) per year
        pval    = mdl.Coefficients.pValue(k);
        ciFull  = coefCI(mdl);  ci_beta = ciFull(k,:);
    end

    % Variance components + residual variance
    [vcTbl, sig2e] = safe_varcomps(mdl);

    % Marginal / conditional R^2 (Nakagawa)
    [R2m, R2c] = nakagawa_r2(mdl);

    % Practical translations
    beta_decade          = beta * 10;
    beta_decade_IQR      = beta * 10 * incIQR_mean;
    delta_1990_2023_IQR  = beta * span_1990_2023 * incIQR_mean;
    delta_2000_2023_IQR  = beta * span_2000_2023 * incIQR_mean;

    % ---------- Random-slope heterogeneity test (ML LRT) ----------
    % Compare models differing in random-effects → must use ML, not REML.
    mdl_full_ML = []; mdl_red_ML = []; pLRT = NaN; varRS = NaN;
    fullRS = 'y ~ inc10k*yearC + (1 + inc10k + yearC + inc10k:yearC | City)';
    redRS  = 'y ~ inc10k*yearC + (1 + inc10k + yearC | City)';
    try
        mdl_full_ML = fitlme(Tall, fullRS, 'FitMethod','ML');  %#ok<*NASGU>
        mdl_red_ML  = fitlme(Tall, redRS,  'FitMethod','ML');
        cmp = compare(mdl_red_ML, mdl_full_ML); % larger model second
        if istable(cmp) && any(strcmpi(cmp.Properties.VariableNames,'pValue'))
            pLRT = cmp.pValue(end);
        end
    catch
        % leave NaN if comparison fails
    end

    % Extract random-slope variance for inc10k:yearC if present (from REML fit)
    try
        varRS = extract_rs_variance(mdl, 'inc10k:yearC');
    catch, end

    % ---------- BLUPs & reliability for city trend heterogeneity ----------
    rho = NaN; noiseShare = NaN; meanSE2 = NaN; sdBLUP = NaN; nBLUP = 0;
    try
        [blup_delta, blup_se] = get_city_rs_blups(mdl, 'inc10k:yearC');  % deltas around FE beta
        nBLUP = numel(blup_delta);
        if nBLUP > 1
            sdBLUP   = std(blup_delta, 'omitnan');
            meanSE2  = mean(blup_se.^2, 'omitnan');
            varBLUP  = sdBLUP.^2;
            rho      = varBLUP / (varBLUP + meanSE2);
            noiseShare = 1 - rho;
        end
    catch
        % keep NaNs
    end

    % ---------- LOCO (Leave-One-City-Out) ----------
    uC   = categories(Tall.City);
    nUC  = numel(uC);
    loco = NaN(nUC,1);

    parfor i = 1:nUC
        ci_excl = categorical(uC(i), uC);
        Tall_i  = Tall(Tall.City ~= ci_excl, :);
        if height(Tall_i) < 50 || numel(categories(Tall_i.City)) < 2
            loco(i) = NaN; 
            continue
        end
        % Try same structure; fall back if needed
        try
            mdl_i = fitlme_df(Tall_i, usedForm, usedDF);
        catch
            try
                mdl_i = fitlme_df(Tall_i, 'y ~ inc10k*yearC + (1 + inc10k | City)', usedDF);
            catch
                mdl_i = fitlme_df(Tall_i, 'y ~ inc10k*yearC + (1 | City)', usedDF);
            end
        end
        cni = string(mdl_i.CoefficientNames);
        ki  = find((cni=="inc10k:yearC") | (cni=="yearC:inc10k") | (contains(cni,"inc10k") & contains(cni,"yearC") & contains(cni,":")),1);
        if ~isempty(ki), loco(i) = mdl_i.Coefficients.Estimate(ki); end
    end

    loco_med = median(loco,'omitnan');
    loco_iqr = iqr(loco(~isnan(loco)));

    % ---------- TWO-STAGE (city-year) confirmation ----------
    G  = findgroups(Tall.City, Tall.year);
    LE = splitapply(@(x,y) local_slope(x,y), Tall.inc10k, Tall.y, G);
    C_first = splitapply(@(c) c(1), Tall.City, G);
    Y_first = splitapply(@(y) y(1),  Tall.year, G);
    TcityYear = table(C_first, Y_first, LE, ...
        'VariableNames',{'City','year','LE'});
    TcityYear = TcityYear(isfinite(TcityYear.LE), :);
    TcityYear.yearC = TcityYear.year - mean(TcityYear.year);
    mdl2 = fitlme_df(TcityYear, 'LE ~ yearC + (1 | City)', usedDF);
    cn2   = string(mdl2.CoefficientNames);
    j2    = find(cn2=="yearC",1);
    b2    = mdl2.Coefficients.Estimate(j2);
    ci2   = coefCI(mdl2); ci2 = ci2(j2,:);
    p2    = mdl2.Coefficients.pValue(j2);
    [R2m2, R2c2] = nakagawa_r2(mdl2);
	
	% ---------- META-ANALYSIS style heterogeneity across cities (no LRT) ----------
	% Build per-city LE trend (slope of LE ~ yearC) + its SE, then compute Q, I2, tau2.
	% Uses the same TcityYear you just built for mdl2.

	% 1) Per-city slope of LE ~ yearC (with SE)
	uCities = categories(TcityYear.City);
	k = numel(uCities);
	beta_c = nan(k,1); se_c = nan(k,1); nObs_c = nan(k,1);

	for i = 1:k
		Ci = uCities{i};
		Ti = TcityYear(TcityYear.City == Ci, :);
		if height(Ti) < 3 || std(Ti.yearC)==0
			continue
		end
		Xi = [ones(height(Ti),1) double(Ti.yearC)];
		yi = double(Ti.LE);
		% OLS slope + robust SE
		betaHat = Xi \ yi;
		resid   = yi - Xi*betaHat;
		s2      = sum(resid.^2) / (height(Ti) - 2);
		covB    = s2 * inv(Xi.'*Xi);
		beta_c(i) = betaHat(2);
		se_c(i)   = sqrt(covB(2,2));
		nObs_c(i) = height(Ti);
	end

	ok = isfinite(beta_c) & isfinite(se_c) & se_c>0;
	beta_c = beta_c(ok); se_c = se_c(ok);
	k_eff  = numel(beta_c);

	if k_eff >= 3
		w_fixed = 1 ./ (se_c.^2);
		beta_FE = sum(w_fixed .* beta_c) / sum(w_fixed);

		% Cochran's Q
		Q = sum(w_fixed .* (beta_c - beta_FE).^2);
		dfQ = k_eff - 1;
		p_Q = 1 - chi2cdf(Q, dfQ);

		% I^2 (proportion of total variance due to between-city heterogeneity)
		I2 = max(0, (Q - dfQ) / max(Q, eps));

		% DerSimonian–Laird tau^2 (between-city variance of slopes)
		cJ   = sum(w_fixed) - (sum(w_fixed.^2) / sum(w_fixed));
		tau2 = max(0, (Q - dfQ) / max(cJ, eps));

		% Random-effects mean (for completeness)
		w_RE   = 1 ./ (se_c.^2 + tau2);
		beta_RE = sum(w_RE .* beta_c) / sum(w_RE);

		fprintf('  META heterogeneity (per-city trend of LE): k=%d\n', k_eff);
		fprintf('    Q=%.3f (df=%d), p=%.3g;  I^2=%.1f%%%%;  tau^2=%.6g\n', Q, dfQ, p_Q, 100*I2, tau2);
		fprintf('    Mean trend (FE)=%.6f, (RE)=%.6f (units per year)\n', beta_FE, beta_RE);
	else
		fprintf('  META heterogeneity: insufficient cities with valid per-city trend SEs (k=%d).\n', k_eff);
	end

	% ---------- BLUP spread sanity check (already printed above) ----------
	% If random-slope variance ~0 and I^2 ~0, the between-city trend differences are negligible.


    % ---------- PRINT RESULTS ----------
    ptxt = 'NA'; if isfinite(pval), ptxt = num2str(pval,'%.3g'); end
    fprintf('\n[%s] PRIMARY tract-level LMM\n', tag);
    fprintf('  Model form    : %s\n', usedForm);
    fprintf('  Fit/DF method : %s / %s\n', upper(mdl.FitMethod), safe_dfmethod(mdl));
    fprintf('  Fixed effect  : d(LE)/yr = %+0.4f per $10k/yr  [95%% CI %+0.4f, %+0.4f],  p = %s\n', ...
        beta, ci_beta(1), ci_beta(2), ptxt);
    fprintf('  Variance comps: %s\n', varcomp_str(vcTbl, sig2e));
    fprintf('  R^2 (Nakagawa): marginal = %.3f, conditional = %.3f\n', R2m, R2c);
    fprintf('  Income IQR    : mean within-city IQR(inc10k) = %.3f (≈ $%.0f)\n', incIQR_mean, 10000*incIQR_mean);
    fprintf('  Translation   : per decade = %+0.4f; per decade × IQR = %+0.4f\n', beta_decade, beta_decade_IQR);
    fprintf('                  1990→2023 × IQR = %+0.4f; 2000→2023 × IQR = %+0.4f\n', ...
        delta_1990_2023_IQR, delta_2000_2023_IQR);
    fprintf('  LOCO          : median d(LE)/yr = %+0.4f; IQR = %.4f\n', loco_med, loco_iqr);

    fprintf('  TWO-STAGE city-year LMM (LE ~ yearC + (1|City)):\n');
    fprintf('    d(LE)/yr = %+0.4f  [95%% CI %+0.4f, %+0.4f],  p = %.3g;  R^2_m=%.3f, R^2_c=%.3f\n', ...
        b2, ci2(1), ci2(2), p2, R2m2, R2c2);

    fprintf('  HETEROGENEITY of city trend differences (inc10k:yearC random slope):\n');
    if ~isnan(pLRT)
        fprintf('    LRT (with RS vs without RS, ML): p = %.3g  %s\n', pLRT, tern(pLRT<0.05,'[evidence of heterogeneity]','[no strong evidence]'));
    else
        fprintf('    LRT: unavailable (model comparison failed)\n');
    end
    if isfinite(varRS)
        fprintf('    Random-slope variance (City on inc10k:yearC): %.6f\n', varRS);
    else
        fprintf('    Random-slope variance: unavailable\n');
    end
    if nBLUP > 1 && isfinite(sdBLUP)
        fprintf('    BLUP SD across cities (deltas around FE): %.6f (n=%d)\n', sdBLUP, nBLUP);
    else
        fprintf('    BLUP SD: unavailable\n');
    end
    if isfinite(rho)
        fprintf('    Reliability rho of city trend estimates: %.3f  (noise share ≈ %.3f)\n', rho, noiseShare);
        fprintf('    Implied **max explainable cross-city R^2** (upper bound): ≈ %.3f\n', max(0,min(1,rho)));
    else
        fprintf('    Reliability rho: unavailable (no SEs/BLUPs)\n');
    end

    % store a single global value per metric (primary LMM fixed effect)
    colB = ['dLE_perYear_LME_' tag '_per10k'];
    colP = ['p_LE_perYear_LME_' tag];
    if ~ismember(colB, A_HYPO_TRENDS.Properties.VariableNames),       A_HYPO_TRENDS.(colB) = NaN(height(A_HYPO_TRENDS),1); end
    if ~ismember(colP, A_HYPO_TRENDS_pvals.Properties.VariableNames), A_HYPO_TRENDS_pvals.(colP) = NaN(height(A_HYPO_TRENDS_pvals),1); end
    A_HYPO_TRENDS.(colB)(:)       = beta;
    A_HYPO_TRENDS_pvals.(colP)(:) = pval;

    lbl = sprintf('[%s] primary LMM', tag);
	if exist('print_lme','file')==2
		% Try user's/other print_lme if it exists and matches this signature
		try
			print_lme(mdl, lbl);
		catch
			local_print_lme(mdl, lbl); % fall back to our local one
		end
	else
		local_print_lme(mdl, lbl);     % no external function → use local one
	end

end

%% ---------------------- helpers ----------------------
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

function varRS = extract_rs_variance(mdl, term)
% Try to find the variance for the requested random-slope 'term' in Psi
    varRS = NaN;
    try
        [Psi, ~] = covarianceParameters(mdl);
        M = [];
        if iscell(Psi) && ~isempty(Psi) && isnumeric(Psi{1})
            M = Psi{1};
        elseif iscell(Psi) && ~isempty(Psi) && iscell(Psi{1}) && ~isempty(Psi{1}) && isnumeric(Psi{1}{1})
            M = Psi{1}{1};
        end
        if ~isempty(M)
            names = [];
            try, names = mdl.RandomInfo(1).Name; catch, end
            if isempty(names) || numel(names) ~= size(M,1)
                names = {'(Intercept)','inc10k','yearC','inc10k:yearC'};
                names = names(1:min(end,size(M,1)));
            end
            idx = find(strcmpi(string(names), string(term)),1);
            if ~isempty(idx), varRS = M(idx,idx); end
        end
    catch
        % leave NaN
    end
end

function [delta, se] = get_city_rs_blups(mdl, term)
% Robustly return per-City BLUP deltas for the requested random-slope 'term'
% and their SEs if available (handles older/newer MATLAB).
    delta = []; se = [];

    % Attempt the 3-output signature first (most robust for SEs)
    try
        [b, names, stats] = randomEffects(mdl);
        Name  = []; Level = []; SE = [];
        if istable(names)
            if any(strcmpi('Name', names.Properties.VariableNames))
                Name = string(names.Name);
            elseif any(strcmpi('Effect', names.Properties.VariableNames))
                Name = string(names.Effect);
            end
            if any(strcmpi('Level', names.Properties.VariableNames))
                Level = string(names.Level);
            elseif any(strcmpi('Group', names.Properties.VariableNames))
                Level = string(names.Group);
            end
        end
        if isstruct(stats) && isfield(stats,'SE')
            SE = stats.SE(:);
        elseif istable(stats) && any(strcmpi('SE', stats.Properties.VariableNames))
            SE = double(stats.SE);
        else
            SE = NaN(size(b(:)));
        end
        if ~isempty(Name)
            want = contains(lower(Name), lower(term));
            delta = b(want);
            se    = SE(want);
            if ~isempty(delta), return; end
        end
    catch
        % fall through
    end

    % Fallback: single-table return
    try
        RE_tbl = randomEffects(mdl);  % table with columns like Name/Level/Estimate/SE
        v = string(RE_tbl.Properties.VariableNames);
        % Harmonize likely column names
        colName  = find(ismember(lower(v), lower(["Name","Effect","EffName","REName"])),1);
        colEst   = find(ismember(lower(v), lower(["Estimate","RE","Value","Est"])),1);
        if isempty(colEst) || isempty(colName), return; end
        Name = string(RE_tbl.(v(colName)));
        Est  = double(RE_tbl.(v(colEst)));
        if any(strcmpi(v,'SE'))
            SE = double(RE_tbl.('SE'));
        else
            SE = NaN(size(Est));
        end
        want = contains(lower(Name), lower(term));
        delta = Est(want);
        se    = SE(want);
    catch
        % leave empty
    end
end

function mdl = fitlme_df(ds, form, dfmeth)
    try
        mdl = fitlme(ds, form, 'FitMethod','REML', 'DFMethod', dfmeth);
    catch
        mdl = fitlme(ds, form, 'FitMethod','REML');
    end
end

function [R2m, R2c] = nakagawa_r2(mdl)
    yhatF = fitted(mdl,'Conditional',false);
    varF  = var(yhatF, 1, 'omitnan');
    [vcTbl, sig2e] = safe_varcomps(mdl);
    varRand = nansum(vcTbl.Var);
    varRes  = sig2e;
    denom   = varF + varRand + varRes;
    R2m = varF   / denom;
    R2c = (varF + varRand) / denom;
    if ~isfinite(R2m), R2m = NaN; end
    if ~isfinite(R2c), R2c = NaN; end
end

function [vcTbl, sig2e] = safe_varcomps(mdl)
    vcTbl = table(); sig2e = NaN;
    try, sig2e = mdl.MSE; catch, try, sig2e = mdl.Sigma^2; catch, end, end
    try
        [Psi, ~] = covarianceParameters(mdl);
        M = [];
        if iscell(Psi) && ~isempty(Psi) && isnumeric(Psi{1})
            M = Psi{1};
        elseif iscell(Psi) && ~isempty(Psi) && iscell(Psi{1}) && ~isempty(Psi{1}) && isnumeric(Psi{1}{1})
            M = Psi{1}{1};
        end
        if ~isempty(M)
            names = [];
            try, names = mdl.RandomInfo(1).Name; catch, end
            if isempty(names) || numel(names) ~= numel(diag(M))
                names = arrayfun(@(k) sprintf('RE%d',k), 1:numel(diag(M)), 'uni', 0);
            end
            vcTbl = table(string(names(:)), diag(M), 'VariableNames', {'Term','Var'});
            return
        end
        CP = covarianceParameters(mdl);
        if iscell(CP) && ~isempty(CP) && iscell(CP{1}) && ~isempty(CP{1}) && isnumeric(CP{1}{1})
            M = CP{1}{1};
            names = [];
            try, names = mdl.RandomInfo(1).Name; catch, end
            if isempty(names) || numel(names) ~= numel(diag(M))
                names = arrayfun(@(k) sprintf('RE%d',k), 1:numel(diag(M)), 'uni', 0);
            end
            vcTbl = table(string(names(:)), diag(M), 'VariableNames', {'Term','Var'});
            return
        end
        [~, names, stats] = randomEffects(mdl);
        lab = string(names.Name);
        vc  = accumarray(double(grp2idx(lab)), stats.SE.^2, [], @mean, NaN);
        u   = unique(lab,'stable');
        vcTbl = table(u, vc, 'VariableNames', {'Term','Var'});
    catch
        vcTbl = table();
    end
end

function s = varcomp_str(vcTbl, sig2e)
    parts = "ε=" + fmt(sig2e);
    if exist('vcTbl','var') && istable(vcTbl) && ~isempty(vcTbl) ...
            && any(strcmpi(vcTbl.Properties.VariableNames,'Term')) ...
            && any(strcmpi(vcTbl.Properties.VariableNames,'Var'))
        for i = 1:height(vcTbl)
            parts = parts + ", " + string(vcTbl.Term(i)) + "=" + fmt(vcTbl.Var(i));
        end
    end
    s = char(parts);
    function t = fmt(x), if ~isfinite(x), t='NA'; else, t = num2str(x,'%.4f'); end, end
end

function b = local_slope(x, y)
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    if numel(x) < 3 || std(x)==0, b = NaN; return; end
    X = [ones(numel(x),1) x];
    beta = X\y;
    b = beta(2);
end

function names = local_varnames(T)
    if istable(T)
        names = T.Properties.VariableNames;
    elseif isa(T,'dataset')
        names = get(T,'VarNames');
    else
        names = {};
    end
end

function write_nan_results(tag)
    if ~exist('A_HYPO_TRENDS','var') || ~exist('A_HYPO_TRENDS_pvals','var'), return; end
    colB = ['dLE_perYear_LME_' tag '_per10k'];
    colP = ['p_LE_perYear_LME_' tag];
    if ~ismember(colB, A_HYPO_TRENDS.Properties.VariableNames),       A_HYPO_TRENDS.(colB) = NaN(height(A_HYPO_TRENDS),1); end
    if ~ismember(colP, A_HYPO_TRENDS_pvals.Properties.VariableNames), A_HYPO_TRENDS_pvals.(colP) = NaN(height(A_HYPO_TRENDS_pvals),1); end
    A_HYPO_TRENDS.(colB)(:)       = NaN;
    A_HYPO_TRENDS_pvals.(colP)(:) = NaN;
end

function mdlOut = try_fit(mdlIn, Tall, form, dfmeth)
    mdlOut = mdlIn;
    try
        if ~strcmp(strrep(char(mdlIn.Formula),' ',''), strrep(form,' ',''))
            mdlOut = fitlme_df(Tall, form, dfmeth);
        end
    catch
        mdlOut = fitlme_df(Tall, form, dfmeth);
    end
end

function s = safe_dfmethod(mdl)
    try
        s = mdl.DFMethod;
        if isempty(s), s = 'N/A'; end
    catch
        s = 'N/A';
    end
end

function local_print_lme(mdl, label)
%LOCAL_PRINT_LME Safe, dependency-free LMM pretty-printer (no output args)
if nargin < 2, label = ''; end
try
    fprintf('\n==================== %s ====================\n', label);

    % Basics
    try, fprintf('Formula    : %s\n', char(mdl.Formula)); end
    try, fprintf('FitMethod  : %s\n', mdl.FitMethod);      end
    try, fprintf('LogLik     : %.4f\n', mdl.LogLikelihood); end

    % Fixed effects
    fprintf('\nFixed effects:\n');
    C = mdl.Coefficients;
    if istable(C)
        names = C.Properties.VariableNames;
    elseif isa(C,'dataset')
        names = get(C,'VarNames');
    else
        names = {};
    end
    keep  = intersect({'Name','Estimate','SE','DF','tStat','pValue'}, names, 'stable');
    if isempty(keep), disp(C); else, disp(C(:, keep)); end

    % Variance components (robust across releases)
    fprintf('\nCovariance parameters (variance components):\n');
    printed = false;
    try
        CP = covarianceParameters(mdl);
        if istable(CP), disp(CP); printed = true; end
    catch, end
    if ~printed
        try
            [psi, grp] = covarianceParameters(mdl); %#ok<ASGLU>
            fprintf('  psi (matrix/cell):\n'); disp(psi);
            printed = true;
        catch, end
    end
    if ~printed, fprintf('  (unavailable in this MATLAB release)\n'); end

    fprintf('=========================================================\n');
catch ME
    warning('local_print_lme failed: %s', ME.message);
end
end


function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end
