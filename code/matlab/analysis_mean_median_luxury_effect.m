%% ===== Script: Mean/Median luxury effect across all cities & years (per $10k) =====
% For each city-year, fit OLS: y ~ (RAW_INCOME_CPI/10000).
% Collect slopes (per $10k) + p-values, then:
%  - report pooled and city-averaged summaries
%  - report % significant (p<0.05) by year
%  - report 1990→2023 change and "percent fewer" significant cities

% ---------------------- USER SETTINGS ----------------------
yearsUse  = [1990 2000 2010 2020 2023];
threshold = 50;          % % overlap cutoff
minTracts = 5;           % min valid tracts per city-year
pctNames  = {'PCT_OVERLAP','PCT_COVER','PCT_COV','PercentAreaNDVI'};
winsorize = false;       % light trimming
winLim    = [1 99];
alphaSig  = 0.05;        % significance threshold
saveSummaryCSV = '';     % e.g., 'luxury_summary_per10k.csv' ('' = no save)
% -----------------------------------------------------------

metrics = struct( ...
  'y',   {'MEDIAN_NDVI','MEAN_NDVI','MEDIAN_LST','MEAN_LST'}, ...
  'name',{'Median NDVI','Mean NDVI','Median LST','Mean LST'});

% ---------- inputs / guards ----------
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

tblNames = string(fieldnames(CENSUS_TABLES_rebuilt));
numcol   = @(x) str2double(regexprep(string(x),'[,\$%]',''));

% ---------- per–city-year slopes (per $10k) ----------
LuxPerCityYear_per10k = table( ...
    strings(0,1), strings(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), false(0,1), ...
    'VariableNames', {'Metric','City','Year','N','Slope_per10k','pValue','Significant'});

for ci = 1:numel(CityListMaster)
    city = CityListMaster(ci);

    for yr = yearsUse
        tname = sprintf('T_%d_%s', yr, city);
        if ~ismember(string(tname), tblNames), continue, end

        T = CENSUS_TABLES_rebuilt.(tname);
        if ~istable(T) || height(T)==0, continue, end
        V = string(T.Properties.VariableNames);

        % overlap
        jPct = find(ismember(lower(V), lower(string(pctNames))), 1, 'first');
        if isempty(jPct), continue, end
        pct = numcol(T.(char(V(jPct))));

        % RAW income → per $10k
        if ~ismember("RAW_INCOME_CPI", V), continue, end
        incRaw = numcol(T.("RAW_INCOME_CPI"));

        for m = 1:numel(metrics)
            jY = find(strcmpi(V, metrics(m).y), 1, 'first');
            if isempty(jY), continue, end
            yv = numcol(T.(char(V(jY))));

            valid = pct >= threshold & isfinite(incRaw) & isfinite(yv);
            if nnz(valid) < minTracts, continue, end

            x = incRaw(valid) / 10000;   % per $10k
            y = yv(valid);

            if winsorize
                loX = prctile(x, winLim(1)); hiX = prctile(x, winLim(2));
                loY = prctile(y, winLim(1)); hiY = prctile(y, winLim(2));
                x = min(max(x, loX), hiX);
                y = min(max(y, loY), hiY);
            end
            if std(x,'omitnan') == 0, continue, end

            mdl = fitlm(x, y);                          % OLS
            b10 = mdl.Coefficients.Estimate(2);         % slope per $10k
            pB  = mdl.Coefficients.pValue(2);           % p-value for slope

            LuxPerCityYear_per10k = [LuxPerCityYear_per10k; ...
                {string(metrics(m).name), string(city), yr, numel(x), b10, pB, isfinite(pB) && pB<alphaSig}]; %#ok<AGROW>
        end
    end
end

% ---------- summaries (pooled & city-averaged) ----------
summaryRows = strings(0,1);
nPooled = []; meanP = []; medP = [];
nCity   = []; meanC = [];  medC = [];
pMeanC  = []; pMedC = [];

for m = 1:numel(metrics)
    rowsM = LuxPerCityYear_per10k(strcmp(LuxPerCityYear_per10k.Metric, metrics(m).name), :);
    v = rowsM.Slope_per10k;
    c = rowsM.City;

    valid = isfinite(v);
    summaryRows(end+1,1) = string(metrics(m).name); %#ok<AGROW>
    nPooled(end+1,1)     = nnz(valid);              %#ok<AGROW>
    meanP(end+1,1)       = mean(v(valid),'omitnan');   %#ok<AGROW>
    medP(end+1,1)        = median(v(valid),'omitnan'); %#ok<AGROW>

    % --- city-averaged values & p-values ---
    if any(valid)
        cityU = unique(c(valid));
        vCity = nan(numel(cityU),1);
        for i = 1:numel(cityU)
            mask = valid & c==cityU(i);
            vCity(i) = mean(v(mask),'omitnan');
        end

        % keep only finite city means
        okC = isfinite(vCity);
        vCity = vCity(okC);

        nCity(end+1,1) = numel(vCity);                %#ok<AGROW>
        meanC(end+1,1) = mean(vCity,'omitnan');       %#ok<AGROW>
        medC(end+1,1)  = median(vCity,'omitnan');     %#ok<AGROW>

        % mean(cityAvg) ≠ 0 (two-sided t-test)
        p_t = NaN;
        if numel(vCity) >= 2
            try
                [~,p_t] = ttest(vCity, 0, 'Alpha', 0.05);
            catch
                % simple fallback: z-test using sample mean/SE
                mu = mean(vCity); se = std(vCity)/sqrt(numel(vCity));
                if isfinite(se) && se>0
                    z = mu/se; p_t = 2*(1 - normcdf(abs(z),0,1));
                end
            end
        end

        % median(cityAvg) ≠ 0 (Wilcoxon signed-rank; fallback to sign test)
        p_sr = NaN;
        try
            if exist('signrank','file')==2 && numel(vCity) >= 5
                p_sr = signrank(vCity, 0, 'method','approx');  % two-sided
            else
                % Binomial sign test as a fallback
                kpos = sum(vCity > 0);
                kneg = sum(vCity < 0);
                nNZ  = kpos + kneg;
                if nNZ > 0
                    % two-sided binomial about 0.5
                    p_one = binocdf(min(kpos,kneg), nNZ, 0.5);
                    p_sr  = 2 * p_one;
                    p_sr  = min(1, p_sr);
                end
            end
        catch
            % leave NaN if tests fail
        end

        pMeanC(end+1,1) = p_t;  %#ok<AGROW>
        pMedC(end+1,1)  = p_sr; %#ok<AGROW>
    else
        nCity(end+1,1) = 0; meanC(end+1,1)=NaN; medC(end+1,1)=NaN;
        pMeanC(end+1,1)=NaN; pMedC(end+1,1)=NaN; %#ok<AGROW>
    end
end

LuxSummary_per10k = table( ...
    summaryRows, nPooled, meanP, medP, nCity, meanC, medC, pMeanC, pMedC, ...
    'VariableNames', {'Metric','N_cityYears','Mean_pooled','Median_pooled', ...
                      'N_cities','Mean_cityAvg','Median_cityAvg', ...
                      'pMean_cityAvg','pMedian_cityAvg'});

disp('— Mean/Median luxury effect across cities & years (units: per $10,000) —');
disp(LuxSummary_per10k);

if ~isempty(saveSummaryCSV)
    writetable(LuxSummary_per10k, saveSummaryCSV);
    fprintf('Saved summary to %s\n', saveSummaryCSV);
end


% ---------- % significant by year (per metric) ----------
yearsAll = yearsUse(:);
metricsNames = string({metrics.name})';
PercentSig_byYear = table(yearsAll, 'VariableNames', {'Year'});

for m = 1:numel(metrics)
    mname = metrics(m).name;
    pctSigCol = nan(numel(yearsAll),1);
    for yi = 1:numel(yearsAll)
        yr = yearsAll(yi);
        rows = LuxPerCityYear_per10k(strcmp(LuxPerCityYear_per10k.Metric, mname) & LuxPerCityYear_per10k.Year==yr, :);
        nTot = height(rows);
        if nTot==0
            pctSigCol(yi) = NaN;
        else
            nSig = nnz(rows.Significant);
            pctSigCol(yi) = 100 * (nSig / nTot);
        end
    end
    PercentSig_byYear.(matlab.lang.makeValidName("pctSig_"+mname)) = pctSigCol;
end

disp('— Percent significant luxury-effect relationships by year (p < 0.05) —');
disp(PercentSig_byYear);

% ---------- 1990→2023 delta & "percent fewer" ----------
getPct = @(col,yr) PercentSig_byYear.(col)(PercentSig_byYear.Year==yr);
Delta_1990_2023 = table(metricsNames, zeros(numel(metricsNames),1), zeros(numel(metricsNames),1), zeros(numel(metricsNames),1), ...
    'VariableNames', {'Metric','PctSig_1990','PctSig_2023','Delta_2023_minus_1990'});

for m = 1:numel(metrics)
    col = matlab.lang.makeValidName("pctSig_"+metrics(m).name);
    p1990 = getPct(col, 1990); if isempty(p1990), p1990 = NaN; end
    p2023 = getPct(col, 2023); if isempty(p2023), p2023 = NaN; end
    Delta_1990_2023.PctSig_1990(m) = p1990;
    Delta_1990_2023.PctSig_2023(m) = p2023;
    Delta_1990_2023.Delta_2023_minus_1990(m) = p2023 - p1990;
end
Delta_1990_2023.PctFewer_2023_vs_1990 = Delta_1990_2023.PctSig_1990 - Delta_1990_2023.PctSig_2023;

disp('— 1990→2023 change in % significant and "% fewer" in 2023 vs 1990 —');
disp(Delta_1990_2023);
