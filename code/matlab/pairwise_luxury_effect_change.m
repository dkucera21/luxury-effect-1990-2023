%% ===================== Script: Pairwise change (RAW income per $10k) =====================
% Compares the income–outcome slope between all year pairs, city-by-city,
% using RAW_INCOME_CPI scaled by $10,000 (per $10k).
% Outputs a struct Pairwise.(metric) with dBeta_*, p_*, used_* columns per pair.

% ---------------------- USER SETTINGS ----------------------
yearsUse          = [1990 2000 2010 2020 2023];
CoverageThreshold = 50;      % % overlap cutoff per tract
MinTractsPerCY    = 5;       % minimum valid tracts per city–year
pctNames          = {'PCT_OVERLAP'};  % coverage candidates (prefer PCT_OVERLAP)
useRobust         = false;   % robust regression option
useWinsor         = false;   % optional winsorization of X/Y within each year (on valid tracts only)
winzPct           = [1 99];  % use [5 95] automatically if n<50
savePairwiseCSV   = '';      % e.g., 'pairwise_summary.csv' ('' = no save)
% -----------------------------------------------------------

metrics = {'MEDIAN_NDVI','MEAN_NDVI','MEDIAN_LST','MEAN_LST'};
labels  = {'Median NDVI','Mean NDVI','Median LST','Mean LST'};

% ---------- inputs / guards ----------
assert(exist('CENSUS_TABLES_rebuilt','var')==1 && isstruct(CENSUS_TABLES_rebuilt), ...
    'CENSUS_TABLES_rebuilt not found or not a struct.');

% City list (mirror repo behavior; uncomment next line to force OLD CityNames)
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
% CityListMaster = string(CityNames(:));  % <-- uncomment to hard-match OLD script selection
CityListMaster = unique(strtrim(CityListMaster));

tblNames = fieldnames(CENSUS_TABLES_rebuilt);
numify   = @(x) str2double(regexprep(string(x),'[,\$%]',''));
winsorF  = @(v,lo,hi) min(max(v, lo), hi);

% ---------- pairwise interaction per city ----------
PAIRS    = nchoosek(yearsUse, 2);
pairTags = arrayfun(@(i) sprintf('%d_%d', PAIRS(i,1), PAIRS(i,2)), 1:size(PAIRS,1), 'uni', 0);

Pairwise = struct();
for m = 1:numel(metrics)
    yName = metrics{m};
    nC = numel(CityListMaster);
    Tm = table(CityListMaster(:), 'VariableNames', {'City'});
    for pj = 1:size(PAIRS,1)
        tag = pairTags{pj};
        Tm.(['dBeta_' tag]) = NaN(nC,1);  % (Δ slope) per $10k
        Tm.(['p_'     tag]) = NaN(nC,1);
        Tm.(['used_'  tag]) = false(nC,1);
    end

    for ci = 1:nC
        city = CityListMaster(ci);

        for pj = 1:size(PAIRS,1)
            a = PAIRS(pj,1); b = PAIRS(pj,2);
            tag = sprintf('%d_%d', a, b);

            % ---- gather year a ----
            tA = sprintf('T_%d_%s', a, city);
            if ~ismember(tA, tblNames), continue, end
            TA = CENSUS_TABLES_rebuilt.(tA); VA = string(TA.Properties.VariableNames);
            if ~ismember('RAW_INCOME_CPI', VA) || ~ismember(yName, VA), continue, end

            % coverage column (A): prefer PCT_OVERLAP, else first hit in pctNames
            vPctA = '';
            if any(strcmpi(VA,'PCT_OVERLAP'))
                vPctA = 'PCT_OVERLAP';
            else
                for k = 1:numel(pctNames)
                    j = find(strcmpi(VA, pctNames{k}), 1);
                    if ~isempty(j), vPctA = char(VA(j)); break, end
                end
            end
            if isempty(vPctA), continue, end

            xA_raw = numify(TA.('RAW_INCOME_CPI'));   % dollars
            xA     = xA_raw / 10000;              % per $10k
            yA     = numify(TA.(yName));
            pctA   = numify(TA.(vPctA));
            % If coverage looks like 0–1, convert to percent
            mxA = max(pctA(isfinite(pctA)));
            if isfinite(mxA) && mxA <= 1.5, pctA = 100*pctA; end

            % ---- gather year b ----
            tB = sprintf('T_%d_%s', b, city);
            if ~ismember(tB, tblNames), continue, end
            TB = CENSUS_TABLES_rebuilt.(tB); VB = string(TB.Properties.VariableNames);
            if ~ismember('RAW_INCOME_CPI', VB) || ~ismember(yName, VB), continue, end

            % coverage column (B): prefer PCT_OVERLAP, else first hit in pctNames
            vPctB = '';
            if any(strcmpi(VB,'PCT_OVERLAP'))
                vPctB = 'PCT_OVERLAP';
            else
                for k = 1:numel(pctNames)
                    j = find(strcmpi(VB, pctNames{k}), 1);
                    if ~isempty(j), vPctB = char(VB(j)); break, end
                end
            end
            if isempty(vPctB), continue, end

            xB_raw = numify(TB.('RAW_INCOME_CPI'));   % dollars
            xB     = xB_raw / 10000;              % per $10k
            yB     = numify(TB.(yName));
            pctB   = numify(TB.(vPctB));
            mxB = max(pctB(isfinite(pctB)));
            if isfinite(mxB) && mxB <= 1.5, pctB = 100*pctB; end

            % ---- validity masks: coverage + finite X/Y ----
            vA = isfinite(xA) & isfinite(yA) & isfinite(pctA) & pctA >= CoverageThreshold;
            vB = isfinite(xB) & isfinite(yB) & isfinite(pctB) & pctB >= CoverageThreshold;
            if nnz(vA) < MinTractsPerCY || nnz(vB) < MinTractsPerCY, continue, end
            if std(xA(vA),'omitnan')==0 || std(xB(vB),'omitnan')==0, continue, end

            % Optional winsorization (within-year, on valid tracts only)
            if useWinsor
                limA = winzPct; if nnz(vA)<50, limA = [5 95]; end
                limB = winzPct; if nnz(vB)<50, limB = [5 95]; end
                xA(vA) = winsorF(xA(vA), prctile(xA(vA),limA(1)), prctile(xA(vA),limA(2)));
                yA(vA) = winsorF(yA(vA), prctile(yA(vA),limA(1)), prctile(yA(vA),limA(2)));
                xB(vB) = winsorF(xB(vB), prctile(xB(vB),limB(1)), prctile(xB(vB),limB(2)));
                yB(vB) = winsorF(yB(vB), prctile(yB(vB),limB(1)), prctile(yB(vB),limB(2)));
                % recheck variance after winsor
                if std(xA(vA),'omitnan')==0 || std(xB(vB),'omitnan')==0, continue, end
            end

            % ---- stacked regression: Y ~ (RAW_INCOME_CPI per $10k) * Year ----
            X    = [xA(vA); xB(vB)];
            Y    = [yA(vA); yB(vB)];
            Year = categorical([repmat(a,nnz(vA),1); repmat(b,nnz(vB),1)], [a b], {num2str(a),num2str(b)});

            % OLD-parity guard: require ≥ 2*MinTractsPerCY total stacked rows
            if numel(Y) < 2*MinTractsPerCY, continue, end

            mdl = stableFit(table(X,Year,Y),'Y ~ X*Year', useRobust);

            C     = mdl.Coefficients;
            rows  = string(C.Properties.RowNames);
            k     = find(rows==("X:Year_"+string(b)) | rows==("Year_"+string(b)+":X"), 1, 'first');
            if isempty(k), continue, end

            Tm.(['dBeta_' tag])(ci) = C.Estimate(k);  % per $10k
            Tm.(['p_'     tag])(ci) = C.pValue(k);
            Tm.(['used_'  tag])(ci) = true;
        end
    end

    Pairwise.(yName) = Tm;
end

disp('— Finished pairwise change tables (RAW income per $10k, coverage-filtered) —');

% ---------- reporting ----------
AR = char(8594); AU = char(8593); AD = char(8595);
VN0   = Pairwise.(metrics{1}).Properties.VariableNames;
pairs = {};
for i = 1:numel(VN0)
    v = VN0{i};
    if strncmpi(v,'dBeta_',6)
        tag = v(7:end);
        if any(strcmpi(['p_'    tag], VN0)) && any(strcmpi(['used_' tag], VN0))
            pairs{end+1} = tag; %#ok<AGROW>
        end
    end
end
if isempty(pairs), fprintf('No pairwise columns detected.\n'); return, end
AB  = cellfun(@(s) sscanf(s,'%d_%d').', pairs, 'UniformOutput', false);
ABm = vertcat(AB{:}); [~, ord] = sortrows(ABm,[1 2]); pairs = pairs(ord);

% optional CSV summary rows
outRows = {};

for p = 1:numel(pairs)
    tag  = pairs{p}; yrs = sscanf(tag,'%d_%d').';
    sVar = ['dBeta_' tag]; pVar = ['p_' tag]; uVar = ['used_' tag];
    fprintf('\n=== %d%s%d (RAW income per $10k; cov>=%.0f%%, n>=%d) ===\n', yrs(1), AR, yrs(2), CoverageThreshold, MinTractsPerCY);
    for i = 1:numel(metrics)
        T = Pairwise.(metrics{i});
        if ~all(ismember({sVar,pVar,uVar}, T.Properties.VariableNames))
            fprintf('%-12s: columns missing for %s\n', labels{i}, tag); continue
        end
        used = logical(T.(uVar)); sl = str2double(string(T.(sVar))); pv = str2double(string(T.(pVar)));
        n = nnz(used);
        if n==0
            fprintf('%-12s:  0.0%% significant (0.0%% %s, 0.0%% %s) [n=0]\n', labels{i}, AU, AD);
            if ~isempty(savePairwiseCSV)
                outRows(end+1,:) = {metrics{i}, tag, 0, 0, 0, 0}; %#ok<AGROW>
            end
            continue
        end
        sig   = used & (pv < 0.05) & isfinite(sl);
        nSig  = nnz(sig); nUp = nnz(sig & sl > 0); nDown = nnz(sig & sl < 0);
        pctSig=100*nSig/n; pctUp=100*nUp/n; pctDown=100*nDown/n;
        fprintf('%-12s:  %4.1f%% significant (%4.1f%% %s, %4.1f%% %s)  [n=%d]\n', ...
            labels{i}, pctSig, pctUp, AU, pctDown, AD, n);

        if ~isempty(savePairwiseCSV)
            outRows(end+1,:) = {metrics{i}, tag, n, pctSig, pctUp, pctDown}; %#ok<AGROW>
        end
    end
end

if ~isempty(savePairwiseCSV)
    Summary = cell2table(outRows, ...
        'VariableNames', {'Metric','YearPair','N_used','Pct_Significant','Pct_Up','Pct_Down'});
    writetable(Summary, savePairwiseCSV);
    fprintf('Saved pairwise summary to %s\n', savePairwiseCSV);
end

% ===================== Helper =====================
function mdl = stableFit(tbl, formula, useRobust)
    % Fit linear model with optional robust regression and a safe fallback.
    try
        if useRobust
            mdl = fitlm(tbl, formula, 'RobustOpts','on');
        else
            mdl = fitlm(tbl, formula);
        end
    catch
        mdl = fitlm(tbl, formula);
    end
end
