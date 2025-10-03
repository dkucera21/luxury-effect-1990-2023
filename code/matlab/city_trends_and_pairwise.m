%% city_trends_and_pairwise.m
% Combines:
%   (A) City-level luxury-effect trends (slopes, FDR, weakening shares)
%   (B) Pairwise year-to-year changes in luxury-effect slopes (per $10k)
%
% Requires these to be in the workspace (as produced by your pipeline):
%   - LUX_MEAN_T and/or LUX_MEDN_T with vars:
%       City, Year, LUX_NDVI_per10k, LUX_LST_per10k
%   - CENSUS_TABLES_rebuilt (struct of city-year tract tables) for (B)
%
% IMPORTANT: No clearvars; nothing is deleted from your workspace.

%% ===================== (A) City-level trends & weakening =====================
alpha   = 0.05;   % FDR target per metric
minYears = 3;     % require >=3 unique years per city
baseYear = 1990;  % year centering for slopes
useFDR  = true;

haveMEAN = exist('LUX_MEAN_T','var')==1 && istable(LUX_MEAN_T);
haveMEDN = exist('LUX_MEDN_T','var')==1 && istable(LUX_MEDN_T);
if ~haveMEAN && ~haveMEDN
    warning('Section (A) skipped: need LUX_MEAN_T and/or LUX_MEDN_T in workspace.');
else
    if haveMEAN
        Trends_MEAN = local_analyze_one(LUX_MEAN_T,'MEAN',alpha,minYears,baseYear,useFDR);
        assignin('base','Trends_MEAN',Trends_MEAN);
    end
    if haveMEDN
        Trends_MEDN = local_analyze_one(LUX_MEDN_T,'MEDIAN',alpha,minYears,baseYear,useFDR);
        assignin('base','Trends_MEDN',Trends_MEDN);
    end
    fprintf('\n(A) City-level trends placed in workspace: %s%s\n', ...
        tern(haveMEAN,'Trends_MEAN ',''), tern(haveMEDN,'Trends_MEDN',''));
end

%% =========== (B) Pairwise change (RAW_INCOME_log per $10k) ===================
% Compares the income–outcome slope between all year pairs, city-by-city,
% using RAW_INCOME_log / 10000 (per $10k) and a stacked interaction model.

yearsUse          = [1990 2000 2010 2020 2023];
CoverageThreshold = 50;      % % overlap cutoff per tract
MinTractsPerCY    = 5;       % minimum valid tracts per city–year
pctNames          = {'PCT_OVERLAP'};  % coverage fallbacks
useRobust         = false;   % robust fit option
useWinsor         = false;   % optional winsorization inside each year
winzPct           = [1 99];  % [5 95] auto if n<50
savePairwiseCSV   = '';      % e.g., 'pairwise_summary.csv' ('' = no save)

metrics = {'MEDIAN_NDVI','MEAN_NDVI','MEDIAN_LST','MEAN_LST'};
labels  = {'Median NDVI','Mean NDVI','Median LST','Mean LST'};

if ~(exist('CENSUS_TABLES_rebuilt','var')==1 && isstruct(CENSUS_TABLES_rebuilt))
    warning('Section (B) skipped: CENSUS_TABLES_rebuilt not found.');
else
    % City list (mirror your repo fallbacks)
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
    numify   = @(x) str2double(regexprep(string(x),'[,\$%]',''));
    winsorF  = @(v,lo,hi) min(max(v, lo), hi);

    % Pairs
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

                % ---- Year a ----
                tA = sprintf('T_%d_%s', a, city);
                if ~ismember(tA, tblNames), continue, end
                TA = CENSUS_TABLES_rebuilt.(tA); VA = string(TA.Properties.VariableNames);
                if ~ismember('RAW_INCOME_log', VA) || ~ismember(yName, VA), continue, end

                vPctA = '';
                if any(strcmpi(VA,'PCT_OVERLAP')), vPctA = 'PCT_OVERLAP';
                else
                    for k = 1:numel(pctNames)
                        j = find(strcmpi(VA, pctNames{k}), 1);
                        if ~isempty(j), vPctA = char(VA(j)); break, end
                    end
                end
                if isempty(vPctA), continue, end

                xA = numify(TA.('RAW_INCOME_log'))/10000;  % per $10k
                yA = numify(TA.(yName));
                pctA = numify(TA.(vPctA));
                mxA = max(pctA(isfinite(pctA)));
                if isfinite(mxA) && mxA <= 1.5, pctA = 100*pctA; end

                % ---- Year b ----
                tB = sprintf('T_%d_%s', b, city);
                if ~ismember(tB, tblNames), continue, end
                TB = CENSUS_TABLES_rebuilt.(tB); VB = string(TB.Properties.VariableNames);
                if ~ismember('RAW_INCOME_log', VB) || ~ismember(yName, VB), continue, end

                vPctB = '';
                if any(strcmpi(VB,'PCT_OVERLAP')), vPctB = 'PCT_OVERLAP';
                else
                    for k = 1:numel(pctNames)
                        j = find(strcmpi(VB, pctNames{k}), 1);
                        if ~isempty(j), vPctB = char(VB(j)); break, end
                    end
                end
                if isempty(vPctB), continue, end

                xB = numify(TB.('RAW_INCOME_log'))/10000;  % per $10k
                yB = numify(TB.(yName));
                pctB = numify(TB.(vPctB));
                mxB = max(pctB(isfinite(pctB)));
                if isfinite(mxB) && mxB <= 1.5, pctB = 100*pctB; end

                % ---- Validity ----
                vA = isfinite(xA) & isfinite(yA) & isfinite(pctA) & pctA >= CoverageThreshold;
                vB = isfinite(xB) & isfinite(yB) & isfinite(pctB) & pctB >= CoverageThreshold;
                if nnz(vA) < MinTractsPerCY || nnz(vB) < MinTractsPerCY, continue, end
                if std(xA(vA),'omitnan')==0 || std(xB(vB),'omitnan')==0, continue, end

                % Optional winsorization
                if useWinsor
                    limA = winzPct; if nnz(vA)<50, limA = [5 95]; end
                    limB = winzPct; if nnz(vB)<50, limB = [5 95]; end
                    xA(vA) = winsorF(xA(vA), prctile(xA(vA),limA(1)), prctile(xA(vA),limA(2)));
                    yA(vA) = winsorF(yA(vA), prctile(yA(vA),limA(1)), prctile(yA(vA),limA(2)));
                    xB(vB) = winsorF(xB(vB), prctile(xB(vB),limB(1)), prctile(xB(vB),limB(2)));
                    yB(vB) = winsorF(yB(vB), prctile(yB(vB),limB(1)), prctile(yB(vB),limB(2)));
                    if std(xA(vA),'omitnan')==0 || std(xB(vB),'omitnan')==0, continue, end
                end

                % ---- Stacked interaction: Y ~ X*Year (Year is {a,b}) ----
                X    = [xA(vA); xB(vB)];
                Y    = [yA(vA); yB(vB)];
                Year = categorical([repmat(a,nnz(vA),1); repmat(b,nnz(vB),1)], [a b], {num2str(a),num2str(b)});

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

    % ---------- Reporting ----------
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
    if isempty(pairs)
        fprintf('No pairwise columns detected.\n');
    else
        AB  = cellfun(@(s) sscanf(s,'%d_%d').', pairs, 'UniformOutput', false);
        ABm = vertcat(AB{:}); [~, ord] = sortrows(ABm,[1 2]); pairs = pairs(ord);

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
    end

    assignin('base','Pairwise',Pairwise);
    fprintf('\n(B) Pairwise struct placed in workspace: Pairwise\n');
end

%% ===================== Local helper functions (do not remove) =====================
function Trends = local_analyze_one(TT, label, alpha, minYears, baseYear, useFDR)
    req = ["City","Year","LUX_NDVI_per10k","LUX_LST_per10k"];
    v = string(TT.Properties.VariableNames);
    missing = setdiff(req, v);
    if ~isempty(missing)
        error('%s is missing columns: %s', label, strjoin(missing,', '));
    end

    TT = TT(:, req);
    TT.City = string(TT.City);
    TT.Year = double(TT.Year);

    metrics = struct('col',{'LUX_NDVI_per10k','LUX_LST_per10k'}, ...
                     'name',{'NDVI','LST'});

    rows = [];

    for m = 1:numel(metrics)
        ycol  = metrics(m).col;
        mname = metrics(m).name;

        cities = unique(TT.City,'stable');
        for i = 1:numel(cities)
            Ci  = cities(i);
            S   = TT(TT.City==Ci & isfinite(TT.(ycol)) & isfinite(TT.Year), :);
            if height(S) < minYears || numel(unique(S.Year)) < minYears
                continue
            end

            S.yearC = S.Year - baseYear;
            mdl = fitlm(S, sprintf('%s ~ yearC', ycol), 'Intercept', true);
            if ~ismember('yearC', string(mdl.Coefficients.Properties.RowNames))
                continue
            end

            slope   = mdl.Coefficients{'yearC','Estimate'};
            seSlope = mdl.Coefficients{'yearC','SE'};
            pSlope  = mdl.Coefficients{'yearC','pValue'};
            r2      = mdl.Rsquared.Ordinary;

            % weakening / strengthening direction
            if mname=="NDVI"
                weaken_dir = slope < 0;
                strengthen_dir = slope > 0;
            else
                weaken_dir = slope > 0;
                strengthen_dir = slope < 0;
            end

            % Percent weakening across observed span (move toward zero)
            yrMin = min(S.Year); yrMax = max(S.Year);
            yhatMin = predict(mdl, table(yrMin - baseYear, 'VariableNames',{'yearC'}));
            yhatMax = predict(mdl, table(yrMax - baseYear, 'VariableNames',{'yearC'}));
            denom = max(abs(yhatMin), 1e-8);
            pctWeak = 1 - (abs(yhatMax)/denom);

            rows = [rows; struct( ...
                'Table', string(label), ...
                'Metric', string(mname), ...
                'City', Ci, ...
                'nYears', numel(unique(S.Year)), ...
                'YearMin', yrMin, ...
                'YearMax', yrMax, ...
                'Slope_perYear', slope, ...
                'SE', seSlope, ...
                'pSlope', pSlope, ...
                'Weaken_dir', weaken_dir, ...
                'Strengthen_dir', strengthen_dir, ...
                'PctWeakeningSpan', double(pctWeak), ...
                'R2', r2 )]; %#ok<AGROW>
        end
    end

    if isempty(rows)
        Trends = TT(:,0); % empty table with no rows
        return
    end

    Trends = struct2table(rows);

    % FDR within each metric
    if useFDR && ~isempty(Trends)
        Trends.pSlope_FDR = NaN(height(Trends),1);
        for mname = ["NDVI","LST"]
            k = Trends.Metric==mname & isfinite(Trends.pSlope);
            if any(k)
                p = Trends.pSlope(k);
                [~,~,pFDR] = local_fdr_bh(p, alpha);
                Trends.pSlope_FDR(k) = pFDR;
            end
        end
    else
        Trends.pSlope_FDR = Trends.pSlope;
    end

    Trends.Weaken_sig     = Trends.Weaken_dir     & (Trends.pSlope_FDR < alpha);
    Trends.Strengthen_sig = Trends.Strengthen_dir & (Trends.pSlope_FDR < alpha);

    % Print summary
    for mname = ["NDVI","LST"]
        K = Trends.Metric==mname;
        if ~any(K), continue; end
        nC = nnz(K);
        pctWeak_dir = 100*mean(Trends.Weaken_dir(K),'omitnan');
        pctWeak_sig = 100*mean(Trends.Weaken_sig(K),'omitnan');
        nStrong_sig = nnz(Trends.Strengthen_sig(K));

        medW = median(Trends.PctWeakeningSpan(K),'omitnan')*100;
        q25  = quantile(Trends.PctWeakeningSpan(K),0.25);
        q75  = quantile(Trends.PctWeakeningSpan(K),0.75);

        fprintf('\n[%s – %s]\n', label, mname);
        fprintf(' Cities analyzed: %d (>= %d years)\n', nC, minYears);
        fprintf(' Weakened (direction only): %.1f%%%%\n', pctWeak_dir);
        fprintf(' Weakened & significant (BH-FDR @ alpha=%.2f): %.1f%%%%\n', alpha, pctWeak_sig);
        if nStrong_sig==0
            fprintf(' Significant strengthening (BH-FDR): none\n');
        else
            fprintf(' Significant strengthening (BH-FDR): %d cities\n', nStrong_sig);
        end
        fprintf(' Median %% weakening over observed span: %.1f%%%%  (IQR %.1f–%.1f%%%%)\n', ...
            medW, q25*100, q75*100);
    end
end

function [h, crit_p, adj_p] = local_fdr_bh(pvals, q)
% Benjamini–Hochberg FDR for a vector of p-values.
    if nargin<2, q = 0.05; end
    p = pvals(:);
    [ps, idx] = sort(p);
    m = numel(p);
    thresh = (1:m)' * (q/m);
    rej = ps <= thresh;
    if any(rej)
        kmax = find(rej,1,'last');
        crit_p = thresh(kmax);
        hvec = p <= crit_p;
    else
        crit_p = 0;
        hvec = false(size(p));
    end
    % adjusted p-values (monotone)
    adj = zeros(size(ps));
    adj(end) = ps(end);
    for i = m-1:-1:1
        adj(i) = min((m/i)*ps(i), adj(i+1));
    end
    adj_p = zeros(size(p)); adj_p(idx) = adj;
    h = hvec;
end

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

function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end
