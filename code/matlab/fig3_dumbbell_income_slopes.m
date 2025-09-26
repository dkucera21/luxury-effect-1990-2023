%% FIGURE 3 — Dumbbell (Barbell) Plot
% Income → MEAN NDVI & MEAN LST slopes (per $10k) for two years, per city.
% Uses CENSUS_TABLES_rebuilt and your master city list. No functions.

% ===================== USER SETTINGS =====================
yearA = 1990;
yearB = 2023;

thresh    = 50;        % % overlap cutoff
minTracts = 5;         % minimum valid tracts per city-year
useWinsor = false;     % trim tails to reduce leverage
winz      = [1 99];    % winsor limits (percentiles)

excludeCities = lower(string({'Honolulu','Anchorage'}));  % optional excludes
savePNG = '';                 % e.g., 'fig3_dumbbell.png' ('' = do not save)
% ==========================================================

% ---------- checks / inputs ----------
assert(exist('CENSUS_TABLES_rebuilt','var')==1 && isstruct(CENSUS_TABLES_rebuilt), ...
    'CENSUS_TABLES_rebuilt not found or not a struct in the workspace.');

% Master city list (prefer CityNames_master)
if exist('CityNames_master','var') && ~isempty(CityNames_master)
    CityListMaster = string(CityNames_master(:));
elseif exist('MasterCities','var') && istable(MasterCities) && ismember('City', MasterCities.Properties.VariableNames)
    CityListMaster = string(MasterCities.City);
elseif exist('CityNames','var') && ~isempty(CityNames)
    warning('Master city list not found; falling back to CityNames (all cities).');
    CityListMaster = string(CityNames(:));
else
    error('No valid city list found (CityNames_master, MasterCities.City, or CityNames).');
end
CityListMaster = unique(strtrim(CityListMaster));

% ---------- config ----------
ovlpCandidates = {'PCT_OVERLAP','PCT_COVER','PCT_COV','PercentAreaNDVI'};

% Correct RAW income headers by year (first good match wins)
incByYear = containers.Map( ...
  {'1990','2000','2010','2020','2023'}, ...
   {{'P080A001'}, ...
    {'DP3_C112'}, ...
    {'S1903_C02_001E'}, ...
    {'S1903_C03_001E'}, ...
    {'S1903_C03_001E'}} );

numify  = @(x) str2double(regexprep(string(x),'[,\$%]',''));   % numeric cleaner
stripYR = @(s) regexprep(string(s),'_(?:19|20)\d{2}$','');      % drop trailing _YYYY
tblNames = fieldnames(CENSUS_TABLES_rebuilt);

% ---------- compute slopes ----------
CityKeep = string.empty(0,1);
NDVI_A = []; NDVI_B = [];
LST_A  = []; LST_B  = [];

for ci = 1:numel(CityListMaster)
    city = string(CityListMaster(ci));
    sND = [NaN NaN];  % [yearA yearB]
    sLS = [NaN NaN];

    for yy = 1:2
        yr = [yearA yearB]; yr = yr(yy);
        tname = sprintf('T_%d_%s', yr, city);
        if ~ismember(tname, tblNames), continue; end

        T = CENSUS_TABLES_rebuilt.(tname);
        if ~istable(T), continue; end
        V = T.Properties.VariableNames;

        % ----- overlap column (first match wins) -----
        vPct = '';
        for k = 1:numel(ovlpCandidates)
            j = find(strcmpi(V, ovlpCandidates{k}), 1, 'first');
            if ~isempty(j), vPct = V{j}; break; end
        end
        if isempty(vPct), continue; end
        pct = numify(T.(vPct));

        % ----- income (raw dollars for that year) → per $10k -----
        vInc = '';
        incCands = incByYear(num2str(yr));
        for k = 1:numel(incCands)
            j = find(strcmpi(V, incCands{k}), 1, 'first');
            if ~isempty(j), vInc = V{j}; break; end
        end
        if isempty(vInc), continue; end
        incRaw = numify(T.(vInc));   % raw $ (not CPI-adjusted)

        % ----- responses -----
        jNd = find(strcmpi(V,'MEAN_NDVI'), 1, 'first');
        jLs = find(strcmpi(V,'MEAN_LST'),  1, 'first');
        if isempty(jNd) || isempty(jLs), continue; end
        yND = numify(T.(V{jNd}));
        yLS = numify(T.(V{jLs}));

        % ----- masks -----
        validND = pct >= thresh & isfinite(incRaw) & isfinite(yND);
        validLS = pct >= thresh & isfinite(incRaw) & isfinite(yLS);
        if nnz(validND) < minTracts || nnz(validLS) < minTracts, continue; end

        % ----- prepare X (per $10k) and Y -----
        xND = incRaw(validND)/10000; yNDv = yND(validND);
        xLS = incRaw(validLS)/10000; yLSv = yLS(validLS);

        if useWinsor
            p = prctile(xND, winz);  xND  = min(max(xND, p(1)), p(2));
            p = prctile(yNDv, winz); yNDv = min(max(yNDv, p(1)), p(2));
            p = prctile(xLS, winz);  xLS  = min(max(xLS, p(1)), p(2));
            p = prctile(yLSv, winz); yLSv = min(max(yLSv, p(1)), p(2));
        end

        % need some variation in X
        if std(xND,'omitnan')==0 || std(xLS,'omitnan')==0, continue; end

        % ----- OLS slopes -----
        mdlND = fitlm(xND, yNDv);  sND(yy) = mdlND.Coefficients.Estimate(2);
        mdlLS = fitlm(xLS, yLSv);  sLS(yy) = mdlLS.Coefficients.Estimate(2);
    end

    % keep only fully observed cities (both metrics in both years)
    if all(isfinite(sND)) && all(isfinite(sLS))
        CityKeep(end+1,1) = city;                 %#ok<AGROW>
        NDVI_A(end+1,1)   = sND(1); NDVI_B(end+1,1) = sND(2); %#ok<AGROW>
        LST_A (end+1,1)   = sLS(1); LST_B (end+1,1) = sLS(2); %#ok<AGROW>
    end
end

% ---------- build plot table ----------
Tplot = table(CityKeep, NDVI_A, NDVI_B, LST_A, LST_B, ...
              'VariableNames', {'City','NDVI_A','NDVI_B','LST_A','LST_B'});

% Exclude cities (robust to case/suffix)
maskExcl = false(height(Tplot),1);
for ex = excludeCities(:)'
    maskExcl = maskExcl | contains(lower(Tplot.City), ex);
end
Tplot(maskExcl,:) = [];

% Sort by earlier NDVI slope (descending) and set y positions
[~, ord] = sort(Tplot.NDVI_A, 'descend');
Tplot = Tplot(ord,:);
nCities = height(Tplot);
y = (nCities:-1:1)';

% ---------- Clean and map labels ----------
labelsLeft = stripYR(flip(Tplot.City));       % remove trailing _YYYY
labelsLeft = strrep(labelsLeft, '_', ' ');    % simple prettification

% ------------ Draw figure ------------
fig = figure('Color','w','Units','pixels','Position',[100 100 1200 800]);

% (a) Income → MEAN NDVI
subplot(1,2,1); hold on;
for i = 1:nCities
    plot([Tplot.NDVI_A(i) Tplot.NDVI_B(i)], [y(i) y(i)], '-', 'Color',[0.7 0.7 0.7], 'LineWidth',1);
    plot(Tplot.NDVI_A(i), y(i), 'o', 'MarkerFaceColor',[0 0.45 0.74], 'MarkerEdgeColor','k', 'MarkerSize',6);  % yearA (blue)
    plot(Tplot.NDVI_B(i), y(i), 'o', 'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','k', 'MarkerSize',6); % yearB (orange)
end
yticks(1:nCities); yticklabels(labelsLeft);
xlabel('Income–NDVI (NDVI per $10,000)'); 
title(sprintf('(a) Income–NDVI slope (%d vs %d)', yearA, yearB));
xline(0,'k-');
xlim([min([Tplot.NDVI_A;Tplot.NDVI_B],[],'omitnan')-0.01, max([Tplot.NDVI_A;Tplot.NDVI_B],[],'omitnan')+0.01]);
ylim([0.5 nCities+0.5]); set(gca,'FontSize',8,'Box','off');

% ---- explicit legend handles so 1990≡blue, 2023≡orange, lines excluded
hA = plot(NaN,NaN,'o','MarkerFaceColor',[0 0.45 0.74],'MarkerEdgeColor','k','MarkerSize',6);
hB = plot(NaN,NaN,'o','MarkerFaceColor',[0.85 0.33 0.10],'MarkerEdgeColor','k','MarkerSize',6);
legend([hA hB], {num2str(yearA), num2str(yearB)}, 'Location','best');
legend boxoff; 
hold off;

% (b) Income → MEAN LST (same ordering; hide y labels)
subplot(1,2,2); hold on;
for i = 1:nCities
    plot([Tplot.LST_A(i) Tplot.LST_B(i)], [y(i) y(i)], '-', 'Color',[0.7 0.7 0.7], 'LineWidth',1);
    plot(Tplot.LST_A(i), y(i), 'o', 'MarkerFaceColor',[0 0.45 0.74], 'MarkerEdgeColor','k', 'MarkerSize',6);   % yearA (blue)
    plot(Tplot.LST_B(i), y(i), 'o', 'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','k', 'MarkerSize',6); % yearB (orange)
end
yticks(1:nCities); yticklabels([]);
xlabel('Income–LST (\circC per $10,000)'); 
title(sprintf('(b) Income–LST slope (%d vs %d)', yearA, yearB));
xline(0,'k-');
xlim([min([Tplot.LST_A;Tplot.LST_B],[],'omitnan')-0.1, max([Tplot.LST_A;Tplot.LST_B],[],'omitnan')+0.1]);
ylim([0.5 nCities+0.5]); set(gca,'FontSize',8,'Box','off'); 
hold off;

fprintf('Dumbbell plot includes %d cities with valid %d and %d slopes for both metrics (master-set + exclusions applied).\n', ...
        nCities, yearA, yearB);

% optional save
if ~isempty(savePNG)
    exportgraphics(fig, savePNG, 'Resolution',300);
    fprintf('Saved figure to %s\n', savePNG);
end
