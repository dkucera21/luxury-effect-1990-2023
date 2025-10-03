%% fig3_dumbbell_income_slopes.m
% ------------------------------------------------------------------------------
% FIGURE 3 — Dumbbell (Barbell) Plot
% Income (RAW_INCOME_CPI, scaled per $10k) → NDVI & LST slopes for two years, per city
% Tract filters: PCT_OVERLAP >= OverlapThresh and at least MinTracts valid tracts
% Choose MEAN_* or MEDIAN_* outcome columns via metricMode.
% ------------------------------------------------------------------------------

%% ===================== USER SETTINGS =====================
yearA        = 1990;
yearB        = 2023;

MinTracts    = 5;            % minimum valid tracts per city-year to fit a slope
OverlapThresh = 50;          % threshold for PCT_OVERLAP (use 0–100 scale; code handles 0–1)
metricMode   = 'MEDIAN';     % 'MEAN' or 'MEDIAN'

useWinsor    = false;        % trim tails to reduce leverage
winz         = [1 99];       % winsor limits (percentiles)

excludeCities = lower(string({'Honolulu','Anchorage'}));  % optional excludes
savePNG       = '';          % e.g., 'fig3_dumbbell.png' ('' = do not save)
%% =========================================================

% ---------- checks / inputs ----------
assert(exist('CENSUS_TABLES_rebuilt','var')==1 && isstruct(CENSUS_TABLES_rebuilt), ...
    'CENSUS_TABLES_rebuilt not found or not a struct in the workspace.');

% Master city list (prefer CityListMaster)
if exist('CityListMaster','var') && ~isempty(CityListMaster)
    CityListMaster = string(CityListMaster(:));
elseif exist('MasterCities','var') && istable(MasterCities) && ismember('City', MasterCities.Properties.VariableNames)
    CityListMaster = string(MasterCities.City);
elseif exist('CityNames','var') && ~isempty(CityNames)
    warning('Master city list not found; falling back to CityNames (all cities).');
    CityListMaster = string(CityNames(:));
else
    error('No valid city list found (CityListMaster, MasterCities.City, or CityNames).');
end
CityListMaster = unique(strtrim(CityListMaster));

tblNames = string(fieldnames(CENSUS_TABLES_rebuilt));
toNum    = @(x) str2double(regexprep(string(x),'[,\$%]',''));   % robust numeric cleaner

% Metric columns based on metricMode
metricMode = upper(string(metricMode));
switch metricMode
    case "MEAN"
        ndviVar = "MEAN_NDVI";
        lstVar  = "MEAN_LST";
        metricLabel = "MEAN";
    case "MEDIAN"
        ndviVar = "MEDIAN_NDVI";
        lstVar  = "MEDIAN_LST";
        metricLabel = "MEDIAN";
    otherwise
        error('metricMode must be ''MEAN'' or ''MEDIAN''.');
end

% ---------- compute slopes ----------
CityKeep = string.empty(0,1);
NDVI_A = []; NDVI_B = [];
LST_A  = []; LST_B  = [];

for ci = 1:numel(CityListMaster)
    city = string(CityListMaster(ci));
    sND = [NaN NaN];  % [yearA, yearB]
    sLS = [NaN NaN];

    for yy = 1:2
        yr = [yearA yearB];
        yr = yr(yy);

        % Accept T_YYYY_city or T_YYYYn_city
        t1 = sprintf('T_%d_%s',  yr, city);
        t2 = sprintf('T_%dn_%s', yr, city);
        if ismember(t1, tblNames)
            T = CENSUS_TABLES_rebuilt.(t1);
        elseif ismember(t2, tblNames)
            T = CENSUS_TABLES_rebuilt.(t2);
        else
            continue
        end
        if ~istable(T), continue; end
        V = string(T.Properties.VariableNames);

        % Require RAW_INCOME_CPI and outcome columns
        need = ["RAW_INCOME_CPI", ndviVar, lstVar];
        if ~all(ismember(need, V)), continue; end

        % Income per $10k
        x = T.("RAW_INCOME_CPI"); if ~isnumeric(x), x = toNum(x); end
        x = x / 10000;  % per $10,000

        % Outcomes
        yND = T.(ndviVar); if ~isnumeric(yND), yND = toNum(yND); end
        yLS = T.(lstVar);  if ~isnumeric(yLS), yLS = toNum(yLS); end

        % Overlap filter (supports either 0–100 or 0–1 scales)
        if any(V=="PCT_OVERLAP")
            ov = T.("PCT_OVERLAP"); if ~isnumeric(ov), ov = toNum(ov); end
            thr = OverlapThresh;
            if all(ov<=1 | isnan(ov)), thr = OverlapThresh/100; end
            keepOv = ov >= thr;
        else
            keepOv = true(size(x));
        end

        % Valid (finite) pairs only after overlap filter
        okND = isfinite(x) & isfinite(yND) & keepOv;
        okLS = isfinite(x) & isfinite(yLS) & keepOv;

        % Require enough tracts and variation in X
        if nnz(okND) >= MinTracts && std(x(okND),'omitnan') > 0
            xND = x(okND); yNDv = yND(okND);
            if useWinsor
                p = prctile(xND, winz);  xND  = min(max(xND, p(1)), p(2));
                p = prctile(yNDv, winz); yNDv = min(max(yNDv, p(1)), p(2));
            end
            mdlND = fitlm(xND, yNDv);
            sND(yy) = mdlND.Coefficients.Estimate(2); % NDVI per $10k
        end

        if nnz(okLS) >= MinTracts && std(x(okLS),'omitnan') > 0
            xLS = x(okLS); yLSv = yLS(okLS);
            if useWinsor
                p = prctile(xLS, winz);  xLS  = min(max(xLS, p(1)), p(2));
                p = prctile(yLSv, winz); yLSv = min(max(yLSv, p(1)), p(2));
            end
            mdlLS = fitlm(xLS, yLSv);
            sLS(yy) = mdlLS.Coefficients.Estimate(2); % °C per $10k
        end
    end

    % keep only fully observed cities (both metrics in both years)
    if all(isfinite(sND)) && all(isfinite(sLS))
        CityKeep(end+1,1) = city;                 %#ok<AGROW>
        NDVI_A(end+1,1)   = sND(1); NDVI_B(end+1,1) = sND(2); %#ok<AGROW>
        LST_A (end+1,1)   = sLS(1); LST_B (end+1,1) = sLS(2); %#ok<AGROW>
    end
end

%% ----- Pretty city labels for publication -----
% Explicit map for special/ambiguous names and those without a trailing state code
cityMap = containers.Map( ...
{
% ---- your condensed input names ----
'Albuquerque','Atlanta','Bakersfield','Baltimore','BangorME','Billings','Birmingham','Bismarck','Bloomington','Boise','BostonMA','Buffalo','BurlingtonVT','Casper','CharlestonSC','CharlestonWV','Chicago','Cincinnati','ClearwaterFL','Cleveland','ColumbusOH','CorpusChristi','Dallas','Denver','DesMoinesIA','Detroit','ElPaso','FargoND','FortCollins','Fresno','GrandJunction','GreenBay','HartfordCT','HoustonTX','IdahoFallsID','JacksonMS','KansasCity','LakeHavasuCityAZ','LasCrucesNM','LasVegas','LexingtonKY','LittleRock','LosAngeles','LouisvilleKY','LubbockTXuse','Manchester','Memphis','Miami','Milwaukee','Missoula','NYC_BOUNDS','Nashville','NewOrleans','NorfolkVA','OklahomaCity','Omaha','Pensacola','Philadelphia','Phoenix','PittsburghPA','PortlandME','PortlandOR','ProvidenceRI','Raleigh','Reno','Riverside','Sacramento','SaltLakeCity','SanAntonio','SanDiego','SanFrancisco','SantaFeNM','SavannahGA','SeattleWA','ShreveportLA','SiouxFalls','Spokane','StLouis','Syracuse','Tacoma','Tallahassee','ThousandOaks','Topeka','Tucson','TwinCities','WashingtonDC','Wichita'
}, ...
{
% ---- pretty printed "City, ST" ----
'Albuquerque, NM','Atlanta, GA','Bakersfield, CA','Baltimore, MD','Bangor, ME','Billings, MT','Birmingham, AL','Bismarck, ND','Bloomington, IN','Boise, ID','Boston, MA','Buffalo, NY','Burlington, VT','Casper, WY','Charleston, SC','Charleston, WV','Chicago, IL','Cincinnati, OH','Clearwater, FL','Cleveland, OH','Columbus, OH','Corpus Christi, TX','Dallas, TX','Denver, CO','Des Moines, IA','Detroit, MI','El Paso, TX','Fargo, ND','Fort Collins, CO','Fresno, CA','Grand Junction, CO','Green Bay, WI','Hartford, CT','Houston, TX','Idaho Falls, ID','Jackson, MS','Kansas City, MO','Lake Havasu City, AZ','Las Cruces, NM','Las Vegas, NV','Lexington, KY','Little Rock, AR','Los Angeles, CA','Louisville, KY','Lubbock, TX','Manchester, NH','Memphis, TN','Miami, FL','Milwaukee, WI','Missoula, MT','New York City, NY','Nashville, TN','New Orleans, LA','Norfolk, VA','Oklahoma City, OK','Omaha, NE','Pensacola, FL','Philadelphia, PA','Phoenix, AZ','Pittsburgh, PA','Portland, ME','Portland, OR','Providence, RI','Raleigh, NC','Reno, NV','Riverside, CA','Sacramento, CA','Salt Lake City, UT','San Antonio, TX','San Diego, CA','San Francisco, CA','Santa Fe, NM','Savannah, GA','Seattle, WA','Shreveport, LA','Sioux Falls, SD','Spokane, WA','St. Louis, MO','Syracuse, NY','Tacoma, WA','Tallahassee, FL','Thousand Oaks, CA','Topeka, KS','Tucson, AZ','Twin Cities, MN','Washington, DC','Wichita, KS'
});

% Optional convenience for a few you listed that are already clean:
% 'Phoenix' => 'Phoenix, AZ', 'Miami' => 'Miami, FL', etc. are covered above.

% --- Helper: insert spaces before interior capitals (FortCollins -> Fort Collins)
camel2space = @(s) regexprep(s, '(?<=\p{Ll})(\p{Lu})', ' $1');

% --- Helper: pretty fallback when cityMap lacks a key
function out = prettyFallback(in)
    s = string(in);
    s = replace(s, "_", " ");  % underscores to spaces
    s = regexprep(s, 'StLouis$', 'St. LouisMO'); % handled by map already, but just in case
    % If trailing 2-letter state exists, split it off
    m = regexp(s, '^(.+?)([A-Z]{2})$', 'tokens','once');
    if ~isempty(m)
        city = strtrim(camel2space(m{1}));
        st   = m{2};
        % Common multi-word fixes
        city = strrep(city,'Fort Collins','Fort Collins');
        city = regexprep(city,'\s+',' ');
        out = city + ", " + st;
    else
        % No trailing state: just space CamelCase and leave as-is
        out = camel2space(s);
    end
end

% --- Helper: get pretty label via map with fallback
function out = prettyCityLabel(rawName, cityMap)
    key = char(string(rawName));
    if isKey(cityMap, key)
        out = string(cityMap(key));
    else
        out = prettyFallback(rawName);
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

if isempty(Tplot)
    error('No cities produced valid slopes for BOTH years. Check RAW_INCOME_CPI, %s_NDVI and %s_LST, overlap, and variation.', metricLabel, metricLabel);
end

% Sort by earlier NDVI slope (descending) and set y positions
[~, ord] = sort(Tplot.NDVI_A, 'descend');
Tplot = Tplot(ord,:);
nCities = height(Tplot);
y = (nCities:-1:1)';

% Labels (prettify underscores)
labelsLeft = arrayfun(@(s) prettyCityLabel(s, cityMap), flip(Tplot.City));

% ------------ Draw figure ------------
fig = figure('Color','w','Units','pixels','Position',[100 100 1200 800]);

% (a) Income → NDVI
subplot(1,2,1); hold on;
for i = 1:nCities
    plot([Tplot.NDVI_A(i) Tplot.NDVI_B(i)], [y(i) y(i)], '-', 'Color',[0.7 0.7 0.7], 'LineWidth',1);
    plot(Tplot.NDVI_A(i), y(i), 'o', 'MarkerFaceColor',[0 0.45 0.74], 'MarkerEdgeColor','k', 'MarkerSize',6);  % yearA (blue)
    plot(Tplot.NDVI_B(i), y(i), 'o', 'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','k', 'MarkerSize',6); % yearB (orange)
end
yticks(1:nCities); yticklabels(labelsLeft);
xlabel(sprintf('Income–%s NDVI (per $10{,}000)', metricLabel));
title(sprintf('(a) Income–NDVI slope (%d vs %d)', yearA, yearB));
xline(0,'k-');

ndVals = [Tplot.NDVI_A; Tplot.NDVI_B];
if all(~isfinite(ndVals))
    xlim([-0.1 0.1]);
else
    mn = min(ndVals,[],'omitnan'); mx = max(ndVals,[],'omitnan');
    pad = (mn==mx) * 0.1 + (mn~=mx) * 0.01;
    if mn==mx, mn = mn - pad; mx = mx + pad; end
    xlim([mn - pad, mx + pad]);
end
ylim([0.5 nCities+0.5]); set(gca,'FontSize',8,'Box','off');

% legend (point-only)
hA = plot(NaN,NaN,'o','MarkerFaceColor',[0 0.45 0.74],'MarkerEdgeColor','k','MarkerSize',6);
hB = plot(NaN,NaN,'o','MarkerFaceColor',[0.85 0.33 0.10],'MarkerEdgeColor','k','MarkerSize',6);
legend([hA hB], {num2str(yearA), num2str(yearB)}, 'Location','best'); legend boxoff;
hold off;

% (b) Income → LST
subplot(1,2,2); hold on;
for i = 1:nCities
    plot([Tplot.LST_A(i) Tplot.LST_B(i)], [y(i) y(i)], '-', 'Color',[0.7 0.7 0.7], 'LineWidth',1);
    plot(Tplot.LST_A(i), y(i), 'o', 'MarkerFaceColor',[0 0.45 0.74], 'MarkerEdgeColor','k', 'MarkerSize',6);   % yearA (blue)
    plot(Tplot.LST_B(i), y(i), 'o', 'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','k', 'MarkerSize',6); % yearB (orange)
end
yticks(1:nCities); yticklabels([]);
xlabel(sprintf('Income–%s LST (\\circC per $10{,}000)', metricLabel));  % TeX degree symbol
title(sprintf('(b) Income–LST slope (%d vs %d)', yearA, yearB));
xline(0,'k-');

lsVals = [Tplot.LST_A; Tplot.LST_B];
if all(~isfinite(lsVals))
    xlim([-0.5 0.5]);
else
    mn = min(lsVals,[],'omitnan'); mx = max(lsVals,[],'omitnan');
    pad = (mn==mx) * 0.5 + (mn~=mx) * 0.1;
    if mn==mx, mn = mn - pad; mx = mx + pad; end
    xlim([mn - pad, mx + pad]);
end
ylim([0.5 nCities+0.5]); set(gca,'FontSize',8,'Box','off'); 
hold off;

fprintf(['Dumbbell plot includes %d cities with valid %d and %d slopes for both metrics ', ...
         '(RAW\\_INCOME per $10k; %s_* outcomes; PCT\\_OVERLAP >= %g; MinTracts=%d).\n'], ...
        nCities, yearA, yearB, upper(metricMode), OverlapThresh, MinTracts);

% optional save
if ~isempty(savePNG)
    exportgraphics(fig, savePNG, 'Resolution',300);
    fprintf('Saved figure to %s\n', savePNG);
end