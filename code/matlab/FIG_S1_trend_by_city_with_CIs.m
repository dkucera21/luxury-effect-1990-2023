%% FIG_S1_trend_by_city_with_CIs.m
% ------------------------------------------------------------------------------
% Robustness figure: City-specific trend of the luxury effect (all years, with CIs)
% For each City×Year, fit a tract-level slope (y ~ RAW_INCOME_CPI per $10k).
% Then, for each city, regress those yearly slopes on time to get a per-year
% trend in the luxury effect with 95% CI. Plots NDVI and LST side-by-side.
%
% Requirements in workspace:
%   - CENSUS_TABLES_rebuilt : struct of tract tables T_<YEAR>_<CITY> (or T_<YEAR>n_<CITY>)
%   - CityListMaster / MasterCities.City / CityNames (any one)
% ------------------------------------------------------------------------------

%% ===================== USER SETTINGS =====================
yearsUse     = [1990 2000 2010 2020 2023];
metricMode   = 'MEAN';    % 'MEAN' or 'MEDIAN'
OverlapThresh = 50;         % PCT_OVERLAP cutoff (supports 0–100 or 0–1 in data)
MinTracts     = 5;          % minimum valid tracts per city–year to fit slope
useWinsor     = false;      % optional robustness
winzPct       = [1 99];     % winsor limits if useWinsor = true

excludeCities = lower(string({'Honolulu','Anchorage'})); % optional exclude
savePNG       = '';          % e.g., 'figS_trend_by_city_CIs.png' ('' = don’t save)
%% =========================================================

assert(exist('CENSUS_TABLES_rebuilt','var')==1 && isstruct(CENSUS_TABLES_rebuilt), ...
    'CENSUS_TABLES_rebuilt not found.');

% ---------- city list ----------
if exist('CityListMaster','var') && ~isempty(CityListMaster)
    CityList = string(CityListMaster(:));
elseif exist('MasterCities','var') && istable(MasterCities) && ismember('City', MasterCities.Properties.VariableNames)
    CityList = string(MasterCities.City);
elseif exist('CityNames','var') && ~isempty(CityNames)
    CityList = string(CityNames(:));
else
    error('No valid city list found (CityListMaster, MasterCities.City, or CityNames).');
end
CityList = unique(strtrim(CityList));
tblNames = string(fieldnames(CENSUS_TABLES_rebuilt));
tonum    = @(x) str2double(regexprep(string(x),'[,\$%]',''));

% ---------- metric columns ----------
metricMode = upper(string(metricMode));
switch metricMode
    case "MEAN"
        ndviVar = "MEAN_NDVI";  lstVar = "MEAN_LST";  metricLabel = "MEAN";
    case "MEDIAN"
        ndviVar = "MEDIAN_NDVI"; lstVar = "MEDIAN_LST"; metricLabel = "MEDIAN";
    otherwise
        error('metricMode must be ''MEAN'' or ''MEDIAN''.');
end

%% ----- Pretty city labels for publication (define BEFORE plotting) -----
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

% ---------- build per-city trends for NDVI and LST ----------
T_ndvi = compute_city_trends(CityList, yearsUse, ndviVar, ...
                             CENSUS_TABLES_rebuilt, tblNames, ...
                             OverlapThresh, MinTracts, useWinsor, winzPct, tonum);
T_lst  = compute_city_trends(CityList, yearsUse, lstVar,  ...
                             CENSUS_TABLES_rebuilt, tblNames, ...
                             OverlapThresh, MinTracts, useWinsor, winzPct, tonum);

% optional excludes
if ~isempty(excludeCities)
    maskN = false(height(T_ndvi),1);
    maskL = false(height(T_lst ),1);
    for ex = excludeCities(:)'
        maskN = maskN | contains(lower(T_ndvi.City), ex);
        maskL = maskL | contains(lower(T_lst.City ), ex);
    end
    T_ndvi(maskN,:) = [];
    T_lst(maskL,:)  = [];
end

% ---------- plotting (CI dumbbells) ----------
f = figure('Color','w','Units','pixels','Position',[120 120 1200 800]);

% (a) NDVI trend (per year)
subplot(1,2,1); hold on; box off;
Tn = sortrows(T_ndvi, 'Slope', 'descend');
y  = (height(Tn):-1:1)';

for i = 1:height(Tn)
    plot([Tn.CI_L(i) Tn.CI_U(i)], [y(i) y(i)], '-', 'Color',[0.7 0.7 0.7], 'LineWidth', 2); % CI
    plot(Tn.Slope(i), y(i), 'o', 'MarkerFaceColor',[0 0.45 0.74], 'MarkerEdgeColor','k', 'MarkerSize', 6);
end
yticks(1:height(Tn));
yticklabels(pretty_labels(flip(Tn.City), cityMap));
xlabel(sprintf('Change in income–%s NDVI per year (per $10k·yr^{-1})', metricLabel));
title(sprintf('(a) City-specific trend of luxury effect: NDVI (%d–%d, all years)', min(yearsUse), max(yearsUse)));
xline(0,'k-');
ylim([0.5 height(Tn)+0.5]);
set(gca,'FontSize',8,'TickDir','out');

% (b) LST trend (per year)
subplot(1,2,2); hold on; box off;
Tl = sortrows(T_lst, 'Slope', 'descend');
y2 = (height(Tl):-1:1)';

for i = 1:height(Tl)
    plot([Tl.CI_L(i) Tl.CI_U(i)], [y2(i) y2(i)], '-', 'Color',[0.7 0.7 0.7], 'LineWidth', 2); % CI
    plot(Tl.Slope(i), y2(i), 'o', 'MarkerFaceColor',[0.85 0.33 0.10], 'MarkerEdgeColor','k', 'MarkerSize', 6);
end
yticks(1:height(Tl)); yticklabels([]);  % left panel has labels
xlabel(sprintf('Change in income–%s LST per year (°C per $10k·yr^{-1})', metricLabel));
title(sprintf('(b) City-specific trend of luxury effect: LST (%d–%d, all years)', min(yearsUse), max(yearsUse)));
xline(0,'k-');
ylim([0.5 height(Tl)+0.5]);
set(gca,'FontSize',8,'TickDir','out');

% announce and optionally save
fprintf('Built city-specific trends using all available years: NDVI n=%d cities, LST n=%d cities.\n', height(Tn), height(Tl));
if ~isempty(savePNG)
    exportgraphics(f, savePNG, 'Resolution', 300);
    fprintf('Saved figure to %s\n', savePNG);
end

% -------- expose tables in workspace for inspection --------
assignin('base','CityTrends_NDVI', T_ndvi);
assignin('base','CityTrends_LST',  T_lst);

%% ===================== local functions =====================

function Tout = compute_city_trends(CityList, yearsUse, yVar, CENSUS, tblNames, OverlapThresh, MinTracts, useWinsor, winzPct, tonum)
    % Build city-year luxury slopes, then city OLS trend with 95% CI.
    pctNames = {'PCT_OVERLAP','PCT_COVER','PCT_COV','PercentAreaNDVI'};
    CY = table('Size',[0 4],'VariableTypes',{'string','double','double','double'}, ...
               'VariableNames',{'City','Year','LE','nTracts'});

    for ci = 1:numel(CityList)
        city = string(CityList(ci));
        for yr = yearsUse(:)'
            % T_YEAR_city or T_YEARn_city
            t1 = sprintf('T_%d_%s',  yr, city);
            t2 = sprintf('T_%dn_%s', yr, city);
            if ismember(t1, tblNames)
                T = CENSUS.(t1);
            elseif ismember(t2, tblNames)
                T = CENSUS.(t2);
            else
                continue
            end
            if ~istable(T) || height(T)==0, continue; end
            V = string(T.Properties.VariableNames);

            % guards: RAW_INCOME_CPI and outcome present
            if ~all(ismember(["RAW_INCOME_CPI", yVar], V)), continue; end

            % overlap (0–100 or 0–1)
            ov = []; vPct = first_match(V, pctNames);
            if vPct ~= ""
                ov = tonum(T.(vPct));
            end
            if ~isempty(ov)
                thr = OverlapThresh; if all(ov<=1 | isnan(ov)), thr = OverlapThresh/100; end
                keepOv = ov >= thr;
            else
                keepOv = true(height(T),1);
            end

            % income per $10k
            x = tonum(T.("RAW_INCOME_CPI")) ./ 10000;
            y = tonum(T.(yVar));

            ok = isfinite(x) & isfinite(y) & keepOv;
            if nnz(ok) < MinTracts || std(x(ok),'omitnan')==0, continue; end
            x = x(ok); y = y(ok);
            if useWinsor
                p = prctile(x, winzPct);  x = min(max(x, p(1)), p(2));
                p = prctile(y, winzPct);  y = min(max(y, p(1)), p(2));
            end

            % tract-level slope within this city-year
            mdl = fitlm(x, y);                     % y ~ x
            b   = mdl.Coefficients.Estimate(2);    % LE (per $10k)
            CY  = [CY; {char(city), double(yr), b, numel(x)}]; %#ok<AGROW>
        end
    end

    % regress city-year LE on time (per-year trend) with 95% CI
    if isempty(CY)
        Tout = table('Size',[0 6],'VariableTypes',{'string','double','double','double','double','double'}, ...
                     'VariableNames',{'City','Slope','SE','CI_L','CI_U','nYears'});
        return
    end
    CY.time = CY.Year - mean(CY.Year,'omitnan');   % center calendar year

    cities = unique(string(CY.City),'stable');
    Slope = NaN(numel(cities),1); SE = NaN(numel(cities),1);
    CI_L  = NaN(numel(cities),1); CI_U = NaN(numel(cities),1);
    nY    = NaN(numel(cities),1);

    for i = 1:numel(cities)
        ci = cities(i);
        D  = CY(string(CY.City)==ci, {'LE','time','Year'});
        yrs = unique(D.Year);
        if numel(yrs) < 3, continue; end
        mdl = fitlm(D, 'LE ~ time');   % LE per year (centered time)
        k   = strcmp(mdl.Coefficients.Properties.RowNames, 'time');
        if any(k)
            Slope(i) = mdl.Coefficients.Estimate(k);     % per year
            SE(i)    = mdl.Coefficients.SE(k);
            ci2      = coefCI(mdl);                      % 95% CI
            CI_L(i)  = ci2(k,1);
            CI_U(i)  = ci2(k,2);
            nY(i)    = numel(yrs);
        end
    end

    Tout = table(cities, Slope, SE, CI_L, CI_U, nY, ...
        'VariableNames', {'City','Slope','SE','CI_L','CI_U','nYears'});

    % drop cities without estimates
    Tout = Tout(isfinite(Tout.Slope), :);
end

function v = first_match(V, cands)
    v = "";
    VL = lower(string(V));
    for i = 1:numel(cands)
        j = find(VL == lower(string(cands{i})), 1, 'first');
        if ~isempty(j), v = string(V(j)); return; end
    end
end

function L = pretty_labels(cities, cityMap)
% Map raw city tokens to "City, ST" using containers.Map cityMap.
% Falls back to the original name if a key is missing.
    L = cellstr(cities);
    if exist('cityMap','var')==1 && isa(cityMap,'containers.Map')
        for i = 1:numel(L)
            k = L{i};
            if isKey(cityMap, k)
                L{i} = cityMap(k);
            end
        end
    end
end
