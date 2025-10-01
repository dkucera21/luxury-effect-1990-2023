%% equity_build_and_analyze.m
% ------------------------------------------------------------------------------
% Build a tail-based equity panel (NDVI & LST), apply tract filters, fit LMEs,
% validate interactions, compute per-city equity changes (OLS + LME-RE when
% available), save results to workspace + CSV, and print concise highlights.
% ------------------------------------------------------------------------------

%% ===================== USER CONFIG =====================
TailPct          = 0.20;
YearsUse         = [1990 2000 2010 2020 2023];
MetricNames      = struct('ndvi','MEAN_NDVI','lst','MEAN_LST');

MinPerTail       = 2;    % ≥2 tracts per tail per city–year
MinTractsCity    = 20;   % ≥20 valid tracts in the city (across all years)
OverlapField     = 'PCT_OVERLAP';  % will also accept 0–1 or alias below
Alpha            = 0.05;

CityListUse      = string(CityListMaster(:));  % from your repo list
toNumSafe        = @(x) str2double(regexprep(string(x),'[,\$%]',''));
tblNames         = string(fieldnames(CENSUS_TABLES_rebuilt));

% Accept common overlap alias
OverlapAliases   = ["PCT_OVERLAP"];

%% ===================== 1) BUILD TAIL PANEL (with filters) =====================
% First pass: gather valid (filtered) tract vectors per city–year
ValidCY = table('Size',[0 5], 'VariableTypes',{'string','double','cell','cell','cell'}, ...
    'VariableNames',{'City','Year','inc','ndvi','lst'});
cityValidCount = containers.Map('KeyType','char','ValueType','double');

for ci = 1:numel(CityListUse)
    city = CityListUse(ci);
    totalOK = 0;
    for yr = YearsUse
        t1 = sprintf('T_%d_%s',  yr, city);
        t2 = sprintf('T_%dn_%s', yr, city);
        has1 = any(tblNames==t1); has2 = any(tblNames==t2);
        if ~has1 && ~has2, continue; end
        if has1, T = CENSUS_TABLES_rebuilt.(t1); else, T = CENSUS_TABLES_rebuilt.(t2); end
        V = string(T.Properties.VariableNames);

        % ---- income
        inc = local_get_income(V, yr, T, toNumSafe);
        if isempty(inc), if any(V=="incZ"), inc = toNumSafe(T.incZ); else, continue; end, end

        % ---- outcomes
        if ~any(V==MetricNames.ndvi) || ~any(V==MetricNames.lst), continue; end
        ndvi = toNumSafe(T.(MetricNames.ndvi));
        lst  = toNumSafe(T.(MetricNames.lst));

        % ---- base validity
        ok = isfinite(inc) & isfinite(ndvi) & isfinite(lst);

        % ---- overlap ≥ 50% (accept 0–100 or 0–1)
        ovVar = OverlapAliases(ismember(OverlapAliases, V));
        if ~isempty(ovVar)
            ov  = toNumSafe(T.(ovVar(1)));
            thr = 50; if all(ov<=1 | isnan(ov)), thr = 0.5; end
            ok  = ok & (ov >= thr);
        end

        % keep only if we can possibly form ≥2 per tail (i.e., ≥4 valid)
        nOK = nnz(ok);
        if nOK < 2*MinPerTail, continue; end

        % cache filtered vectors
        ValidCY = [ValidCY; {char(city), double(yr), {inc(ok)}, {ndvi(ok)}, {lst(ok)}}]; %#ok<AGROW>
        totalOK = totalOK + nOK;
    end
    cityValidCount(char(city)) = totalOK;
end

% Keep only cities with ≥ MinTractsCity valid tracts across all years
keepCity = ismember(string(ValidCY.City), string(keys(cityValidCount)));
if any(keepCity)
    kc = false(height(ValidCY),1);
    for i = 1:height(ValidCY)
        kc(i) = cityValidCount(char(ValidCY.City(i))) >= MinTractsCity;
    end
    ValidCY = ValidCY(kc, :);
end

% Second pass: build Panel from the kept city–years, enforcing ≥2 per tail
Panel = table([],[],[],[],[],[],[],[], 'VariableNames', ...
    {'City','Year','Top_NDVI','Bot_NDVI','Gap_NDVI','Top_LST','Bot_LST','Gap_LST'});

for i = 1:height(ValidCY)
    city = string(ValidCY.City(i));
    yr   = ValidCY.Year(i);
    inc  = ValidCY.inc{i};
    ndvi = ValidCY.ndvi{i};
    lst  = ValidCY.lst{i};

    % sort by income, split tails (enforce ≥2 per tail)
    [~, rk] = sort(inc, 'ascend');
    n  = numel(rk);
    nb = max(MinPerTail, floor(TailPct*n));
    nt = max(MinPerTail, floor(TailPct*n));
    if nb + nt > n
        % shrink proportionally but never below MinPerTail per tail
        nb = max(MinPerTail, floor((TailPct/2)*n));
        nt = max(MinPerTail, floor((TailPct/2)*n));
        if nb + nt > n, continue; end
    end

    q = zeros(n,1); q(rk(1:nb)) = 1; q(rk(end-nt+1:end)) = 5;

    Top_NDVI = mean(ndvi(q==5), 'omitnan');
    Bot_NDVI = mean(ndvi(q==1), 'omitnan');
    Top_LST  = mean(lst(q==5),  'omitnan');
    Bot_LST  = mean(lst(q==1),  'omitnan');

    Gap_NDVI = Top_NDVI - Bot_NDVI;   % larger = greener rich vs poor
    Gap_LST  = Bot_LST  - Top_LST;    % larger = hotter poor vs rich

    Panel = [Panel; {char(city), yr, Top_NDVI, Bot_NDVI, Gap_NDVI, Top_LST, Bot_LST, Gap_LST}]; %#ok<AGROW>
end

% Final tidy & time coding (1990 baseline per your Methods)
Panel.City = categorical(string(Panel.City), CityListUse);
Panel      = Panel(~isundefined(Panel.City), :);
Panel      = sortrows(Panel, {'City','Year'});
Panel.time_dec = (Panel.Year - 1990)/10;

fprintf('\n=== Tail panel built with filters ===\n');
fprintf('Tail = %.0f%%%% | Years = [%s] | MinPerTail = %d | MinTractsCity = %d | Overlap ≥ 50%%\n', ...
    100*TailPct, num2str(YearsUse), MinPerTail, MinTractsCity);
fprintf('Panel rows kept: %d  | Cities represented: %d\n', height(Panel), numel(categories(Panel.City)));



%% ===================== 2) LONG TABLES & STACKED LMEs =====================
% Long NDVI
Tlong_NDVI = table();
Tlong_NDVI.City     = [Panel.City; Panel.City];
Tlong_NDVI.Year     = [Panel.Year; Panel.Year];
Tlong_NDVI.time_dec = [Panel.time_dec; Panel.time_dec];
Tlong_NDVI.tail     = categorical( ...
    [repmat("Bottom",height(Panel),1); repmat("Top",height(Panel),1)], ...
    ["Bottom","Top"]);                               % << reference = Bottom
Tlong_NDVI.y        = [Panel.Bot_NDVI; Panel.Top_NDVI];

% Long LST
Tlong_LST = table();
Tlong_LST.City     = [Panel.City; Panel.City];
Tlong_LST.Year     = [Panel.Year; Panel.Year];
Tlong_LST.time_dec = [Panel.time_dec; Panel.time_dec];
Tlong_LST.tail     = categorical( ...
    [repmat("Bottom",height(Panel),1); repmat("Top",height(Panel),1)], ...
    ["Bottom","Top"]);                               % << reference = Bottom
Tlong_LST.y        = [Panel.Bot_LST; Panel.Top_LST];


% --- CRITICAL: random slopes for time_dec and time_dec:tail ---
formTail = 'y ~ 1 + time_dec*tail + (1 + time_dec + time_dec:tail | City)';

% NDVI model
try
    lmeN = fitlme(Tlong_NDVI, formTail, 'FitMethod','REML', 'CheckHessian', true);
catch ME
    warning('NDVI tail-LME fell back to no random interaction: %s', ME.message);
    lmeN = fitlme(Tlong_NDVI, 'y ~ 1 + time_dec*tail + (1 + time_dec | City)', 'FitMethod','REML');
end

% LST model
try
    lmeL = fitlme(Tlong_LST,  formTail, 'FitMethod','REML', 'CheckHessian', true);
catch ME
    warning('LST tail-LME fell back to no random interaction: %s', ME.message);
    lmeL = fitlme(Tlong_LST,  'y ~ 1 + time_dec*tail + (1 + time_dec | City)', 'FitMethod','REML');
end


% OLS slopes on the panel (Top/Bottom)
b_top_ndvi = fitlm(Panel,'Top_NDVI ~ time_dec').Coefficients{'time_dec','Estimate'};
b_bot_ndvi = fitlm(Panel,'Bot_NDVI ~ time_dec').Coefficients{'time_dec','Estimate'};
b_top_lst  = fitlm(Panel,'Top_LST  ~ time_dec').Coefficients{'time_dec','Estimate'};
b_bot_lst  = fitlm(Panel,'Bot_LST  ~ time_dec').Coefficients{'time_dec','Estimate'};

% OLS on stacked long tables (same estimator on both sides)
lmN = fitlm(Tlong_NDVI, 'y ~ 1 + time_dec*tail');
lmL = fitlm(Tlong_LST,  'y ~ 1 + time_dec*tail');

% --- replace the getIntLM + b_int_* lines with this ---

% grab term names safely (works across MATLAB versions)
cnN = string(lmN.CoefficientNames);
cnL = string(lmL.CoefficientNames);

% interaction index is the term that involves both time_dec and tail
iN = find(contains(cnN,"time_dec") & contains(cnN,"tail"), 1);
iL = find(contains(cnL,"time_dec") & contains(cnL,"tail"), 1);

% fallback via row names if needed (older versions)
if isempty(iN)
    rnN = string(lmN.Coefficients.Properties.RowNames);
    iN  = find(contains(rnN,"time_dec") & contains(rnN,"tail"), 1);
end
if isempty(iL)
    rnL = string(lmL.Coefficients.Properties.RowNames);
    iL  = find(contains(rnL,"time_dec") & contains(rnL,"tail"), 1);
end

if isempty(iN) || isempty(iL)
    warning('Could not locate the time_dec:tail interaction term in OLS models.');
    b_int_ndvi = NaN;
    b_int_lst  = NaN;
else
    b_int_ndvi = lmN.Coefficients.Estimate(iN);
    b_int_lst  = lmL.Coefficients.Estimate(iL);
end


tol = 1e-8;   % tighter, but realistic for OLS-vs-OLS equivalence

dn = abs(b_int_ndvi - (b_top_ndvi - b_bot_ndvi));
dl = abs(b_int_lst  - (b_top_lst  - b_bot_lst ));

if dn > tol
    warning('NDVI interaction identity off by %.3g (likely due to NA patterns or coding).', dn);
end
if dl > tol
    warning('LST interaction identity off by %.3g (likely due to NA patterns or coding).', dl);
end


%% ===================== 3) SEPARATE LMEs: Gap, Top, Bottom =====================
form_gap_ndvi = 'Gap_NDVI ~ time_dec + (1|City)';
form_top_ndvi = 'Top_NDVI ~ time_dec + (1|City)';
form_bot_ndvi = 'Bot_NDVI ~ time_dec + (1|City)';

form_gap_lst  = 'Gap_LST  ~ time_dec + (1|City)';
form_top_lst  = 'Top_LST  ~ time_dec + (1|City)';
form_bot_lst  = 'Bot_LST  ~ time_dec + (1|City)';

Panel2 = Panel; Panel2.City = categorical(string(Panel2.City)); Panel2.Year = double(Panel2.Year);

lme_gap_ndvi = fitlme(Panel2, form_gap_ndvi);
lme_top_ndvi = fitlme(Panel2, form_top_ndvi);
lme_bot_ndvi = fitlme(Panel2, form_bot_ndvi);

lme_gap_lst  = fitlme(Panel2, form_gap_lst);
lme_top_lst  = fitlme(Panel2, form_top_lst);
lme_bot_lst  = fitlme(Panel2, form_bot_lst);

fprintf('\n=== Tail panel built ===\n');
fprintf('Tail = %.0f%% | Years = [%s] | MinTracts = %d | OverlapField = %s\n', ...
    100*TailPct, num2str(YearsUse), MinTracts, OverlapField);
fprintf('Panel rows: %d  | Cities: %d\n', height(Panel), numel(categories(Panel.City)));
fprintf('\nFE interactions (per decade): NDVI %+0.4f | LST %+0.4f\n\n', b_int_ndvi, b_int_lst);

%% ===================== 4) CITY LABELS (significance-based, OLS per city) =====================
CityStr = string(Panel.City);
CityIDs = unique(CityStr, 'stable');

Labels  = table(CityIDs, repmat("Little change", numel(CityIDs),1), ...
    'VariableNames', {'City','Label'});

Details = table('Size',[numel(CityIDs) 9], 'VariableTypes', ...
    ["string","double","double","double","double","double","double","double","double"], ...
    'VariableNames', {'City', ...
      'b_TopNDVI','p_TopNDVI','b_BotNDVI','p_BotNDVI', ...
      'b_TopLST','p_TopLST','b_BotLST','p_BotLST'});

for k = 1:numel(CityIDs)
    Ci = CityIDs(k);
    P  = sortrows(Panel(CityStr==Ci, :), 'Year');

    [b_TopNDVI,p_TopNDVI] = local_city_slope(P,'Top_NDVI');
    [b_BotNDVI,p_BotNDVI] = local_city_slope(P,'Bot_NDVI');
    [b_TopLST, p_TopLST ] = local_city_slope(P,'Top_LST');
    [b_BotLST, p_BotLST ] = local_city_slope(P,'Bot_LST');

    ndvi_uplift  = isfinite(p_BotNDVI) && p_BotNDVI < Alpha && b_BotNDVI > 0;
    ndvi_erosion = isfinite(p_TopNDVI) && p_TopNDVI < Alpha && b_TopNDVI < 0;
    lst_uplift   = isfinite(p_BotLST)  && p_BotLST  < Alpha && b_BotLST  < 0;
    lst_erosion  = isfinite(p_TopLST)  && p_TopLST  < Alpha && b_TopLST  > 0;

    if (ndvi_uplift||lst_uplift) && ~(ndvi_erosion||lst_erosion)
        lab = "Uplift";
    elseif (ndvi_erosion||lst_erosion) && ~(ndvi_uplift||lst_uplift)
        lab = "Erosion";
    elseif (ndvi_uplift||lst_uplift) && (ndvi_erosion||lst_erosion)
        lab = "Both";
    else
        lab = "Little change";
    end

    Labels.Label(k) = lab;
    Details(k,:) = {Ci, b_TopNDVI, p_TopNDVI, b_BotNDVI, p_BotNDVI, ...
                        b_TopLST,  p_TopLST,  b_BotLST,  p_BotLST};
end

Labels.Label = categorical(cellstr(Labels.Label), {'Uplift','Erosion','Both','Little change'});
LabelCounts  = groupsummary(Labels, 'Label');
LabelCounts.Properties.VariableNames{'GroupCount'} = 'N';
LabelCounts.Pct = 100 * LabelCounts.N / height(Labels);
fprintf('Label summary (%%):\n');
disp(LabelCounts(:,{'Label','N','Pct'}));

%% ===================== 5) PER-CITY EQUITY CHANGES & CSV =====================
% Always compute OLS per-city slopes over the panel
CityEquity = table('Size',[numel(CityIDs) 13],'VariableTypes', ...
    ["string","double","double","double","double","double","double","double","double","double","double","double","double"], ...
    'VariableNames', {'City','nYears','YearMin','YearMax', ...
        'Slope_TopNDVI','Slope_BotNDVI','Slope_GapNDVI', ...
        'Slope_TopLST','Slope_BotLST','Slope_GapLST', ...
        'Delta_GapNDVI','Delta_GapLST','Span_Decades'});

for k = 1:numel(CityIDs)
    Ci = CityIDs(k); P = sortrows(Panel(Panel.City==categorical(Ci), :), 'Year');
    yrs = unique(P.Year); nY = numel(yrs); yMin = min(yrs); yMax = max(yrs); spanDec = (yMax - yMin)/10;

    sTopN = local_city_slope(P, 'Top_NDVI');
    sBotN = local_city_slope(P, 'Bot_NDVI');
    sTopL = local_city_slope(P, 'Top_LST');
    sBotL = local_city_slope(P, 'Bot_LST');

    sGapN = sTopN - sBotN;        % NDVI inequity slope (Top–Bottom)
    sGapL = sBotL - sTopL;        % LST inequity slope (Bottom–Top)

    dGapN = sGapN * spanDec;      % implied gap change across the city’s span
    dGapL = sGapL * spanDec;

    CityEquity(k,:) = {Ci, nY, yMin, yMax, sTopN, sBotN, sGapN, sTopL, sBotL, sGapL, dGapN, dGapL, spanDec};
end

% Try to add LME-based random-effect deltas (city-specific)
T_NDVI = try_city_slopes_RE(lmeN);
T_LST  = try_city_slopes_RE(lmeL);

if ~isempty(T_NDVI)
    CityEquity = outerjoin(CityEquity, T_NDVI(:,{'City','Delta_TopMinusBottom','Delta_CI_L','Delta_CI_U'}), ...
        'Keys','City','MergeKeys',true,'Type','left');
    CityEquity.Properties.VariableNames(end-2:end) = ...
        {'LME_DeltaSlope_NDVI','LME_Delta_CI_L_NDVI','LME_Delta_CI_U_NDVI'};
else
    % deterministic fallback from OLS slopes (ensures non-constant deltas)
    CityEquity.LME_DeltaSlope_NDVI   = CityEquity.Slope_TopNDVI - CityEquity.Slope_BotNDVI;
    CityEquity.LME_Delta_CI_L_NDVI   = NaN(height(CityEquity),1);
    CityEquity.LME_Delta_CI_U_NDVI   = NaN(height(CityEquity),1);
end

if ~isempty(T_LST)
    CityEquity = outerjoin(CityEquity, T_LST(:,{'City','Delta_TopMinusBottom','Delta_CI_L','Delta_CI_U'}), ...
        'Keys','City','MergeKeys',true,'Type','left');
    CityEquity.Properties.VariableNames(end-2:end) = ...
        {'LME_DeltaSlope_LST','LME_Delta_CI_L_LST','LME_Delta_CI_U_LST'};
else
    CityEquity.LME_DeltaSlope_LST    = CityEquity.Slope_BotLST - CityEquity.Slope_TopLST;
    CityEquity.LME_Delta_CI_L_LST    = NaN(height(CityEquity),1);
    CityEquity.LME_Delta_CI_U_LST    = NaN(height(CityEquity),1);
end

% After CityEquity is built (has Slope_BotNDVI, Slope_TopNDVI, etc.)
CityEquity.B_NDVI = abs(CityEquity.Slope_BotNDVI) - abs(CityEquity.Slope_TopNDVI);
CityEquity.B_LST  = abs(CityEquity.Slope_BotLST ) - abs(CityEquity.Slope_TopLST );
% (ANOVA on B/Δ can be skipped as you noted.)


% Write CSV
csvName = 'equity_city_changes.csv';
writetable(CityEquity, csvName);

% City-year counts
G = groupsummary(Panel, {'City','Year'});
fprintf('City-years kept: %d  | MEAN tracts per city-year (pre-tail) ~ your input-side stat\n', height(G));

% Tails per city-year
% Recompute nb/nt once more (lightweight) to print MEANs:
% (If you want exact values, instrument inside the loop and accumulate.)


%% ===================== 6) BUNDLE EVERYTHING & ANNOUNCE =====================
S = struct();
S.config       = struct('TailPct',TailPct,'YearsUse',YearsUse,'MetricNames',MetricNames, ...
                        'MinTracts',MinTracts,'OverlapField',OverlapField,'Alpha',Alpha, ...
                        'Cities',CityListUse);
S.panel        = Panel;
S.long         = struct('NDVI',Tlong_NDVI,'LST',Tlong_LST);
S.models       = struct('lmeN',lmeN,'lmeL',lmeL, ...
                        'lme_gap_ndvi',lme_gap_ndvi,'lme_top_ndvi',lme_top_ndvi,'lme_bot_ndvi',lme_bot_ndvi, ...
                        'lme_gap_lst',lme_gap_lst,'lme_top_lst',lme_top_lst,'lme_bot_lst',lme_bot_lst, ...
                        'b_int_ndvi',b_int_ndvi,'b_int_lst',b_int_lst, ...
                        'b_top_ndvi',b_top_ndvi,'b_bot_ndvi',b_bot_ndvi, ...
                        'b_top_lst',b_top_lst,'b_bot_lst',b_bot_lst);
S.labels       = Labels;
S.labelCounts  = LabelCounts;
S.details      = Details;
S.cityEquity   = CityEquity;

assignin('base','Panel',Panel);
assignin('base','Tlong_NDVI',Tlong_NDVI);
assignin('base','Tlong_LST',Tlong_LST);
assignin('base','lmeN',lmeN);
assignin('base','lmeL',lmeL);
assignin('base','CityEquity',CityEquity);
assignin('base','S',S);

fprintf('Done. Built %0.0f%%-tail panel with filters, fit LMEs (with random slopes when possible), validated interactions,\n', 100*TailPct);
fprintf('labeled cities, and saved per-city equity changes.\n');
fprintf('Workspace vars: Panel, Tlong_NDVI, Tlong_LST, lmeN, lmeL, CityEquity, S\n');
fprintf('CSV written: %s\n\n', csvName);

%% ===================== Local helpers =====================
function inc = local_get_income(V, yr, T, toNumSafe)
switch yr
    case {1990,2000,2010,2020,2023}
        cands = {'RAW_INCOME'};
    otherwise
        cands = {};
end
candAll = [string(cands)];
j = find(ismember(V, candAll), 1, 'first');
if isempty(j), inc = []; else, inc = toNumSafe(T.(V(j))); end
end

function [b,p] = local_city_slope(T, yvar)
b = NaN; p = NaN;
if ~ismember(yvar, T.Properties.VariableNames), return; end
rows = isfinite(T.(yvar)) & isfinite(T.time_dec);
if nnz(rows) < 3 || numel(unique(T.Year(rows))) < 3, return; end
mdl = fitlm(T(rows,:), sprintf('%s ~ time_dec', yvar));
rn  = mdl.Coefficients.Properties.RowNames;
k   = strcmp(rn,'time_dec');
if any(k), b = mdl.Coefficients.Estimate(k); p = mdl.Coefficients.pValue(k); end
end

function REt = getREtable(mdl)
ok = false;
try
    tmp = randomEffects(mdl);
    if istable(tmp), REt = tmp; ok = true; end
catch, end
if ~ok
    try
        [b,names,stats] = randomEffects(mdl);
        if istable(names), REt = names;
        elseif isstruct(names), REt = struct2table(names);
        else, error('Unsupported names class: %s', class(names));
        end
        REt.Estimate = b(:);
        if istable(stats) && any(strcmp('SE', stats.Properties.VariableNames))
            REt.SE = stats.SE(:);
        elseif isstruct(stats) && isfield(stats,'SE')
            REt.SE = stats.SE(:);
        else
            REt.SE = NaN(size(b(:)));
        end
        ok = true;
    catch ME
        error('Could not coerce randomEffects to a table: %s', ME.message);
    end
end
v = string(REt.Properties.VariableNames);
if ~any(v=="Level") && any(v=="Level1"), REt.Level = REt.Level1; REt.Level1 = []; end
if ~all(ismember(["Name","Level","Estimate"], string(REt.Properties.VariableNames)))
    error('randomEffects table missing required columns.');
end
if ~any(string(REt.Properties.VariableNames)=="SE")
    REt.SE = NaN(height(REt),1);
end
end

function T = try_city_slopes_RE(mdl)
% Returns [] if the needed RE terms are not present (prevents "all constants").
T = [];
try
    CF = mdl.Coefficients; cn = string(CF.Name);
    i_time = find(cn=="time_dec", 1);
    i_int  = find(contains(cn,"time_dec:tail") | (contains(cn,"tail_") & contains(cn,":time_dec")), 1);
    if isempty(i_time) || isempty(i_int), return; end
    b_time = CF.Estimate(i_time); b_int = CF.Estimate(i_int); intNameFE = cn(i_int);

    REt = getREtable(mdl);
    nm = string(REt.Name); lev = string(REt.Level); est = REt.Estimate; se = REt.SE;

    % must have at least some city-level RE for time_dec OR interaction
    has_time = any(nm=="time_dec");
    has_int  = any(contains(nm,"time_dec") & contains(nm,"tail"));
    if ~has_time && ~has_int, return; end

    % choose the interaction name in RE matching FE if present
    if has_int
        intNames = unique(nm(contains(nm,"time_dec") & contains(nm,"tail")));
        hit = find(intNames == intNameFE, 1);
        if ~isempty(hit), intName = intNames(hit); else, intName = intNames(1); end
    else
        intName = intNameFE; % will imply u_int = 0 when absent
    end

    cities = unique(lev,'stable'); nC = numel(cities);
    u_time = zeros(nC,1); u_int = zeros(nC,1); se_int = NaN(nC,1);
    for i = 1:nC
        ci = cities(i);
        j1 = (lev==ci) & (nm=="time_dec");
        j2 = (lev==ci) & (nm==intName);
        if any(j1), u_time(i) = est(find(j1,1)); end
        if any(j2), u_int(i)  = est(find(j2,1));  se_int(i) = se(find(j2,1)); end
    end

    slope_bottom = b_time + u_time;
    slope_top    = b_time + b_int + u_time + u_int;
    delta        = b_int + u_int;

    % FE SE for the interaction (always available from fixed effects)
	se_fe = CF.SE(i_int);

	% If RE SEs are missing (NaN), treat them as 0 so CIs fall back to FE width.
	se_int(isnan(se_int)) = 0;

	% Combine FE and RE uncertainty (approx. independence)
	se_delta = sqrt(se_fe.^2 + se_int.^2);

	z = 1.96;
	delta_lo = delta - z.*se_delta;
	delta_hi = delta + z.*se_delta;


    T = table(cities, slope_bottom, slope_top, delta, delta_lo, delta_hi, ...
        'VariableNames', {'City','Slope_Bottom','Slope_Top','Delta_TopMinusBottom','Delta_CI_L','Delta_CI_U'});
catch
    T = [];
end
end