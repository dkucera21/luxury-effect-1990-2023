%% build_census_tables.m  —  HARMONIZED EDITION (with TRACT_ID)
% Rebuilds CENSUS_TABLES_rebuilt from harmonized CSVs:
%   E:\LUXURY_NATL_FINAL\TRACTS\RAW_CENSUS\NEW_MAT_TABLES\T_YYYY_CITY.csv
%
% Adds/keeps a year-specific tract ID column and creates a unified TRACT_ID.
% Year → ID column:
%   1990: JOIN_TRACT
%   2000: GEOID
%   2010: GEOID
%   2020: AFFGEOID
%   2023: GEOIDFQ
%
% Inputs are expected to already include your harmonized vars; if PCT_OVERLAP
% is missing but LANDSAT_AREA is present, it will be computed.
%
% We:
%   • read a minimal column set (skip Shape_*),
%   • coerce numeric, map sentinel -222222222 → NaN,
%   • (optionally) compute PCT_OVERLAP if missing,
%   • drop rows with RAW_INCOME==0/NaN OR HOUSING_DENSITY==0/NaN,
%   • build CityListMaster by coverage rules.

%% ---------------- User settings ----------------
csvFolder  = 'E:\LUXURY_NATL_FINAL\TRACTS\RAW_CENSUS\NEW_MAT_TABLES'; % harmonized CSVs
CoverageThreshold   = 50;     % require PCT_OVERLAP >= this for a tract to be "valid"
MinTracts           = 5;      % minimum valid tracts per city-year
YearsUse            = [1990 2000 2010 2020 2023];
StrictMode          = true;   % true: all requested years for a city must pass; false: any year passes
SkipMissingCoverage = true;   % if PCT_OVERLAP missing, pass (true) or fail (false)
ExcludeCities       = ["Anchorage","Fairbanks","HiloHI","Honolulu","Ponce","SanJuan"]; % optional

% === ADD: CPI constants & mapping (base = 2023) ===
CPI = struct('y2023',304.702, 'y2020',258.811, 'y2010',218.056, 'y1999',166.600, 'y1989',124.000);
CPI_TARGET = CPI.y2023;
cpiByDataYear = containers.Map( ...
    {'1990','2000','2010','2020','2023'}, ...
    {CPI.y1989, CPI.y1999, CPI.y2010, CPI.y2020, CPI.y2023} );


% Columns guaranteed/expected in valid city-years
keepVars = [ ...
  "PCT_OVERLAP","CENSUS_TRACT_AREA", ...
  "MEAN_NDVI","MEDIAN_NDVI","MEAN_LST","MEDIAN_LST", ...
  "HOUSING_DENSITY","POP_DENSITY", ...
  "PERCENT_OCCUPIED_HOUSING","PERCENT_OWNED_HOUSING", ...
  "PERCENT_HISPANIC","PERCENT_WHITE","PERCENT_BLACK","PERCENT_ASIAN", ...
  "PERCENT_BACHELOR","PERCENT_GRADUATE","PERCENT_HS", ...
  "RAW_INCOME","RAW_INCOME_CPI"];

% Year-specific ID variable mapping
idMap = containers.Map( ...
  {'1990','2000','2010','2020','2023'}, ...
  {'JOIN_TRACT','GEOID','GEOID','AFFGEOID','GEOIDFQ'});

% Reasonable fallbacks if a file has a different label (kept if present)
idFallbacks = ["JOIN_TRACT","GEOID","AFFGEOID","GEOIDFQ","GEOID10","TRACT","CTID"];

%% ---------------- basic setup ----------------
assert(isfolder(csvFolder), 'Input folder not found: %s', csvFolder);
files = dir(fullfile(csvFolder, 'T_*.csv'));
if isempty(files), error('No T_*.csv files found in %s', csvFolder); end

toDbl = @(x) double(str2double(regexprep(string(x),'[,\s%$]','')));
CENSUS_TABLES_rebuilt = struct();

fprintf('Reading %d harmonized CSV(s) from %s\n', numel(files), csvFolder);

%% ---------------- 1) Load with selected columns ----------------
nRead = 0;
for k = 1:numel(files)
    fp = fullfile(files(k).folder, files(k).name);
    [~, base] = fileparts(files(k).name);

    % Expect name: T_YYYY_CITY  (YYYY may have had 'n' previously; accept both)
    tok = regexp(base, '^T_(\d{4})n?_(.+?)$', 'tokens', 'once');
    if isempty(tok)
        tok = regexp(base, '^T_(\d{4})n?_(.+)$', 'tokens', 'once');
        if isempty(tok)
            fprintf('  - Skip (name parse): %s\n', files(k).name);
            continue;
        end
    end
    yr = str2double(tok{1});
    cy = string(tok{2});
    fieldName = sprintf('T_%d_%s', yr, cy);  % normalized field name

    % Which ID column should we import for this year?
    idWanted = "";
    key = string(yr);
    if isKey(idMap, char(key))
        idWanted = string(idMap(char(key)));
    end

    % Build import options
    opts = detectImportOptions(fp);
    allVars = string(opts.VariableNames);

    % Drop Shape_* if present
    kill = ismember(allVars, ["Shape_Length","Shape_Area"]);
    allVars = allVars(~kill);

    % Collect desired columns: keepVars + year-specific ID + fallbacks + LANDSAT_AREA (for PCT_OVERLAP)
    wanted = unique([keepVars, idWanted, idFallbacks, "LANDSAT_AREA", "TRACT_ID"], 'stable');
    sel = intersect(allVars, wanted, 'stable');
    if isempty(sel)
        % Read all, we’ll prune later (extremely unlikely)
        sel = allVars;
    end
    opts.SelectedVariableNames = cellstr(sel);

    % Read
    T = readtable(fp, opts);

    % --- Ensure the target ID column exists or pick a fallback
    V = string(T.Properties.VariableNames);
    idColName = "";
    if idWanted ~= "" && ismember(idWanted, V)
        idColName = idWanted;
    else
        fb = intersect(idFallbacks, V, 'stable');
        if ~isempty(fb), idColName = fb(1); end
    end

    % --- Coerce numerics & clean sentinels; keep ID as string
    for v = 1:width(T)
        vn = string(T.Properties.VariableNames{v});
        col = T.(v);
        if vn == idColName || vn == "TRACT_ID"
            % force to string and trim
            T.(v) = strtrim(string(col));
        else
            if isnumeric(col)
                col(col == -222222222) = NaN;
                T.(v) = double(col);
            else
                s = string(col);
                s = regexprep(s, "^['""]+|['""]+$", "");
                s = strtrim(s);
                x = toDbl(s);
                x(x == -222222222) = NaN;
                if any(isfinite(x))
                    T.(v) = x;
                else
                    T.(v) = s; % non-numeric text – keep as string
                end
            end
        end
    end

	% === ADD: create RAW_INCOME_CPI (2023$) ===
	V = string(T.Properties.VariableNames);
	if ismember("RAW_INCOME", V)
		raw = T.RAW_INCOME;
		if ~isnumeric(raw), raw = toDbl(raw); end
		raw(raw == -222222222) = NaN;  % honor sentinel
		key = char(string(yr));
		if isKey(cpiByDataYear, key)
			factor = CPI_TARGET / cpiByDataYear(key);
			T.RAW_INCOME_CPI = raw .* factor;   % 2023 dollars
		else
			T.RAW_INCOME_CPI = NaN(height(T),1); % unmapped year → NA column
		end
	else
		% RAW_INCOME missing → still create placeholder so schema is stable
		T.RAW_INCOME_CPI = NaN(height(T),1);
	end


    % --- Create unified TRACT_ID as string (if we found an ID column)
    if idColName ~= ""
        T.TRACT_ID = strtrim(string(T.(idColName)));
    else
        warning('No ID column found for %s (year %d). Creating empty TRACT_ID.', files(k).name, yr);
        T.TRACT_ID = strings(height(T),1);
    end

    % --- Compute PCT_OVERLAP if missing but LANDSAT_AREA & CENSUS_TRACT_AREA present
    V = string(T.Properties.VariableNames);
    if ~ismember("PCT_OVERLAP", V) && all(ismember(["LANDSAT_AREA","CENSUS_TRACT_AREA"], V))
        den = double(T.CENSUS_TRACT_AREA);
        num = double(T.LANDSAT_AREA);
        pct = NaN(height(T),1);
        ok = isfinite(num) & isfinite(den) & den>0;
        pct(ok) = 100 * (num(ok) ./ den(ok));
        T.PCT_OVERLAP = pct;
    end

    % --- Ensure all desired columns exist (create missing as NaN) — except ID, which is string
    for v = keepVars
        if ~ismember(v, string(T.Properties.VariableNames))
            T.(char(v)) = NaN(height(T),1);
        end
    end
    if ~ismember("TRACT_ID", string(T.Properties.VariableNames))
        T.TRACT_ID = strings(height(T),1);
    end
    % Keep the year-specific ID column too (create empty if missing so downstream joins don’t break)
    if idWanted ~= "" && ~ismember(idWanted, string(T.Properties.VariableNames))
        T.(char(idWanted)) = strings(height(T),1);
    end

	% Drop rows if ANY of the following are 0 or NaN:
	raw = T.RAW_INCOME;        if ~isnumeric(raw), raw = toDbl(raw); end
	pop = T.POP_DENSITY;       if ~isnumeric(pop), pop = toDbl(pop); end
	pct = T.PCT_OVERLAP;       if ~isnumeric(pct), pct = toDbl(pct); end

	bad = (~isfinite(raw) | raw==0) ...
	   |  (~isfinite(pop) | pop==0) ...
	   |  (~isfinite(pct));          % no Landsat overlap → drop

	if any(bad), T(bad,:) = []; end


    % --- Final column order: TRACT_ID, year-specific ID, keepVars (only those present)
    finalVars = ["TRACT_ID"];
    if idWanted ~= "" && ismember(idWanted, string(T.Properties.VariableNames))
        finalVars = [finalVars, idWanted];
    end
    have = intersect(keepVars, string(T.Properties.VariableNames), 'stable');
    finalVars = [finalVars, have];
    % Append any other imported ID fallbacks (if they exist and aren’t duplicates)
    otherIDs = setdiff(intersect(idFallbacks, string(T.Properties.VariableNames), 'stable'), finalVars);
    finalVars = [finalVars, otherIDs];
    % Keep order that exists in table:
    finalVars = intersect(finalVars, string(T.Properties.VariableNames), 'stable');
    T = T(:, cellstr(finalVars));

    % --- Store
    CENSUS_TABLES_rebuilt.(fieldName) = T;
    nRead = nRead + 1;
    fprintf('  • Loaded %-40s → %s (%d rows) | ID: %s → TRACT_ID\n', ...
        files(k).name, fieldName, height(T), ternStr(idColName~="", idColName, "none"));
end
fprintf('Loaded %d table(s).\n', nRead);

%% ---------------- 2) Build CityListMaster from coverage rules ----------------
YearPattern = '^T_(\d{4})_(.+)$';
fns = string(fieldnames(CENSUS_TABLES_rebuilt));

City_col = {}; Year_col = []; Pass_col = []; Reason_col = {};
toDblLocal = @(x) double(str2double(regexprep(string(x),'[,\s%$]','')));

for i = 1:numel(fns)
    tname = fns(i);
    tk = regexp(tname, YearPattern, 'tokens', 'once');
    if isempty(tk), continue; end
    yr   = str2double(tk{1});  if ~ismember(yr, YearsUse), continue; end
    city = string(tk{2});

    T = CENSUS_TABLES_rebuilt.(char(tname));
    if ~istable(T) || height(T)==0
        City_col{end+1} = city; Year_col(end+1) = yr; Pass_col(end+1) = false; Reason_col{end+1} = "empty"; %#ok<AGROW>
        continue
    end

    if ~ismember("PCT_OVERLAP", string(T.Properties.VariableNames))
        if SkipMissingCoverage
            City_col{end+1} = city; Year_col(end+1) = yr; Pass_col(end+1) = true;  Reason_col{end+1} = "no coverage col (ignored)"; %#ok<AGROW>
        else
            City_col{end+1} = city; Year_col(end+1) = yr; Pass_col(end+1) = false; Reason_col{end+1} = "no coverage col (fail)";    %#ok<AGROW>
        end
        continue
    end

    pct = T.PCT_OVERLAP; if ~isnumeric(pct), pct = toDblLocal(pct); end
    pass = nnz(isfinite(pct) & pct >= CoverageThreshold) >= MinTracts;
    City_col{end+1} = city; Year_col(end+1) = yr; Pass_col(end+1) = pass; Reason_col{end+1} = ternStr(pass,"ok","too few valid tracts"); %#ok<AGROW>
end

CityPassDiag = table();
CityPassDiag.City   = string(City_col(:));
CityPassDiag.Year   = double(Year_col(:));
CityPassDiag.Pass   = logical(Pass_col(:));
CityPassDiag.Reason = string(Reason_col(:));

G = findgroups(CityPassDiag.City);
if StrictMode
    passByCity = splitapply(@all, CityPassDiag.Pass, G);
else
    passByCity = splitapply(@any, CityPassDiag.Pass, G);
end
uCities = splitapply(@(c) c(1), CityPassDiag.City, G);

CityListMaster = uCities(passByCity);
CityListMaster = setdiff(CityListMaster, ExcludeCities, 'stable');

fprintf('\nCity inclusion rule: PCT_OVERLAP >= %.0f%% & MinTracts >= %d; %s mode.\n', ...
    CoverageThreshold, MinTracts, ternStr(StrictMode,'STRICT','LENIENT'));
fprintf('Years considered: %s\n', strjoin(string(YearsUse), ', '));
fprintf('Included cities: %d\n', numel(CityListMaster));

%% ---------------- 3) Push to base ----------------
assignin('base','CENSUS_TABLES_rebuilt', CENSUS_TABLES_rebuilt);
assignin('base','CityListMaster',        CityListMaster);
assignin('base','CityPassDiag',          CityPassDiag);

%% ---------------- helpers ----------------
function out = ternStr(cond, a, b)
    if cond, out = a; else, out = b; end
end
