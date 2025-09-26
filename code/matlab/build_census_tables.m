%% build_census_tables.m
% Rebuild CENSUS_TABLES_new from raw CSVs and apply standardized renames/fields
% - Reads all *.csv from csvFolder
% - Adds PCT_OVERLAP = 100 * LANDSAT_AREA / AREAn
% - Normalizes table names "T_2020n_*" / "T_2023n_*" → "T_2020_*" / "T_2023_*"
% - 2000: renames + POP_DENSITY & HOUSING_DENSITY
% - 1990: HOUSING_DENSITY from H0010001 (total housing units)
% - 2010: ACS-style renames + POP/HOUSING densities from canonical vars
% - 2020/2023: ACS-style renames + POP/HOUSING densities (fallback-aware)
% - 2023-only extra renames
% - Zero→NaN for densities in all tables
%
% Output:
%   - Struct `CENSUS_TABLES_new` in workspace
%   - MAT file saved to `outMatPath`
%
% Notes:
%   - No `clear`, `eval`, or `assignin` calls
%   - Robust to case/commas in numeric fields
%   - Safe divisions guard against zero/NaN areas

%% ---------------- User settings ----------------
csvFolder  = 'E:\LUXURY_NATL_FINAL\TRACTS\RAW_CENSUS\NEW_MAT_TABLES';  % <-- update to your extracted CSV folder
outMatPath = fullfile(pwd, 'outputs', 'CENSUS_TABLES_new.mat');        % change if you prefer

if ~isfolder(csvFolder)
    error('Input folder not found: %s', csvFolder);
end
outDir = fileparts(outMatPath);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end

%% ---------------- 1) Load all CSVs into a struct ----------------
files = dir(fullfile(csvFolder, '*.csv'));
if isempty(files)
    error('No CSV files found in: %s', csvFolder);
end

CENSUS_TABLES_new = struct();
fprintf('Reading %d CSV file(s) from %s\n', numel(files), csvFolder);

for k = 1:numel(files)
    filePath = fullfile(csvFolder, files(k).name);
    [~, name] = fileparts(files(k).name);
    fieldName = matlab.lang.makeValidName(name);

    % Robust import (let MATLAB infer types, keep text as strings if needed)
    opts = detectImportOptions(filePath, 'NumHeaderLines', 0);
    T    = readtable(filePath, opts);

    CENSUS_TABLES_new.(fieldName) = T;
    fprintf('  • Loaded %s → CENSUS_TABLES_new.%s\n', files(k).name, fieldName);
end

%% ---------------- 2) Add PCT_OVERLAP in each table ----------------
fn = fieldnames(CENSUS_TABLES_new);
for k = 1:numel(fn)
    tblName = fn{k};
    T = CENSUS_TABLES_new.(tblName);
    V = T.Properties.VariableNames;

    vAreaN   = findVarCase(V, 'AREAn');
    vLandsat = findVarCase(V, 'LANDSAT_AREA');

    if isempty(vAreaN) || isempty(vLandsat)
        fprintf('  - %s: missing AREAn or LANDSAT_AREA — skipped PCT_OVERLAP\n', tblName);
    else
        areaVals = toDbl(T.(vAreaN));
        lsVals   = toDbl(T.(vLandsat));
        pct      = 100 .* (lsVals ./ areaVals);
        bad      = ~isfinite(areaVals) | areaVals <= 0;
        pct(bad) = NaN;
        T.PCT_OVERLAP = pct;
        fprintf('  • %s: PCT_OVERLAP added\n', tblName);
    end

    CENSUS_TABLES_new.(tblName) = T;
end

%% ---------------- 3) Normalize table names: 2020n/2023n → 2020/2023 ----------------
fn = fieldnames(CENSUS_TABLES_new);
for k = 1:numel(fn)
    oldName = fn{k};
    newName = regexprep(oldName, '^T_(2020|2023)n_', 'T_$1_');
    if ~strcmp(oldName, newName)
        CENSUS_TABLES_new.(newName) = CENSUS_TABLES_new.(oldName);
        CENSUS_TABLES_new = rmfield(CENSUS_TABLES_new, oldName);
        fprintf('  • Renamed table: %s → %s\n', oldName, newName);
    end
end

%% ---------------- 4) YEAR 2000: renames + densities ----------------
tblNames = fieldnames(CENSUS_TABLES_new);
for i = 1:numel(tblNames)
    tname = tblNames{i};
    if ~startsWith(tname,'T_2000_'), continue; end
    T = CENSUS_TABLES_new.(tname);
    V = T.Properties.VariableNames;

    T = renameVar(T, 'PERCENT_OCC_HOUSING', 'PERCENT_OCCUPIED_HOUSING');
    T = renameVar(T, 'PERCENT_GRAD',        'PERCENT_GRADUATE');
    T = renameVar(T, 'PERCENT_BACH',        'PERCENT_BACHELOR');
    T = renameVar(T, 'PERCENT_OWN',         'PERCENT_OWNED_HOUSING');

    aName   = findVarCase(V, 'AREAn');
    occName = findVarCase(V, 'OCCUPIED_HOUSING');
    rentName= findVarCase(V, 'RENTER_HOUSING');
    popName = findVarCase(V, 'TOTAL_POP');

    area_km2 = [];
    if ~isempty(aName)
        area_km2 = toDbl(T.(aName)) ./ 1e6;
    end

    if ~isempty(area_km2) && ~isempty(occName) && ~isempty(rentName)
        units2000 = toDbl(T.(occName)) + toDbl(T.(rentName));
        T.HOUSING_DENSITY = safeDiv(units2000, area_km2);
    end
    if ~isempty(area_km2) && ~isempty(popName)
        T.POP_DENSITY = safeDiv(toDbl(T.(popName)), area_km2);
    end

    CENSUS_TABLES_new.(tname) = T;
    fprintf('  • Updated 2000 table: %s\n', tname);
end

%% ---------------- 5) YEAR 1990: housing density from H0010001 ----------------
tblNames = fieldnames(CENSUS_TABLES_new);
for i = 1:numel(tblNames)
    tname = tblNames{i};
    if ~startsWith(tname,'T_1990_'), continue; end
    T = CENSUS_TABLES_new.(tname);
    V = T.Properties.VariableNames;

    aName  = findVarCase(V, 'AREAn');     % m^2
    totHou = findVarCase(V, 'H0010001');  % total housing units

    if ~isempty(aName) && ~isempty(totHou)
        area_km2 = toDbl(T.(aName));
        T.HOUSING_DENSITY = safeDiv(toDbl(T.(totHou)), area_km2);
        CENSUS_TABLES_new.(tname) = T;
        fprintf('  • 1990 HOUSING_DENSITY (total units): %s\n', tname);
    else
        fprintf('  - 1990 skip (AREAn/H0010001 missing): %s\n', tname);
    end
end

%% ---------------- 6) YEAR 2010: ACS-style renames + densities ----------------
tblNames = fieldnames(CENSUS_TABLES_new);
for i = 1:numel(tblNames)
    tname = tblNames{i};
    if ~startsWith(tname,'T_2010_'), continue; end
    T = CENSUS_TABLES_new.(tname);
    V = T.Properties.VariableNames;

    T = renameVar(T, 'Percent_HispanicOrigin', 'PERCENT_HISPANIC');
    T = renameVar(T, 'Percent_White',          'PERCENT_WHITE');
    T = renameVar(T, 'Percent_Black',          'PERCENT_BLACK');
    T = renameVar(T, 'Percent_Indian',         'PERCENT_NATIVE');
    T = renameVar(T, 'Percent_Asian',          'PERCENT_ASIAN');
    T = renameVar(T, 'Percent_Other',          'PERCENT_OTHER_RACE');
    T = renameVar(T, 'Percent_Ownerocc',       'PERCENT_OWNED_HOUSING');
    T = renameVar(T, 'Percent_Occupied',       'PERCENT_OCCUPIED_HOUSING');
    T = renameVar(T, 'PERCENT_BACH',           'PERCENT_BACHELOR');
    T = renameVar(T, 'UNDER_5_percent',        'PERCENT_UNDER_5');

    V = T.Properties.VariableNames; % refresh
    aName = findVarCase(V, 'AREAn');
    if isempty(aName)
        CENSUS_TABLES_new.(tname) = T;
        fprintf('  - 2010 no AREAn; densities skipped: %s\n', tname);
        continue;
    end
    area_km2 = toDbl(T.(aName)) ./ 1e6;

    % Canonical ACS totals (prefer these; otherwise fallback)
    popName = pickFirstField(V, {'B02001_001E','B01003_001E','B01001_001E','P001001','P0010001','TOTAL_POP','TOTALPOP','POP','POPULATION'});
    houName = pickFirstField(V, {'B25003_001E','H001001','H0010001','TOTAL_HOUSING','TOTALHOUSING','HOUSING_UNITS','HOUSINGUNITS'});

    if ~isempty(popName)
        T.POP_DENSITY = safeDiv(toDbl(T.(popName)), area_km2);
    end
    if ~isempty(houName)
        T.HOUSING_DENSITY = safeDiv(toDbl(T.(houName)), area_km2);
    end

    CENSUS_TABLES_new.(tname) = T;
    fprintf('  • Updated 2010 table: %s\n', tname);
end

%% ---------------- 7) YEARS 2020/2023: ACS-style renames + densities ----------------
tblNames = fieldnames(CENSUS_TABLES_new);
norm = @(s) lower(regexprep(string(s),'[^A-Za-z0-9]+',''));

renamePairsACS = {
  'Percent_HispanicOrigin', 'PERCENT_HISPANIC';
  'Percent_White',           'PERCENT_WHITE';
  'Percent_Black',           'PERCENT_BLACK';
  'Percent_Indian',          'PERCENT_NATIVE';
  'Percent_Asian',           'PERCENT_ASIAN';
  'Percent_Other',           'PERCENT_OTHER_RACE';
  'Percent_Ownerocc',        'PERCENT_OWNED_HOUSING';
  'Percent_Occupied',        'PERCENT_OCCUPIED_HOUSING';
  'PERCENT_BACH',            'PERCENT_BACHELOR';
  'UNDER_5_percent',         'PERCENT_UNDER_5';
};

popACS_candidates = { ...
  'B01003_001E','B01001_001E','B02001_001E', ...
  'P001001','P0010001', ...
  'TOTAL_POP','TOTALPOP','POP','POPULATION'};

housingACS_candidates = { ...
  'B25003_001E', ...
  'H001001','H0010001', ...
  'TOTAL_HOUSING','TOTALHOUSING','HOUSING_UNITS','HOUSINGUNITS'};

yearPrefixes = {'T_2020_','T_2023_'};
nDone = 0;

for i = 1:numel(tblNames)
    tname = tblNames{i};
    if ~any(startsWith(tname, yearPrefixes)), continue; end
    T = CENSUS_TABLES_new.(tname);
    V = T.Properties.VariableNames;
    Vn = norm(V);

    % Renames
    didRename = false;
    for r = 1:size(renamePairsACS,1)
        src = renamePairsACS{r,1};
        dst = renamePairsACS{r,2};
        j = find(Vn == norm(src), 1, 'first');
        if ~isempty(j) && ~strcmp(V{j}, dst)
            V{j} = dst;
            didRename = true;
        end
    end
    if didRename
        T.Properties.VariableNames = V;
        V = T.Properties.VariableNames; Vn = norm(V);
    end

    % Densities
    jArea = find(Vn == norm('AREAn'), 1, 'first');
    if isempty(jArea)
        CENSUS_TABLES_new.(tname) = T;
        fprintf('  - %s: AREAn not found; densities skipped\n', tname);
        continue;
    end
    area_km2 = toDbl(T.(V{jArea})) ./ 1e6;

    popName = pickFirstField(V, popACS_candidates);
    if ~isempty(popName)
        T.POP_DENSITY = safeDiv(toDbl(T.(popName)), area_km2);
    end

    housName = pickFirstField(V, housingACS_candidates);
    if ~isempty(housName)
        T.HOUSING_DENSITY = safeDiv(toDbl(T.(housName)), area_km2);
    end

    CENSUS_TABLES_new.(tname) = T;
    nDone = nDone + 1;
end
fprintf('  • Processed %d table(s) for 2020/2023.\n', nDone);

%% ---------------- 8) 2023-only extra renames ----------------
renamePairs2023 = {
  'Percent_Bachelor',   'PERENT_BACHELOR'; % (keep as-is per your note)
  'Percent_GradDegree', 'PERCENT_GRAD';
  'Percent_HighSchool', 'PERCENT_HS';
  'PercentUnder5',      'PERCENT_UNDER_5';
};
tblNames = fieldnames(CENSUS_TABLES_new);
for i = 1:numel(tblNames)
    tname = tblNames{i};
    if ~startsWith(tname,'T_2023_'), continue; end
    T = CENSUS_TABLES_new.(tname);
    V = T.Properties.VariableNames;
    Vn = norm(V);

    didRename = false;
    for r = 1:size(renamePairs2023,1)
        src = renamePairs2023{r,1};
        dst = renamePairs2023{r,2};
        j = find(Vn == norm(src), 1, 'first');
        if ~isempty(j) && ~strcmp(V{j}, dst)
            V{j} = dst;
            didRename = true;
        end
    end
    if didRename
        T.Properties.VariableNames = V;
        CENSUS_TABLES_new.(tname) = T;
        fprintf('  • 2023 renames applied: %s\n', tname);
    end
end

%% ---------------- 9) Zero→NaN for densities (all tables) ----------------
tblNames = fieldnames(CENSUS_TABLES_new);
touchedPop = 0; touchedHou = 0;
for i = 1:numel(tblNames)
    tname = tblNames{i};
    T = CENSUS_TABLES_new.(tname);
    V = T.Properties.VariableNames;
    Vn = norm(V);

    jPop = find(Vn == norm('POP_DENSITY'), 1, 'first');
    if ~isempty(jPop)
        col = V{jPop};
        x = toDbl(T.(col));
        wasZero = isfinite(x) & (x == 0);
        if any(wasZero)
            x(wasZero) = NaN;
            T.(col) = x;
            touchedPop = touchedPop + 1;
        end
    end

    jHou = find(Vn == norm('HOUSING_DENSITY'), 1, 'first');
    if ~isempty(jHou)
        col = V{jHou};
        x = toDbl(T.(col));
        wasZero = isfinite(x) & (x == 0);
        if any(wasZero)
            x(wasZero) = NaN;
            T.(col) = x;
            touchedHou = touchedHou + 1;
        end
    end

    CENSUS_TABLES_new.(tname) = T;
end
fprintf('  • Zero→NaN applied in POP_DENSITY for %d table(s)\n', touchedPop);
fprintf('  • Zero→NaN applied in HOUSING_DENSITY for %d table(s)\n', touchedHou);

%% Clean tokens in CENSUS_TABLES_new: '', '-', '**', '***' -> missing
% PLUS: set any numeric value == -222222222 to NaN
% Robust to quoted tokens like "'-'" or '"**"'.
% Leaves logical variables untouched.

if ~exist('CENSUS_TABLES_new','var') || ~isstruct(CENSUS_TABLES_new)
    error('CENSUS_TABLES_new not found or not a struct.');
end

TOK_STR = ["", "-", "**", "***", "-222222222"];  % treat these strings as missing
SENTINEL = -222222222;                           % numeric sentinel to convert to NaN

tblNames = fieldnames(CENSUS_TABLES_new);
grandTotal = 0;

for ti = 1:numel(tblNames)
    tname = tblNames{ti};
    T = CENSUS_TABLES_new.(tname);
    if ~istable(T), continue; end

    n = height(T);
    changed = 0;

    for v = 1:width(T)
        colName = T.Properties.VariableNames{v};
        col = T.(v);

        % --- Fast path: numeric columns → replace sentinel with NaN
        if isnumeric(col)
            mask = (col == SENTINEL);
            nFix = nnz(mask);
            if nFix > 0
                col(mask) = NaN;
                T.(v) = col;
                changed = changed + nFix;
            end
            continue  % numeric handled; move to next variable
        end

        % Skip logical
        if islogical(col)
            continue
        end

        % --- Coerce to a string array (same NUMEL as original), robust across types ---
        try
            if isstring(col)
                s = col;
            elseif iscell(col)
                % cell array -> string array of same size
                s = strings(size(col));
                for k = 1:numel(col)
                    x = col{k};
                    if isstring(x)
                        s(k) = strjoin(x(:).', "");
                    elseif ischar(x)
                        s(k) = string(x);
                    elseif isempty(x)
                        s(k) = "";
                    else
                        s(k) = string(x);
                    end
                end
            elseif ischar(col)
                % char matrix (n×m) -> cellstr (n×1) -> string (n×1)
                s = string(cellstr(col));
            elseif iscategorical(col)
                s = string(col);
            else
                % last resort
                s = string(col);
            end
        catch
            fprintf('Skipping %s.%s: could not convert to string safely.\n', tname, colName);
            continue
        end

        % Ensure length matches rows
        if numel(s) ~= n
            if size(s,1) == n
                s = s(:,1);
            elseif size(s,2) == n && size(s,1) == 1
                s = s.';
            else
                fprintf('Skipping %s.%s: size %s not compatible with %d rows.\n', ...
                        tname, colName, mat2str(size(s)), n);
                continue
            end
        end

        % Strip surrounding quotes and whitespace
        s = regexprep(s, "^['""]+", "");
        s = regexprep(s, "['""]+$", "");
        s = strtrim(s);

        % Count & standardize missing tokens (includes "-222222222" as text)
        willMiss = ismember(s, TOK_STR);
        nMiss = nnz(willMiss);
        s = standardizeMissing(s, TOK_STR);

        % Write back
        try
            T.(v) = s;
            changed = changed + nMiss;
        catch ME
            fprintf('Skipping write-back %s.%s: %s\n', tname, colName, ME.message);
        end
    end

    CENSUS_TABLES_new.(tname) = T;
    if changed > 0
        fprintf('Table %-30s : set %6d entries to missing/NaN.\n', tname, changed);
        grandTotal = grandTotal + changed;
    end
end

fprintf('Done. Total replacements across all tables: %d\n', grandTotal);

%% ── Derive 1990 variables on all T_1990_* tables in CENSUS_TABLES_new ──
% Creates: POP_DENSITY, HOUSING_DENSITY, PERCENT_UNDER_5,
%          PERCENT_OCCUPIED_HOUSING, PERCENT_OWNED_HOUSING,
%          PERCENT_{HISPANIC,WHITE,BLACK,NATIVE,ASIAN,OTHER_RACE},
%          PERCENT_{BACHELOR,GRADUATE,HS}

tblNames = fieldnames(CENSUS_TABLES_new);
is1990   = startsWith(tblNames, 'T_1990_', 'IgnoreCase', true);

missingLog = struct();    % collect missing-required fields per table
nDone = 0;

for ii = 1:numel(tblNames)
    if ~is1990(ii), continue; end
    tname = tblNames{ii};
    T = CENSUS_TABLES_new.(tname);
    n = height(T);

    % --- helpers for this table ---
    V   = T.Properties.VariableNames;
    get = @(nm) getNumCol(T, V, nm);     % numeric column (NaN if missing / unparsable)

    % --- base inputs ---
    A     = get('AREAn');          % area
    POP   = get('P0010001');       % total population
    Htot  = get('H0010001');       % total housing units
    U5    = get('P0130001') + get('P0130002') + get('P0130003') + get('P0130004');
    Hocc  = get('H0040001');       % occupied housing units
    Hown  = get('H0080001');       % owner-occupied units

    Hisp  = get('P0100001');
    White = get('P0120001');
    Black = get('P0120002');
    Native= get('P0120003');
    Asian = get('P0120004');
    Other = get('P0120005');

    % Education 25+ denominator
    Eden = get('P0570001') + get('P0570002') + get('P0570003') + ...
           get('P0570004') + get('P0570005') + get('P0570006') + get('P0570007');
    EBach = get('P0570006');
    EGrad = get('P0570007');
    EHS   = get('P0570003');

    % --- derived metrics ---
    T.POP_DENSITY             = safeDiv(POP,  A);                 % per AREA unit
    T.HOUSING_DENSITY         = safeDiv(Htot, A);                 % per AREA unit
    T.PERCENT_UNDER_5         = cap01( 100 * safeDiv(U5,   POP) );
    T.PERCENT_OCCUPIED_HOUSING= cap01( 100 * safeDiv(Hocc, Htot) );
    T.PERCENT_OWNED_HOUSING   = cap01( 100 * safeDiv(Hown, Hocc) );

    T.PERCENT_HISPANIC        = cap01( 100 * safeDiv(Hisp,  POP) );
    T.PERCENT_WHITE           = cap01( 100 * safeDiv(White, POP) );
    T.PERCENT_BLACK           = cap01( 100 * safeDiv(Black, POP) );
    T.PERCENT_NATIVE          = cap01( 100 * safeDiv(Native,POP) );
    T.PERCENT_ASIAN           = cap01( 100 * safeDiv(Asian, POP) );
    T.PERCENT_OTHER_RACE      = cap01( 100 * safeDiv(Other, POP) );

    T.PERCENT_BACHELOR        = cap01( 100 * safeDiv(EBach, Eden) );
    T.PERCENT_GRADUATE        = cap01( 100 * safeDiv(EGrad, Eden) );
    T.PERCENT_HS              = cap01( 100 * safeDiv(EHS,   Eden) );

    % Log which required inputs were missing (if any)
    required = {'AREAn','P0010001','H0010001','P0130001','P0130002','P0130003','P0130004', ...
                'H0040001','H0080001','P0100001','P0120001','P0120002','P0120003','P0120004','P0120005', ...
                'P0570001','P0570002','P0570003','P0570004','P0570005','P0570006','P0570007'};
    miss = required(~isMemberCI(required, V));
    if ~isempty(miss)
        missingLog.(tname) = miss;
    end

    % write back
    CENSUS_TABLES_new.(tname) = T;
    nDone = nDone + 1;
end

fprintf('✓ Updated %d T_1990_* tables with derived variables.\n', nDone);

% Report any missing fields (so you can patch upstream if needed)
if ~isempty(fieldnames(missingLog))
    fprintf('\nTables with missing source fields (skipped in affected formulas):\n');
    tn = fieldnames(missingLog);
    for k = 1:numel(tn)
        fprintf('  %s: %s\n', tn{k}, strjoin(missingLog.(tn{k}), ', '));
    end
else
    fprintf('All required source fields were found in every processed table.\n');
end

%% Convert 1990 densities from per m^2 to per km^2
tblNames = fieldnames(CENSUS_TABLES_new);
targets  = {'POP_DENSITY','HOUSING_DENSITY'};

nTablesChanged = 0;
logRows = [];

for i = 1:numel(tblNames)
    tname = tblNames{i};
    if ~startsWith(tname, 'T_1990_', 'IgnoreCase', true), continue; end

    T = CENSUS_TABLES_new.(tname);
    V = T.Properties.VariableNames;

    didAny = false;
    for t = 1:numel(targets)
        want = targets{t};
        j = find(strcmpi(V, want), 1, 'first');
        if isempty(j), continue; end
        vn = V{j};

        x = T.(vn);
        % coerce to numeric if needed
        if ~isnumeric(x)
            x = str2double(regexprep(string(x), '[,\s%$]', ''));
        else
            x = double(x);
        end

        beforeMed = median(x, 'omitnan');

        % convert: people (or units) per m^2 -> per km^2
        x = x * 1e6;

        afterMed = median(x, 'omitnan');
        T.(vn) = x;
        didAny = true;

        logRows = [logRows; {tname, vn, beforeMed, afterMed}]; %#ok<AGROW>
    end

    if didAny
        CENSUS_TABLES_new.(tname) = T;
        nTablesChanged = nTablesChanged + 1;
    end
end

fprintf('✓ Converted densities to per km^2 in %d T_1990_* tables.\n', nTablesChanged);
if ~isempty(logRows)
    fprintf('Examples (median before → after):\n');
    % show up to 10 lines
    showN = min(10, size(logRows,1));
    for k = 1:showN
        fprintf('  %s :: %s   %.6g → %.6g\n', logRows{k,1}, logRows{k,2}, logRows{k,3}, logRows{k,4});
    end
end

fn = fieldnames(CENSUS_TABLES_new);
for i = 1:numel(fn)
    m = regexp(fn{i}, '^T_(\d{4})n_(.+)$', 'tokens', 'once');
    if isempty(m), continue; end
    newName = sprintf('T_%s_%s', m{1}, m{2});
    CENSUS_TABLES_new.(newName) = CENSUS_TABLES_new.(fn{i});
    CENSUS_TABLES_new = rmfield(CENSUS_TABLES_new, fn{i});
end

% Rename income fields to MED_INCOME in CENSUS_TABLES_new
assert(exist('CENSUS_TABLES_new','var')==1 && isstruct(CENSUS_TABLES_new), ...
    'CENSUS_TABLES_new not found or not a struct.');

tbls = fieldnames(CENSUS_TABLES_new);
nDone = 0;

for i = 1:numel(tbls)
    tname = tbls{i};
    tok = regexp(tname, '^T_(\d{4})n?_', 'tokens', 'once');
    if isempty(tok), continue; end
    yr = str2double(tok{1});

    % Pick the source column name by year
    switch yr
        case 2010
            src = 'S1903_C02_001E';
        case {2020, 2023}
            src = 'S1903_C03_001E';
        otherwise
            continue;
    end

    T = CENSUS_TABLES_new.(tname);
    if ~istable(T), continue; end

    vnames = T.Properties.VariableNames;
    j = find(strcmpi(vnames, src), 1, 'first');
    if isempty(j), continue; end

    % Rename: create/overwrite MED_INCOME, then remove the source column
    T.MED_INCOME = T.(vnames{j});
    T.(vnames{j}) = [];
    CENSUS_TABLES_new.(tname) = T;

    nDone = nDone + 1;
    fprintf('Renamed %s → MED_INCOME in %s\n', src, tname);
end

fprintf('Done. Renamed MED_INCOME in %d table(s).\n', nDone);

% Add income10k and incRank to CENSUS_TABLES_rebuilt
threshold  = 50;
pctNames   = {'PCT_OVERLAP','PCT_COVER','PCT_COV','PercentAreaNDVI'};
numify     = @(x) str2double(regexprep(string(x),'[,\$%]',''));

rawIncByYear = containers.Map( ...
  {'1990','2000','2010','2020','2023'}, ...
  { {'P080A001','MED_INCOME'}, ...
    {'DP3_C112','MED_INCOME'}, ...
    {'S1903_C02_001E','MED_INCOME'}, ...
    {'S1903_C03_001E','S1903_003_001E','MED_INCOME'}, ...
    {'S1903_C03_001E','S1903_003_001E','MED_INCOME'} });

tblNames = fieldnames(CENSUS_TABLES_rebuilt);

for i = 1:numel(tblNames)
    tname = tblNames{i};
    T = CENSUS_TABLES_rebuilt.(tname); if ~istable(T), continue, end
    tok = regexp(tname,'^T_(\d{4})_','tokens','once'); if isempty(tok), continue, end
    yr = tok{1};  V = T.Properties.VariableNames;

    % overlap column
    vPct=''; for k=1:numel(pctNames), j=find(strcmpi(V,pctNames{k}),1); if ~isempty(j), vPct=V{j}; break, end, end
    if isempty(vPct), continue, end
    pct = numify(T.(vPct));

    % income field (prefer canonical raw)
    vInc=''; C = rawIncByYear(yr);
    for k=1:numel(C), j=find(strcmpi(V,C{k}),1); if ~isempty(j), vInc=V{j}; break, end, end
    if isempty(vInc), continue, end
    inc = numify(T.(vInc));
    inc(inc==0 | inc==-222222222) = NaN;

    valid = pct>=threshold & isfinite(inc);

    % income per $10k
    income10k = NaN(height(T),1);
    income10k(valid) = inc(valid)/10000;

    % incRank in [0,1]
    incRank = NaN(height(T),1);
    if nnz(valid) >= 5
        [~,ord] = sort(inc(valid));
        r = NaN(nnz(valid),1);
        r(ord) = ((1:nnz(valid))' - 0.5) / nnz(valid);  % fractional mid-rank
        incRank(valid) = r;
    end

    T.income10k = income10k;
    T.incRank   = incRank;
    CENSUS_TABLES_rebuilt.(tname) = T;
end

fprintf('Added income10k and incRank where possible.\n');



clearvars -except CENSUS_TABLES_new

%% ===================== helpers (local functions) =====================
function nm = findVarCase(V, target)
% Case-insensitive exact match; returns '' if not found
    idx = find(strcmpi(V, target), 1, 'first');
    if isempty(idx), nm = ''; else, nm = V{idx}; end
end

function T = renameVar(T, fromName, toName)
% Case-insensitive rename; no-op if source missing
    V = T.Properties.VariableNames;
    idx = find(strcmpi(V, fromName), 1, 'first');
    if ~isempty(idx)
        V{idx} = toName;
        T.Properties.VariableNames = V;
    end
end

function name = pickFirstField(V, candidates)
% Return the first V member that matches any of the candidates (case/format tolerant)
    name = '';
    if isempty(V), return; end
    Vn = lower(regexprep(string(V),'[^A-Za-z0-9]+',''));
    for k = 1:numel(candidates)
        want = lower(regexprep(string(candidates{k}),'[^A-Za-z0-9]+',''));
        j = find(Vn == want, 1, 'first');
        if ~isempty(j)
            name = V{j};
            return;
        end
    end
end

function x = toDbl(x)
% Robust numeric conversion for vectors/columns; strips commas/spaces
    if istable(x), error('toDbl received a table, not a column.'); end
    if isnumeric(x), x = double(x); return; end
    x = str2double(regexprep(string(x),'[,\s]',''));
end

function z = safeDiv(num, den)
% Elementwise safe division with guards for zeros/NaNs
    num = toDbl(num); den = toDbl(den);
    z = nan(size(num));
    ok = isfinite(num) & isfinite(den) & den > 0;
    z(ok) = num(ok) ./ den(ok);
end

function out = getNumCol(T, V, target)
% Return numeric column from table T by case-insensitive name (or NaNs if missing).
    idx = find(strcmpi(V, target), 1, 'first');
    if isempty(idx)
        out = nan(height(T),1);
        return
    end
    x = T.(V{idx});
    if isnumeric(x)
        out = double(x);
    else
        out = str2double(regexprep(string(x), '[,\s%$]', ''));
    end
end

function y = cap01(y)
% Clamp percentages to [0,100] where rounding/quirks push slightly outside.
    y(y < 0 & isfinite(y))   = 0;
    y(y > 100 & isfinite(y)) = 100;
end

function tf = isMemberCI(need, V)
% Case-insensitive membership for variable names.
    VN = upper(string(V));
    tf = ismember(upper(string(need)), VN);
end
