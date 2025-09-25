/************************************************************
 * File: 05_terraclimate_precip_cumulative.js
 * Purpose: For each polygon in a user-supplied FeatureCollection,
 *          export a monthly time series table of cumulative
 *          precipitation from TerraClimate (PR), with columns:
 *          date, cumulative_00_months, …, cumulative_12_months.
 *
 * Dataset: IDAHO_EPSCOR/TERRACLIMATE (monthly), band 'pr' (mm)
 * Notes:
 *  - cumulative_00_months == the current month’s PR (mm)
 *  - cumulative_01_months == PR of current month + previous 1 month
 *  - …
 *  - cumulative_12_months == PR of current month + previous 12 months
 *  - Uses native TerraClimate resolution (~4638 m).
 ************************************************************/


// ====================== USER PARAMETERS ======================

// 1) Polygons layer: upload your cities.gpkg/geojson to GEE Assets and paste its Asset ID.
var POLYGONS_ASSET = 'users/<your-username>/cities311';   // <-- EDIT

// 2) Polygon name/ID field (used in output filenames).
var NAME_FIELD     = 'name';                               // <-- EDIT

// 3) Date window (inclusive of start, exclusive of end for filterDate end bound).
var START_DATE     = '1989-01-01';
var END_DATE       = '2024-12-31';

// 4) Months to include (set 1..12 to include all months).
var MONTH_START    = 1;
var MONTH_END      = 12;

// 5) Export settings
var OUTPUT_FOLDER  = 'TC_PR_CUM_TS';   // Google Drive folder for CSVs
var SCALE_METERS   = 4638.3;           // TerraClimate native grid
var MAX_PIXELS     = 1e13;

// 6) Optional geometry buffering for very small polygons (to stabilize means)
var MIN_AREA_M2    = 5e6;   // buffer if polygon area < 5 km²
var BUFFER_METERS  = 2000;  // 2 km buffer

// 7) Run all polygons or a single one by name?
var RUN_ALL        = true;
var POLY_TO_RUN    = 'Los_Angeles';     // used only if RUN_ALL=false


// ====================== LOAD INPUTS ==========================

var polys = ee.FeatureCollection(POLYGONS_ASSET);
var fc    = RUN_ALL ? polys : polys.filter(ee.Filter.eq(NAME_FIELD, POLY_TO_RUN));

var dateFilter  = ee.Filter.date(START_DATE, END_DATE);
var monthFilter = ee.Filter.calendarRange(MONTH_START, MONTH_END, 'month');

var tcPR = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
              .select('pr')
              .filter(dateFilter)
              .filter(monthFilter)
              .sort('system:time_start');

print('Polygons to process:', fc.size());
print('TerraClimate PR monthly images selected:', tcPR.size());


// =================== HELPERS (server-safe) ===================

// Buffer geometry if area < threshold
function bufferedGeom(geom){
  var gArea = geom.area();
  return ee.Algorithms.If(gArea.lte(MIN_AREA_M2), geom.buffer(BUFFER_METERS), geom);
}

// Add cumulative precipitation bands (00..12 months) to a given month image.
// cumulative_00_months == current month only (renamed from 'pr').
function addCumulativeBands(img, collection){
  var dStart = ee.Date(img.get('system:time_start'));
  var base   = img.select('pr').rename('cumulative_00_months');

  // Build cumulative bands for 1..12 months
  var cumImg = ee.Image(base);
  for (var m = 1; m <= 12; m++){
    var from = dStart.advance(-m, 'month');
    var to   = dStart.advance(1, 'month');  // inclusive current month
    var sumM = collection.filterDate(from, to).sum().rename('cumulative_' + (m < 10 ? ('0'+m) : m) + '_months');
    cumImg = cumImg.addBands(sumM);
  }
  // Keep a formatted date string YYYY-MM for table export
  return cumImg.set('date', dStart.format('YYYY-MM'));
}


// ====================== EXPORT PER POLYGON ===================

var list = fc.toList(fc.size());
var N    = list.size().getInfo();

for (var i = 0; i < N; i++) {
  var feat = ee.Feature(list.get(i));
  var name = ee.String(feat.get(NAME_FIELD));
  var safe = name.replace(' ', '_').replace('/', '_'); // filename-safe
  var geom = ee.FeatureCollection([feat]).geometry();
  var geomUse = ee.Geometry(ee.Algorithms.If(geom.area().lte(MIN_AREA_M2), geom.buffer(BUFFER_METERS), geom));

  // Build a per-month image with cumulative bands, then reduceRegion to a row.
  var perMonthFeatures = tcPR.map(function(img){
    var cum = addCumulativeBands(img, tcPR).clip(geomUse);
    var dict = cum.reduceRegion({
      reducer:   ee.Reducer.mean(),
      geometry:  geomUse,
      scale:     SCALE_METERS,
      maxPixels: MAX_PIXELS,
      bestEffort: true
    });
    // Ensure we carry the date string through
    var dateStr = ee.String(cum.get('date'));
    return ee.Feature(null, dict.set('date', dateStr));
  })
  // Keep only months that yielded data (non-null for current month band)
  .filter(ee.Filter.notNull(['cumulative_00_months']));

  // Column order for export
  var selectors = [
    'date',
    'cumulative_00_months',
    'cumulative_01_months',
    'cumulative_02_months',
    'cumulative_03_months',
    'cumulative_04_months',
    'cumulative_05_months',
    'cumulative_06_months',
    'cumulative_07_months',
    'cumulative_08_months',
    'cumulative_09_months',
    'cumulative_10_months',
    'cumulative_11_months',
    'cumulative_12_months'
  ];

  var desc = ee.String('TC_PR_CUM_TS_').cat(safe).getInfo();

  Export.table.toDrive({
    collection:  perMonthFeatures,
    description: desc,
    fileNamePrefix: desc,
    folder:      OUTPUT_FOLDER,
    fileFormat:  'CSV',
    selectors:   selectors
  });

  print('Queued export:', desc);
}


// ============================ NOTES ==========================
// • Units are in mm (TerraClimate 'pr').
// • The date column is YYYY-MM for each monthly row in the requested window.
// • If you want seasonal/annual cumulative PR, aggregate downstream by grouping
//   the rows (e.g., with Pandas) rather than recomputing in GEE.
// • For large batches, run subsets (filter fc) to avoid task limits.
