/************************************************************
 * File: 06_spei_monthly_timeseries.js
 * Purpose: Export monthly SPEI time series for each polygon
 *          in a user-supplied FeatureCollection (e.g., cities).
 *
 * Dataset: CSIC/SPEI/2_10  (multi-scale SPEI, monthly)
 * Output:  One CSV per polygon with columns:
 *          date (YYYY-MM), <selected SPEI bands...>
 *
 * Notes:
 *  - This script can export ALL SPEI bands (e.g., SPEI_01_month …)
 *    or just a subset (e.g., SPEI_12_month).
 *  - Units are standardized index values (dimensionless).
 ************************************************************/


// ====================== USER PARAMETERS ======================

// 1) Polygons layer: upload your cities.gpkg/geojson to Assets and paste ID
var POLYGONS_ASSET = 'users/<your-username>/cities311';   // <-- EDIT

// 2) Name/ID field in the polygons table (used in filenames)
var NAME_FIELD     = 'name';                               // <-- EDIT

// 3) Time window and months
var START_DATE     = '1981-01-01';     // SPEI starts ~1981 in many products
var END_DATE       = '2024-12-31';
var MONTH_START    = 1;                // set 1..12 for all months
var MONTH_END      = 12;

// 4) Which SPEI bands to export?
//    Leave EMPTY [] to export all available SPEI_* bands,
//    or specify a subset, e.g.: ['SPEI_12_month']
var SPEI_BANDS_WHITELIST = [];         // <-- EDIT (empty = all)

// 5) Export settings
var OUTPUT_FOLDER  = 'SPEI_TS_EXPORTS'; // Google Drive folder
var REDUCE_SCALE_M = 25000;             // ~25 km; adjust as desired
var MAX_PIXELS     = 1e13;

// 6) Run all polygons or a single one by name?
var RUN_ALL        = true;
var POLY_TO_RUN    = 'Los_Angeles';     // used only if RUN_ALL=false


// ========================= LOAD DATA =========================

var polys = ee.FeatureCollection(POLYGONS_ASSET);
var fc    = RUN_ALL ? polys : polys.filter(ee.Filter.eq(NAME_FIELD, POLY_TO_RUN));

var dateFilter  = ee.Filter.date(START_DATE, END_DATE);
var monthFilter = ee.Filter.calendarRange(MONTH_START, MONTH_END, 'month');

var speiColRaw = ee.ImageCollection('CSIC/SPEI/2_10')
  .filter(dateFilter)
  .filter(monthFilter)
  .sort('system:time_start');

print('Polygons to process:', fc.size());
print('SPEI monthly images selected:', speiColRaw.size());

if (speiColRaw.size().eq(0)) {
  print('WARNING: SPEI collection is empty with current filters.');
}

// Determine which bands to select
var firstImg    = ee.Image(speiColRaw.first());
var allBands    = ee.List(ee.Algorithms.If(firstImg, firstImg.bandNames(), ee.List([])));
var speiBands   = allBands.filter(ee.Filter.stringStartsWith('item', 'SPEI')); // e.g., SPEI_12_month
var selectedBands = (SPEI_BANDS_WHITELIST.length > 0)
  ? ee.List(SPEI_BANDS_WHITELIST)
  : speiBands;

print('SPEI bands to export:', selectedBands);


// ====================== EXPORT PER POLYGON ===================

var list = fc.toList(fc.size());
var N    = list.size().getInfo();

for (var i = 0; i < N; i++) {
  var feat = ee.Feature(list.get(i));
  var geom = ee.FeatureCollection([feat]).geometry();
  var name = ee.String(feat.get(NAME_FIELD));
  var safe = name.replace(' ', '_').replace('/', '_'); // filename-safe

  var speiCol = speiColRaw.select(selectedBands);

  // Build monthly rows: reduce means over polygon for the selected bands
  var rows = speiCol.map(function(img){
    var dateStr = ee.Date(img.get('system:time_start')).format('YYYY-MM');
    var stats = img.reduceRegion({
      reducer:   ee.Reducer.mean(),
      geometry:  geom,
      scale:     REDUCE_SCALE_M,
      maxPixels: MAX_PIXELS,
      bestEffort: true
    });
    return ee.Feature(null, stats.set('date', dateStr));
  })
  // Ensure at least one non-null band per row; check the first selected band
  .filter(ee.Filter.notNull([ee.String(selectedBands.get(0))]));

  // Build column order: date first, then all selected bands
  var cols = ee.List(['date']).cat(selectedBands);
  var desc = ee.String('SPEI_TS_').cat(safe).getInfo();

  Export.table.toDrive({
    collection:    rows,
    description:   desc,
    fileNamePrefix: desc,
    folder:        OUTPUT_FOLDER,
    fileFormat:    'CSV',
    selectors:     cols.getInfo()
  });

  print('Queued export:', desc);
}


// ============================ NOTES ==========================
// • Set SPEI_BANDS_WHITELIST to ['SPEI_12_month'] if you only want that timescale.
// • Keep REDUCE_SCALE_M consistent with the dataset resolution you need.
// • For large batches, run subsets (filter fc) to avoid task limits.
// • Post-process in Python/R/Matlab to join with other climate tables by 'date'.
