/************************************************************
 * File: terraclimate_aridity_timeseries.js
 * Purpose: Export monthly aridity index (AI = PR / PET) time
 *          series (1990–2024) for each polygon in a user-
 *          supplied FeatureCollection (e.g., cities).
 * Inputs (user provides):
 *   - A single Asset (table) with all polygons and a name/id field
 *   - TerraClimate (public)
 * Outputs:
 *   - One CSV per polygon: columns = [date, AI]
 * Notes:
 *   - PR: precipitation (mm)   | band 'pr'
 *   - PET: potential ET (mm)   | band 'pet' (scale by 0.1 per TerraClimate doc)
 *   - AI = PR / PET
 *   - Months: configurable (e.g., May–Sep)
 *   - Native TerraClimate scale ~ 4638 m (EPSG:4326 grid)
 ************************************************************/


// ====================== USER PARAMETERS ======================

// 1) After uploading your polygons (e.g., cities.gpkg/geojson) to GEE,
//    copy the Asset ID and paste here:
var POLYGONS_ASSET = 'users/<your-username>/cities311';   // <-- EDIT

// 2) Attribute in your polygons table that holds the polygon name/ID:
var NAME_FIELD     = 'name';                               // <-- EDIT

// 3) Time window (inclusive of start, exclusive of end)
var DATE_START     = '1990-01-01';
var DATE_END       = '2024-12-31';

// 4) Months to include (set to [1,12] for all months)
var MONTH_START    = 5;   // May
var MONTH_END      = 9;   // Sep

// 5) Export location on Google Drive
var OUTPUT_FOLDER  = 'AI_timeseries_exports';              // <-- EDIT if desired

// 6) Run all polygons or a single one by name?
var RUN_ALL        = true;
var POLY_TO_RUN    = 'Los_Angeles';                        // used only if RUN_ALL=false


// ====================== LOAD INPUTS ==========================

// Load polygons
var polys = ee.FeatureCollection(POLYGONS_ASSET);
var fc = RUN_ALL ? polys : polys.filter(ee.Filter.eq(NAME_FIELD, POLY_TO_RUN));

// TerraClimate collections
var dateFilter  = ee.Filter.date(DATE_START, DATE_END);
var monthFilter = ee.Filter.calendarRange(MONTH_START, MONTH_END, 'month');

var terraclimate = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
                      .filter(dateFilter)
                      .filter(monthFilter);

// PET (scaled by 0.1 → mm)
var PETall = terraclimate.select('pet').map(function(img){
  // Scale PET and keep properties/time
  return img.multiply(0.1)
            .copyProperties(img, img.propertyNames());
});

// PR (mm)
var PRall  = terraclimate.select('pr');

// Build monthly AI images with a clean time stamp
var monthlyAI = PRall.map(function(prImg){
  var start = ee.Date(prImg.get('system:time_start'));
  // PET for the same month (Terraclimate is monthly; filter by the exact month)
  var petImg = PETall.filterDate(start, start.advance(1, 'month')).first();
  var ai = prImg.divide(petImg).rename('AI');
  return ai.set({
    'system:time_start': start.millis(),
    'year':  start.get('year'),
    'month': start.get('month')
  });
}).sort('system:time_start');


// ====================== EXPORT PER POLYGON ===================

// Convert to a client-side list to control the export loop
var list = fc.toList(fc.size());
var N    = list.size().getInfo();

print('Exporting AI time series for', N, 'polygons.');

for (var i = 0; i < N; i++) {
  var feat = ee.Feature(list.get(i));
  var geom = ee.FeatureCollection([feat]).geometry();
  var name = ee.String(feat.get(NAME_FIELD));
  var safe = name.replace(' ', '_').replace('/', '_');   // filename-safe

  // Reduce each monthly AI image over the polygon to a mean value
  var tsFeatures = monthlyAI.map(function(img){
    var y  = ee.Number(img.get('year'));
    var m  = ee.Number(img.get('month'));
    var ym = y.format()                             // YYYY
              .cat(ee.String(m.format()).padStart(2, '0')); // MM
    var meanAI = img.reduceRegion({
      reducer:   ee.Reducer.mean(),
      geometry:  geom,
      scale:     4638.3,     // TerraClimate native
      maxPixels: 1e13,
      bestEffort: true
    }).get('AI');

    return ee.Feature(null, {
      date: ym,              // e.g., 199005
      AI:   meanAI
    });
  }).filter(ee.Filter.notNull(['AI']));

  // Export CSV to Drive
  var desc = ee.String('AI_timeseries_').cat(safe).getInfo();
  Export.table.toDrive({
    collection: tsFeatures,
    description: desc,
    fileNamePrefix: desc,
    folder: OUTPUT_FOLDER,
    fileFormat: 'CSV',
    selectors: ['date', 'AI']
  });

  print('Queued export:', desc);
}


// ====================== NOTES ================================
// • If you need all months, set MONTH_START=1 and MONTH_END=12.
// • To aggregate to seasonal/annual AI inside GEE, replace the per-image
//   reduceRegion with a grouped reducer by year/month, or post-process
//   the CSVs in Python/MATLAB.
// • Make sure your polygons layer includes the NAME_FIELD specified above.
// • Large batches may require running in chunks (e.g., split by country).
