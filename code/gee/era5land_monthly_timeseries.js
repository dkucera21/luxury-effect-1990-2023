/************************************************************
 * File: 04_era5land_monthly_timeseries.js
 * Purpose: Export ERA5-Land MONTHLY_AGGR variables as a
 *          monthly time series (one CSV per polygon).
 *
 * Dataset: ECMWF/ERA5_LAND/MONTHLY_AGGR (monthly aggregates)
 * Notes:
 *  - Variables include *_sum, *_max, *_mean variants per ERA5-Land.
 *  - This script filters months (e.g., May–Sep) and years.
 *  - Outputs one CSV per polygon with columns: date, <vars...>
 *  - Units/definitions follow the dataset; document in your README.
 ************************************************************/

// ====================== USER PARAMETERS ======================

// 1) Polygons layer (upload your cities.gpkg/geojson to Assets and paste its ID)
var POLYGONS_ASSET = 'users/<your-username>/cities311';   // <-- EDIT

// 2) Name/ID field in the polygons table (used in filenames)
var NAME_FIELD     = 'name';                               // <-- EDIT

// 3) Time window and months
var START_DATE     = '1989-01-01';
var END_DATE       = '2024-12-31';
var MONTH_START    = 5;    // May
var MONTH_END      = 9;    // September (inclusive)

// 4) ERA5-Land variables to export (adjust as needed)
var ERA5_VARS = [
  'surface_latent_heat_flux_sum',
  'leaf_area_index_high_vegetation',
  'leaf_area_index_low_vegetation',
  'soil_temperature_level_1_max',
  'surface_sensible_heat_flux_sum',
  'surface_net_solar_radiation_sum',
  'forecast_albedo',
  'soil_temperature_level_1',
  'soil_temperature_level_2',
  'temperature_2m',
  'skin_temperature',
  'dewpoint_temperature_2m',
  'evaporation_from_bare_soil_min',
  'evaporation_from_bare_soil_max',
  'potential_evaporation_max'
];

// 5) Export location and scale
var OUTPUT_FOLDER  = 'ERA5LAND_TS_EXPORTS'; // Google Drive folder
var REDUCE_SCALE_M = 10000;                 // ~0.1° ~ 9–11 km; adjust if desired

// 6) Run all polygons or a single one by name?
var RUN_ALL        = true;
var POLY_TO_RUN    = 'Los_Angeles';         // used only if RUN_ALL=false


// ========================= LOAD DATA =========================

var polys = ee.FeatureCollection(POLYGONS_ASSET);
var fc    = RUN_ALL ? polys : polys.filter(ee.Filter.eq(NAME_FIELD, POLY_TO_RUN));

var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
  .filterDate(START_DATE, END_DATE)
  .filter(ee.Filter.calendarRange(MONTH_START, MONTH_END, 'month'))
  .select(ERA5_VARS)
  .sort('system:time_start');

print('Polygons to process:', fc.size());
print('ERA5-Land monthly images selected:', era5.size());

// ====================== EXPORT PER POLYGON ===================

var list = fc.toList(fc.size());
var N    = list.size().getInfo();

for (var i = 0; i < N; i++) {
  var feat = ee.Feature(list.get(i));
  var geom = ee.FeatureCollection([feat]).geometry();
  var name = ee.String(feat.get(NAME_FIELD));
  var safe = name.replace(' ', '_').replace('/', '_'); // filename-safe

  // reduceRegion over each month; keep date as YYYY-MM
  var ts = era5.map(function(img) {
    var dateStr = ee.Date(img.get('system:time_start')).format('YYYY-MM');
    var stats = img.reduceRegion({
      reducer:   ee.Reducer.mean(),
      geometry:  geom,
      scale:     REDUCE_SCALE_M,
      maxPixels: 1e13,
      bestEffort: true
    });
    return ee.Feature(null, stats.set('date', dateStr));
  })
  // keep only rows that have data for at least one variable
  .filter(ee.Filter.notNull(ERA5_VARS.slice(0, 1)));

  var desc = ee.String('ERA5LAND_TS_').cat(safe).getInfo();

  Export.table.toDrive({
    collection:  ts,
    description: desc,
    fileNamePrefix: desc,
    folder:      OUTPUT_FOLDER,
    fileFormat:  'CSV',
    selectors:   ['date'].concat(ERA5_VARS)
  });

  print('Queued export:', desc);
}

// ============================ NOTES ==========================
// • If you prefer a single combined CSV for all polygons, map over the
//   polygons and append a 'name' attribute before exporting one table.
// • Some variables are in K (temperatures); convert to °C downstream if needed.
// • Sums (e.g., flux/radiation sums) are monthly aggregates per dataset docs.
// • Large batches can hit task limits; run subsets (filter fc) if necessary.
