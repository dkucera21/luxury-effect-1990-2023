/************************************************************
 * File: 03_climate_terraclimate_era5.js
 * Purpose: Export seasonal (e.g., May–Sep) multi-year composites
 *          from TerraClimate and ERA5 for each polygon in a user-
 *          supplied FeatureCollection (e.g., cities).
 * Outputs: One GeoTIFF per variable per polygon.
 *
 * TerraClimate variables (native units as provided by dataset):
 *   vpd, srad, aet, def, pet, pr, tmmn, tmmx, vap
 * ERA5 monthly variables (K):
 *   mean_2m_air_temperature, maximum_2m_air_temperature,
 *   minimum_2m_air_temperature, dewpoint_2m_temperature
 *
 * Notes:
 * - This script exports MEDIANS across the selected months/years.
 * - Region is each polygon geometry; CRS = EPSG:4326.
 * - Scales: TerraClimate ~4638 m; ERA5 ~25 km (0.25°). Adjust if needed.
 ************************************************************/


// ====================== USER PARAMETERS ======================

// After uploading your polygons (e.g., cities.gpkg/geojson) to GEE,
// copy the Asset ID and paste here:
var POLYGONS_ASSET = 'users/<your-username>/cities311';   // <-- EDIT

// Attribute in your polygons table that holds the polygon name/ID:
var NAME_FIELD     = 'name';                               // <-- EDIT

// Time windows (inclusive start, exclusive end)
var TC_DATE_START  = '1990-01-01';
var TC_DATE_END    = '2024-12-31';  // TerraClimate

var ERA5_DATE_START = '2019-01-01';
var ERA5_DATE_END   = '2024-12-31'; // ERA5 Monthly

// Months to include (set 1..12 for all months)
var MONTH_START    = 5;   // May
var MONTH_END      = 9;   // Sep

// Export location on Google Drive
var OUTPUT_FOLDER  = 'CLIMATE_EXPORTS';                    // <-- EDIT if desired

// Run all polygons or a single one by name?
var RUN_ALL        = true;
var POLY_TO_RUN    = 'Los_Angeles';                        // used only if RUN_ALL=false

// Export scales (meters)
var TERRACLIMATE_SCALE = 4638.3;   // native TerraClimate grid
var ERA5_SCALE         = 25000;    // ~0.25°; adjust if preferred


// ====================== LOAD INPUTS ==========================

var polys = ee.FeatureCollection(POLYGONS_ASSET);
var fc    = RUN_ALL ? polys : polys.filter(ee.Filter.eq(NAME_FIELD, POLY_TO_RUN));

// Filters
var tcDateFilter   = ee.Filter.date(TC_DATE_START, TC_DATE_END);
var era5DateFilter = ee.Filter.date(ERA5_DATE_START, ERA5_DATE_END);
var monthFilter    = ee.Filter.calendarRange(MONTH_START, MONTH_END, 'month');

// TerraClimate (monthly)
var tc = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
            .filter(tcDateFilter)
            .filter(monthFilter);

// Collections per variable (no scaling applied here; keep native units)
var TC_VPD   = tc.select('vpd');
var TC_SRAD  = tc.select('srad');
var TC_AET   = tc.select('aet');
var TC_DEF   = tc.select('def');
var TC_PET   = tc.select('pet');
var TC_PR    = tc.select('pr');
var TC_TMIN  = tc.select('tmmn');
var TC_TMAX  = tc.select('tmmx');
var TC_VAP   = tc.select('vap');

// ERA5 Monthly (K)
var ERA5 = ee.ImageCollection('ECMWF/ERA5/MONTHLY')
              .filter(era5DateFilter)
              .filter(monthFilter);

var ERA5_TMEAN = ERA5.select('mean_2m_air_temperature');
var ERA5_TMAX  = ERA5.select('maximum_2m_air_temperature');
var ERA5_TMIN  = ERA5.select('minimum_2m_air_temperature');
var ERA5_DEW   = ERA5.select('dewpoint_2m_temperature');


// ==================== EXPORT PER POLYGON =====================

var list = fc.toList(fc.size());
var N    = list.size().getInfo();

print('Exporting climate composites for', N, 'polygons.');

for (var i = 0; i < N; i++) {
  var feat = ee.Feature(list.get(i));
  var geom = ee.FeatureCollection([feat]).geometry();
  var name = ee.String(feat.get(NAME_FIELD));
  var safe = name.replace(' ', '_').replace('/', '_'); // filename-safe

  // Clip & median composites over period/months
  function med(ic) { return ic.map(function(img){return img.clip(geom);}).median(); }

  var VPD_med   = med(TC_VPD);
  var SRAD_med  = med(TC_SRAD);
  var AET_med   = med(TC_AET);
  var DEF_med   = med(TC_DEF);
  var PET_med   = med(TC_PET);
  var PR_med    = med(TC_PR);
  var TMIN_med  = med(TC_TMIN);
  var TMAX_med  = med(TC_TMAX);
  var VAP_med   = med(TC_VAP);

  var TMEAN_med = med(ERA5_TMEAN);
  var TMAXe_med = med(ERA5_TMAX);
  var TMINe_med = med(ERA5_TMIN);
  var DEW_med   = med(ERA5_DEW);

  // Helper to export images
  function exportImage(img, prefix, scale){
    var desc = ee.String(prefix).cat('_').cat(safe).getInfo();
    Export.image.toDrive({
      image: img.toFloat(),          // keep as float; document units in README
      description: desc,
      fileNamePrefix: desc,
      folder: OUTPUT_FOLDER,
      region: geom,
      crs: 'EPSG:4326',
      scale: scale,
      maxPixels: 1e13,
      fileFormat: 'GeoTIFF',
      formatOptions: {cloudOptimized: true}
    });
    print('Queued export:', desc);
  }

  // ---- TerraClimate exports (median May–Sep 1990–2024) ----
  exportImage(VPD_med,  'TC_VPD_MJJA1990_2024',  TERRACLIMATE_SCALE);
  exportImage(SRAD_med, 'TC_SRAD_MJJA1990_2024', TERRACLIMATE_SCALE);
  exportImage(AET_med,  'TC_AET_MJJA1990_2024',  TERRACLIMATE_SCALE);
  exportImage(DEF_med,  'TC_DEF_MJJA1990_2024',  TERRACLIMATE_SCALE);
  exportImage(PET_med,  'TC_PET_MJJA1990_2024',  TERRACLIMATE_SCALE);
  exportImage(PR_med,   'TC_PR_MJJA1990_2024',   TERRACLIMATE_SCALE);
  exportImage(TMIN_med, 'TC_TMIN_MJJA1990_2024', TERRACLIMATE_SCALE);
  exportImage(TMAX_med, 'TC_TMAX_MJJA1990_2024', TERRACLIMATE_SCALE);
  exportImage(VAP_med,  'TC_VAP_MJJA1990_2024',  TERRACLIMATE_SCALE);

  // ---- ERA5 exports (median May–Sep 2019–2024) ----
  exportImage(TMEAN_med,'ERA5_TMEAN_MJJA2019_2024', ERA5_SCALE);
  exportImage(TMAXe_med,'ERA5_TMAX_MJJA2019_2024',  ERA5_SCALE);
  exportImage(TMINe_med,'ERA5_TMIN_MJJA2019_2024',  ERA5_SCALE);
  exportImage(DEW_med,  'ERA5_DEW_MJJA2019_2024',   ERA5_SCALE);
}


// ========================== NOTES ============================
// • If you need annual or per-year stats instead of a multi-year median,
//   group by year/month or export per-year images.
// • TerraClimate variable units/scales follow dataset docs; if you prefer
//   scaled-to-physical units (e.g., temperatures), convert explicitly and
//   document it in your data README.
// • ERA5 variables are in Kelvin. Convert to °C in post-processing if needed.
// • Consider batching polygons (e.g., by country) if you hit task limits.
