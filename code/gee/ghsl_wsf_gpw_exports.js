/************************************************************
 * File: 08_ghsl_wsf_gpw_exports.js
 * Purpose: Export GHSL/WSF/GPW layers per polygon with
 *          explicit epoch/year selections and tidy filenames.
 *
 * Datasets used:
 *  • DLR/WSF/WSF2015/v1                    (settlement footprint, ~10 m)
 *  • JRC/GHSL/P2016/BUILT_LDSMT_GLOBE_V1   (built-up, ~38 m)
 *  • JRC/GHSL/P2023A/GHS_BUILT_C           (built characteristics, epochs)
 *  • JRC/GHSL/P2023A/GHS_BUILT_V           (built volume total, epochs)
 *  • JRC/GHSL/P2023A/GHS_BUILT_S           (built surface total/nres, epochs)
 *  • JRC/GHSL/P2023A/GHS_POP               (population count, epochs)
 *  • JRC/GHSL/P2023A/GHS_SMOD              (degree of urbanization, epochs)
 *  • CIESIN/GPWv411/GPW_UNWPP-Adjusted_Population_Density (years)
 *  • CIESIN/GPWv411/GPW_Population_Count   (years)
 *
 * Outputs: One GeoTIFF per layer per epoch/year per polygon.
 * Notes:
 *  - Filenames include <LAYER>_<EPOCHorYEAR>_<CITYNAME>.
 *  - All exports use EPSG:4326; include units/meaning in your README.
 ************************************************************/


// ====================== USER PARAMETERS ======================

// Polygons (upload your cities.gpkg/geojson; paste Asset ID)
var POLYGONS_ASSET = 'users/<your-username>/cities311';  // <-- EDIT
var NAME_FIELD     = 'name';                              // <-- EDIT

// Epochs for GHSL P2023A (must match the collection 'epoch' property)
var GHSL_EPOCHS = [1975, 1990, 2000, 2015, 2018, 2020];

// GPW years available: 2000, 2005, 2010, 2015, 2020
var GPW_YEARS = [2000, 2005, 2010, 2015, 2020];

// Export settings
var OUTPUT_FOLDER = 'GHSL_WSF_GPW_EXPORTS';
var MAX_PIXELS    = 1e13;

// Recommended pixel scales (meters)
var SCALE_WSF2015        = 10;       // DLR/WSF 2015
var SCALE_BUILT_LDSMT    = 38;       // GHSL P2016
var SCALE_GHSL_FINE      = 10;       // GHS_BUILT_C (native per epoch product; 10 m or resampling)
var SCALE_GHSL_DEFAULT   = 100;      // GHS_BUILT_V / GHS_BUILT_S / GHS_POP (~100 m)
var SCALE_SMOD           = 1000;     // GHS_SMOD (~1 km)
var SCALE_GPW            = 927.67;   // GPWv4 (~30 arcsec)

// Run all polygons, or just one
var RUN_ALL     = true;
var POLY_TO_RUN = 'Los_Angeles';     // used only if RUN_ALL=false


// ========================= LOAD DATA =========================

var polys = ee.FeatureCollection(POLYGONS_ASSET);
var fc    = RUN_ALL ? polys : polys.filter(ee.Filter.eq(NAME_FIELD, POLY_TO_RUN));

// Single images
var WSF2015   = ee.Image('DLR/WSF/WSF2015/v1');                 // settlement footprint (binary)
var BUILT16   = ee.Image('JRC/GHSL/P2016/BUILT_LDSMT_GLOBE_V1') // built-up (age/detect); band 'built'
                 .select('built');

// GHSL 2023A image collections
var GHS_BUILT_C = ee.ImageCollection('JRC/GHSL/P2023A/GHS_BUILT_C');   // band: 'built_characteristics'
var GHS_BUILT_V = ee.ImageCollection('JRC/GHSL/P2023A/GHS_BUILT_V');   // band: 'built_volume_total'
var GHS_BUILT_S = ee.ImageCollection('JRC/GHSL/P2023A/GHS_BUILT_S');   // bands: 'built_surface', 'built_surface_nres'
var GHS_POP     = ee.ImageCollection('JRC/GHSL/P2023A/GHS_POP');       // band: 'population_count'
var GHS_SMOD    = ee.ImageCollection('JRC/GHSL/P2023A/GHS_SMOD');      // band: 'smod_code'

// GPWv4 image collections (years)
var GPW_DENS = ee.ImageCollection('CIESIN/GPWv411/GPW_UNWPP-Adjusted_Population_Density')
                 .select('unwpp-adjusted_population_density');
var GPW_POP  = ee.ImageCollection('CIESIN/GPWv411/GPW_Population_Count')
                 .select('population_count');


// ====================== HELPERS ===============================

function safeName(s) { return ee.String(s).replace(' ', '_').replace('/', '_'); }

function exportImage(img, prefix, citySafe, scale, regionGeom) {
  var desc = ee.String(prefix).cat('_').cat(citySafe).getInfo();
  Export.image.toDrive({
    image: img.toFloat(),
    description: desc,
    fileNamePrefix: desc,
    folder: OUTPUT_FOLDER,
    region: regionGeom,
    crs: 'EPSG:4326',
    scale: scale,
    maxPixels: MAX_PIXELS,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true}
  });
  print('Queued export:', desc);
}

// Filter GHSL 2023A by epoch property (int or string tolerant)
function filterByEpoch(ic, epoch) {
  var eStr = ee.String(ee.Number(epoch).format());
  return ic.filter(ee.Filter.or(
    ee.Filter.eq('epoch', epoch),
    ee.Filter.eq('epoch', eStr)
  ));
}


// ==================== EXPORT PER POLYGON ======================

var list = fc.toList(fc.size());
var N    = list.size().getInfo();

print('Polygons to process:', N);

for (var i = 0; i < N; i++) {
  var feat = ee.Feature(list.get(i));
  var geom = ee.FeatureCollection([feat]).geometry();
  var name = ee.String(feat.get(NAME_FIELD));
  var city = safeName(name);

  // ---------- WSF 2015 (single image) ----------
  exportImage(WSF2015.clip(geom), 'WSF2015', city, SCALE_WSF2015, geom);

  // ---------- GHSL P2016 BUILT_LDSMT (single image; band "built") ----------
  exportImage(BUILT16.clip(geom), 'GHSL_BUILT_LDSMT_P2016', city, SCALE_BUILT_LDSMT, geom);

  // ---------- GHSL P2023A epoch-based exports ----------
  GHSL_EPOCHS.forEach(function(epoch) {
    var epochStr = ee.String(ee.Number(epoch).format());
    var tag = '_' + epochStr.getInfo();

    // Built characteristics (single-band mosaic per epoch)
    var builtC = filterByEpoch(GHS_BUILT_C, epoch).select('built_characteristics').mosaic().clip(geom);
    exportImage(builtC, 'GHS_BUILT_C' + tag, city, SCALE_GHSL_FINE, geom);

    // Built volume (single band)
    var builtV = filterByEpoch(GHS_BUILT_V, epoch).select('built_volume_total').mosaic().clip(geom);
    exportImage(builtV, 'GHS_BUILT_V' + tag, city, SCALE_GHSL_DEFAULT, geom);

    // Built surface (two bands: total + non-residential)
    var builtS_total = filterByEpoch(GHS_BUILT_S, epoch).select('built_surface').mosaic().clip(geom);
    exportImage(builtS_total, 'GHS_BUILT_S_total' + tag, city, SCALE_GHSL_DEFAULT, geom);

    var builtS_nres  = filterByEpoch(GHS_BUILT_S, epoch).select('built_surface_nres').mosaic().clip(geom);
    exportImage(builtS_nres, 'GHS_BUILT_S_nres' + tag, city, SCALE_GHSL_DEFAULT, geom);

    // Population (count)
    var pop = filterByEpoch(GHS_POP, epoch).select('population_count').mosaic().clip(geom);
    exportImage(pop, 'GHS_POP' + tag, city, SCALE_GHSL_DEFAULT, geom);

    // Degree of urbanization (smod_code)
    var smod = filterByEpoch(GHS_SMOD, epoch).select('smod_code').mosaic().clip(geom);
    exportImage(smod, 'GHS_SMOD' + tag, city, SCALE_SMOD, geom);
  });

  // ---------- GPWv4 (by year) ----------
  GPW_YEARS.forEach(function(y) {
    var yStr = ee.String(ee.Number(y).format());
    var tag = '_' + yStr.getInfo();

    var dens = GPW_DENS.filter(ee.Filter.eq('year', y)).first(); // one image per year
    var pop  = GPW_POP .filter(ee.Filter.eq('year', y)).first();

    // Guard against missing years in some regions
    dens = ee.Image(ee.Algorithms.If(dens, dens, ee.Image(0).rename('unwpp-adjusted_population_density')));
    pop  = ee.Image(ee.Algorithms.If(pop,  pop,  ee.Image(0).rename('population_count')));

    exportImage(dens.clip(geom), 'GPWv4_Density' + tag, city, SCALE_GPW, geom);
    exportImage(pop.clip(geom),  'GPWv4_PopCount' + tag, city, SCALE_GPW, geom);
  });
}


// ========================== NOTES ============================
// • We export per epoch/year per city to keep task sizes moderate.
// • If you prefer one multi-band export per layer, you can toBands()
//   within a single epoch (or year set), but beware very large assets.
// • Double-check dataset docs for band units & semantics and document
//   them in your repo README.
// • If you need fewer outputs, trim GHSL_EPOCHS / GPW_YEARS above.
// • If exports hit task limits, run subsets (filter polygons or epochs).
