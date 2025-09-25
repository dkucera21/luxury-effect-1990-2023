// ================== LOAD BOUNDARIES (repository-friendly) ==================
// USER EDITS THESE:
var CITIES_ASSET  = 'users/<your-username>/cities311';   // <- after uploading cities.gpkg/geojson
var NAME_FIELD    = 'name';                               // attribute that holds the city name/ID
var OUTPUT_FOLDER = 'LUXURY_CITIES_PIPELINE_1yr';         // Google Drive export folder

// Load all cities from the single FeatureCollection the user uploads
var cities = ee.FeatureCollection(CITIES_ASSET);

// Run for all cities, or set RUN_ALL=false and pick one by name
var RUN_ALL     = true;
var CITY_TO_RUN = 'Los_Angeles'; // used only if RUN_ALL === false

var fc = RUN_ALL ? cities : cities.filter(ee.Filter.eq(NAME_FIELD, CITY_TO_RUN));

// Convert to a client-side list for controlled exporting
var cityList = fc.toList(fc.size());
var N = cityList.size().getInfo();

for (var i = 0; i < N; i++) {
  var city = ee.Feature(cityList.get(i));
  var shp  = ee.FeatureCollection([city]);

  // Get a safe boundary name for filenames (replace spaces/slashes)
  var boundaryName = ee.String(city.get(NAME_FIELD));
  var safeName     = boundaryName.replace(' ', '_').replace('/', '_');

  var Boundary = assetList.get(i)
  var String_Bound = ee.String(Boundary)
  var String_Split = String_Bound.split('/')
  var Boundary_Name = String_Split.get(4)
  print(Boundary_Name)
  var Elevation_name = ee.String('Elevation_').cat(safeName).getInfo() //Can't use a computed object for the Layer Name.
  var Elevation_slope_name = ee.String('Elevation_Slope_').cat(safeName).getInfo()


var dataset = ee.Image('USGS/SRTMGL1_003')
var dataset = dataset.reproject('EPSG:4326', null, 30);

var clipped = dataset.clip(shp);
var elevation = clipped.select('elevation');
var slope = ee.Terrain.slope(elevation);

// Map.addLayer(slope, {min: 0, max: 60}, 'slope');
// Map.addLayer(elevation, {min: 0, max: 1000}, 'elevation');

Export.image.toDrive({
  image: elevation,
  description: Elevation_name,
  scale: 30,
  folder: 'LUX_HYPOS_USE',
  region: shp
});

Export.image.toDrive({
  image: slope,
  description: Elevation_slope_name,
  scale: 30,
  folder: OUTPUT_FOLDER,
  region: shp
});

//SRTM DIVERSITY
var SRTM_diversity_name = ee.String('SRTM_DIVERSITY_').cat(safeName).getInfo()
var SRTM_Diversity_1 = ee.Image("CSP/ERGo/1_0/Global/SRTM_topoDiversity")
var srtmTopographicDiversity = SRTM_Diversity_1.select('constant');

var DIVERSITY_CLIPPED = srtmTopographicDiversity.clip(shp)

Export.image.toDrive({
  image: DIVERSITY_CLIPPED,
  description: SRTM_diversity_name,
  folder: OUTPUT_FOLDER ,
  scale: 270,
  region: shp
});

}