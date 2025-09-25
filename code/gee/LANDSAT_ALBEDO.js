// ================== LOAD BOUNDARIES (repository-friendly) ==================
// USER EDITS THESE:
var CITIES_ASSET  = 'users/<your-username>/cities311';   // <- after uploading cities.gpkg/geojson
var NAME_FIELD    = 'name';                               // attribute that holds the city name/ID
var OUTPUT_FOLDER = 'LUXURY_CITIES_PIPELINE_1yr';         // Google Drive export folder
var RESAMPLE_SCALE = 30;
var commonCrs = 'EPSG:4326'

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
  print('Processing:', safeName);

  var NAME_albedo = ee.String('ALBEDO_').cat(safeName).getInfo();
  var NAME_albedo_mean = ee.String('ALBEDO_').cat(safeName).getInfo();
  var DATES_albedo = ee.String('DATES_ALBEDO_').cat(safeName).getInfo();
  var mean_albedo_name = ee.String('ALBEDO_TABLE_').cat(safeName).getInfo();
  
  // Define band names for Landsat sensors
  var bname = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'ST', 'QA_PIXEL', 'QA_RADSAT', 'ST_QA'];
  var bl4_7 = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6', 'QA_PIXEL', 'QA_RADSAT', 'ST_QA'];
  var bl8 = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10', 'QA_PIXEL', 'QA_RADSAT', 'ST_QA'];

  // Function to rename bands for consistency
  var rename_bands = function(img, input) {
    return img.select(input, bname);
  };

  // Function to create a strict mask for Landsat images based on QA_PIXEL and QA_RADSAT
  function maskLANDSAT(image) {
    var qa = image.select('QA_PIXEL');
    var mask = qa.bitwiseAnd(1 << 6).neq(0) // Clear
      .and(qa.bitwiseAnd(1 << 3).eq(0)) // Cloud
      .and(qa.bitwiseAnd(1 << 4).eq(0)) // Cloud shadow
      .and(qa.bitwiseAnd(1 << 5).eq(0)) // Snow
      .and(qa.bitwiseAnd(1 << 7).eq(0)) // Water
      .and(qa.bitwiseAnd(1 << 1).eq(0)) // Dilated cloud
      .and(qa.bitwiseAnd(1 << 2).eq(0)) // Cirrus
      .and(qa.bitwiseAnd(parseInt('110000000', 2)).rightShift(8).lte(1)) // Cloud confidence: None or Low
      .and(qa.bitwiseAnd(parseInt('11000000000', 2)).rightShift(10).lte(1)) // Cloud shadow confidence: None or Low
      .and(qa.bitwiseAnd(parseInt('11000000000000', 2)).rightShift(12).lte(1)) // Snow/Ice confidence: None or Low
      .and(qa.bitwiseAnd(parseInt('1100000000000000', 2)).rightShift(14).lte(1)); // Cirrus confidence: None or Low

    var saturationMask = image.select('QA_RADSAT').eq(0); // Non-saturated pixels
    return image.updateMask(mask).updateMask(saturationMask)
                .copyProperties(image, ["system:time_start"]);
  }
  
  var bandsNeeded = ['BLUE', 'RED', 'NIR', 'SWIR1', 'SWIR2'];

  // Load Landsat collections and apply QA masking and band renaming.
  var collection_L8 = rename_bands(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'), bl8)
    .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
    .filter(ee.Filter.eq('IMAGE_QUALITY_TIRS', 9))
    .filter(ee.Filter.eq('IMAGE_QUALITY_OLI', 9))
    .filter(ee.Filter.eq('COLLECTION_CATEGORY', 'T1'))
    .filterBounds(shp)
    .map(function(image) { return image.clip(shp); })
    .map(maskLANDSAT)
    .select(bandsNeeded);

  var collection_L9 = rename_bands(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'), bl8)
    .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
    .filter(ee.Filter.eq('IMAGE_QUALITY_TIRS', 9))
    .filter(ee.Filter.eq('IMAGE_QUALITY_OLI', 9))
    .filter(ee.Filter.eq('COLLECTION_CATEGORY', 'T1'))
    .filterBounds(shp)
    .map(function(image) { return image.clip(shp); })
    .map(maskLANDSAT)
    .select(bandsNeeded);

  var collection_L7 = rename_bands(ee.ImageCollection("LANDSAT/LE07/C02/T1_L2"), bl4_7)
    .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
    // .filterDate('1999-01-01', '2003-05-31') //Skip dates after the scan line corrector failure
    .filter(ee.Filter.eq('IMAGE_QUALITY', 9))
    .filter(ee.Filter.eq('COLLECTION_CATEGORY', 'T1'))
    .filterBounds(shp)
    .map(function(image) { return image.clip(shp); })
    .map(maskLANDSAT)
    .select(bandsNeeded);

  var collection_L5 = rename_bands(ee.ImageCollection("LANDSAT/LT05/C02/T1_L2"), bl4_7)
    .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
    .filter(ee.Filter.eq('IMAGE_QUALITY', 9))
    .filter(ee.Filter.eq('COLLECTION_CATEGORY', 'T1'))
    .filterBounds(shp)
    .map(function(image) { return image.clip(shp); })
    .map(maskLANDSAT)
    .select(bandsNeeded);

  var collection_L4 = rename_bands(ee.ImageCollection("LANDSAT/LT04/C02/T1_L2"), bl4_7)
    .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
    .filter(ee.Filter.eq('IMAGE_QUALITY', 9))
    .filter(ee.Filter.eq('COLLECTION_CATEGORY', 'T1'))
    .filterBounds(shp)
    .map(function(image) { return image.clip(shp); })
    .map(maskLANDSAT)
    .select(bandsNeeded);
    
  // KEEPING LANDSAT 4-7 and Landsat 8-9 separate are necessary at this point due to inter-satellite calibrations.
  var collection_L4_L7 = collection_L4.merge(collection_L5).merge(collection_L7);
  var collection_L8_L9 = collection_L8.merge(collection_L9);

  // Function to scale Landsat bands
  function scaleLandsatBands(img) {
    return img.select(['BLUE', 'RED', 'NIR', 'SWIR1', 'SWIR2'])
              .multiply(0.0000275).subtract(0.2)
              .copyProperties(img, img.propertyNames());
  }

  // Coefficients for Landsat Spectral Harmonization
  var coefficients = {
    itcps: ee.Image.constant([0.0003, 0.0061, 0.0412, 0.0254, 0.0172]),
    slopes: ee.Image.constant([0.8474, 0.9047, 0.8462, 0.8937, 0.9071])
  };

  // Harmonize spectral bands of Landsat 4-7 with Landsat 8-9
  function Harmonize(img) {
    return img.select(['BLUE', 'RED', 'NIR', 'SWIR1', 'SWIR2'])
              .multiply(coefficients.slopes)
              .add(coefficients.itcps)
              .copyProperties(img, img.propertyNames());
  }
  var harmonized_L4_L7 = collection_L4_L7.map(scaleLandsatBands).map(Harmonize);
  var scaled_L8_L9 = collection_L8_L9.map(scaleLandsatBands)
                                   .select(['BLUE', 'RED', 'NIR', 'SWIR1', 'SWIR2']);

  // Combine harmonized collections and filter by date and season.
  var allLandsat = harmonized_L4_L7.merge(scaled_L8_L9)
                        .filterDate('1990-01-01', '2024-12-31')
                        .filter(ee.Filter.calendarRange(5,9, 'month'));
                        
  // ===== Dynamic Filtering Step Based on Valid Pixels =====
  // For each image, compute the number of valid (nonmasked) pixels in the geometry.
  // (Here we use the 'BLUE' band; any band with the proper mask will work.)
  var allLandsatWithValid = allLandsat.map(function(image) {
    var valid = ee.Number(
      image.select('BLUE').reduceRegion({
        reducer: ee.Reducer.count(),
        geometry: shp,
        scale: 30,
        maxPixels: 1e9
      }).get('BLUE')
    );
    return image.set('validCount', valid);
  });

  // Compute the reference valid pixel count.
  // Here we use the maximum valid count across the timeseries as the “full‐coverage” benchmark.
  var maxValidCount = ee.Number(allLandsatWithValid.aggregate_max('validCount'));
  // print('Maximum valid pixel count', maxValidCount);

  // Filter out images that have less than 80% of the max valid pixels.
  var filteredLandsat = allLandsatWithValid.filter(
    ee.Filter.gte('validCount', maxValidCount.multiply(0.80))
  );
  // print('Number of images after dynamic filtering', filteredLandsat.size());
  // ===== End Dynamic Filtering =====

  // Function to calculate albedo
  var Albedo = function(image) {
    var alb = image.expression(
      "((0.356*blue) + (0.130*red) + (0.373*nir) + (0.085*swir) + (0.072*swir2) - 0.018) / 1.016", {
        'blue': image.select('BLUE'),
        'red': image.select('RED'),
        'nir': image.select('NIR'),
        'swir': image.select('SWIR1'),
        'swir2': image.select('SWIR2')
      });
    alb = alb.clamp(-0.001, 1.001).toFloat(); // Ensuring consistent range and type
    return image.addBands(alb.rename('albedo'))
                .copyProperties(image, image.propertyNames());
  };
  
  // Function to mask albedo values outside the valid range [0, 1] and retain 'DATE_ACQUIRED'
  var maskInvalidAlbedo = function(image) {
    var validRangeMask = image.select('albedo').gte(0).and(image.select('albedo').lte(1));
    return image.updateMask(validRangeMask)
                .set('DATE_ACQUIRED', image.get('DATE_ACQUIRED'));
  };

  // Apply albedo calculation to the filtered collection
  var landsatWithAlbedo = filteredLandsat
    .map(Albedo)
    .map(maskInvalidAlbedo)
    .select(['albedo']);
  
  var AlbedoSorted = landsatWithAlbedo.sort('DATE_ACQUIRED');

  // Reduce the image collection to a single image (mean albedo)
  var albedo_mean = AlbedoSorted.select('albedo').mean();

  // Function to extract mean albedo value and include the associated date
  var extractMeanAlbedo = function(image) {
    var geometry = shp.geometry();
    var mean_albedo = ee.Image(image).reduceRegion({
      reducer: ee.Reducer.median(),
      geometry: shp,
      scale: 30,
      maxPixels: 1e9
    }).get('albedo');
    return ee.Feature(null, {
      'mean_albedo': mean_albedo,
      'DATE_ACQUIRED': image.get('DATE_ACQUIRED')
    });
  };

  // Apply the function to extract mean albedo and date
  var Mean_Albedo = AlbedoSorted.map(extractMeanAlbedo);

  var resampled_ALBEDO = AlbedoSorted.map(function(image) {
    return image.reproject({
      crs: commonCrs,
      scale: RESAMPLE_SCALE,
    });
  });
  
// print(resampled_ALBEDO,'resampled')
  // Export the median albedo image to Google Drive (GeoTIFF)
  // Export.image.toDrive({
  //   image: resampled_ALBEDO.median(),
  //   description: NAME_albedo_mean,
  //   scale: 30,
  //   region: shp,
  //   crs: commonCrs,
  //   folder: OUTPUT_FOLDER,
  //   fileFormat: 'GeoTIFF'
  // });

  // Export the mean albedo values and associated dates to a single CSV
  Export.table.toDrive({
    collection: Mean_Albedo, 
    description: mean_albedo_name, 
    fileFormat: "CSV", 
    folder: OUTPUT_FOLDER,
    selectors: ["DATE_ACQUIRED", "mean_albedo"]
  });

}