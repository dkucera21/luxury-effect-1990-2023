// ================== LOAD BOUNDARIES (repository-friendly) ==================
// USER EDITS THESE:
var CITIES_ASSET  = 'users/<your-username>/cities311';   // <- after uploading cities.gpkg/geojson
var NAME_FIELD    = 'name';                               // attribute that holds the city name/ID
var OUTPUT_FOLDER = 'GLOBAL_CITIES_LANDSAT_PIPELINE';         // Google Drive export folder

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

  // Plain JS strings for export descriptions
  var NDVIname = ee.String('NDVI_').cat(safeName).getInfo();
  var LSTname  = ee.String('LST_').cat(safeName).getInfo();

  print('Processing:', safeName);

  var RESAMPLE_SCALE = 30;
  var commonCrs = 'EPSG:4326';
  
  // ========== LANDSAT IMAGERY PRE-PROCESSING ==========
  // Define band names
  var bname = ['BLUE','GREEN','RED','NIR','SWIR1','SWIR2','ST','QA_PIXEL','QA_RADSAT','ST_QA'];
  var bl4_7 = ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7','ST_B6','QA_PIXEL','QA_RADSAT','ST_QA'];
  var bl8_9 = ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','ST_B10','QA_PIXEL','QA_RADSAT','ST_QA'];
  
  // Rename bands for consistency
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
    var saturationMask = image.select('QA_RADSAT').eq(0);
    return image.updateMask(mask).updateMask(saturationMask)
                .copyProperties(image, ["system:time_start"]);
  }
  
  // For these collections we want bands needed for indices.
  var bandsNeeded = ['RED','NIR','ST','QA_PIXEL','QA_RADSAT'];
  
  var collection_L8 = rename_bands(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'), bl8_9)
    .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
    .filter(ee.Filter.eq('IMAGE_QUALITY_TIRS', 9))
    .filter(ee.Filter.eq('IMAGE_QUALITY_OLI', 9))
    .filter(ee.Filter.eq('COLLECTION_CATEGORY', 'T1'))
    .filterBounds(shp)
    .map(function(image) { return image.clip(shp); })
    .map(maskLANDSAT);
  
  var collection_L9 = rename_bands(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'), bl8_9)
    .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
    .filter(ee.Filter.eq('IMAGE_QUALITY_TIRS', 9))
    .filter(ee.Filter.eq('IMAGE_QUALITY_OLI', 9))
    .filter(ee.Filter.eq('COLLECTION_CATEGORY', 'T1'))
    .filterBounds(shp)
    .map(function(image) { return image.clip(shp); })
    .map(maskLANDSAT);
  
  var collection_L7 = rename_bands(ee.ImageCollection("LANDSAT/LE07/C02/T1_L2"), bl4_7)
    .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
    .filterDate('1999-01-01', '2003-05-31') //Ignore L7 after scan line corrector failure 
    .filter(ee.Filter.eq('IMAGE_QUALITY', 9))
    .filter(ee.Filter.eq('COLLECTION_CATEGORY', 'T1'))
    .filterBounds(shp)
    .map(function(image) { return image.clip(shp); })
    .map(maskLANDSAT);
  
  var collection_L5 = rename_bands(ee.ImageCollection("LANDSAT/LT05/C02/T1_L2"), bl4_7)
    .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
    .filter(ee.Filter.eq('IMAGE_QUALITY', 9))
    .filter(ee.Filter.eq('COLLECTION_CATEGORY', 'T1'))
    .filterBounds(shp)
    .map(function(image) { return image.clip(shp); })
    .map(maskLANDSAT);
  
  var collection_L4 = rename_bands(ee.ImageCollection("LANDSAT/LT04/C02/T1_L2"), bl4_7)
    .filter(ee.Filter.eq('PROCESSING_LEVEL', 'L2SP'))
    .filter(ee.Filter.eq('IMAGE_QUALITY', 9))
    .filter(ee.Filter.eq('COLLECTION_CATEGORY', 'T1'))
    .filterBounds(shp)
    .map(function(image) { return image.clip(shp); })
    .map(maskLANDSAT);
    
  // Merge collections into sensor groups and overall collection for LST.
  var collection_L4_L7 = collection_L4.merge(collection_L5).merge(collection_L7)
    .filter(ee.Filter.calendarRange(1990,2012,'year'))
    .filter(ee.Filter.calendarRange(5,9,'month'));
  
  var collection_L8_L9 = collection_L8.merge(collection_L9)
    .filter(ee.Filter.calendarRange(2012,2024,'year'))
    .filter(ee.Filter.calendarRange(5,9,'month'));
    
  var collection_ST = collection_L4.merge(collection_L5)
                          .merge(collection_L7)
                          .merge(collection_L8)
                          .merge(collection_L9)
                          .filter(ee.Filter.calendarRange(1990,2024,'year'))
                          .filter(ee.Filter.calendarRange(5,9,'month'))
                          .map(function(img){ 
        return img.select(['ST']).copyProperties(img, img.propertyNames());
    });
    
  // ===== 1) MOSAICKING BY DATE ==========
function mosaicByDate(imcol) {
  // get list of all unique acquire‐dates
  var uniqueDates = imcol
    .distinct('DATE_ACQUIRED')
    .aggregate_array('DATE_ACQUIRED');
  
  // build one mosaic per date
  var mosaics = ee.ImageCollection(
    uniqueDates.map(function(d) {
      // d is a string "YYYY-MM-dd"
      var daily = imcol.filter(ee.Filter.eq('DATE_ACQUIRED', d));
      var first = daily.first();
      return daily
        .mean()
        .set({
          // preserve date/time metadata
          'system:time_start': first.get('system:time_start'),
          'DATE_ACQUIRED':    first.get('DATE_ACQUIRED')
        })
        .copyProperties(first, first.propertyNames());
    })
  );
  return mosaics;
}

// do the mosaicking on each raw collection
var mosaicked_L4_L7 = mosaicByDate(collection_L4_L7);
var mosaicked_L8_L9 = mosaicByDate(collection_L8_L9);
var mosaicked_ST   = mosaicByDate(collection_ST);

// ===== 2) DYNAMIC MASKING ON MOSAICS - Filters based on the remaining pixels after masking
// compared to all valid pixels for the ROI. Dates with >X% masked pixels are excluded from the analysis. =====
// helper: count non‐masked RED pixels over your city geometry

function addValidCount(image) {
  var valid = ee.Number(
    image.select('RED').reduceRegion({
      reducer:    ee.Reducer.count(),
      geometry:   shp.geometry(),
      scale:      30,
      maxPixels:  1e9
    }).get('RED')
  );
  return image.set('validCount', valid);
}


function addValidCountLST(image) {
  var valid = ee.Number(
    image.select('ST').reduceRegion({
      reducer:    ee.Reducer.count(),
      geometry:   shp,
      scale:      30,
      maxPixels:  1e9
    }).get('ST')
  );
  return image.set('validCount', valid);
}

// apply to each mosaicked collection
var fracKeep = 0.85; //Only keeps dates that have no more than 15% pixels that are masked within the ROI.
// apply to each mosaicked collection
var valid_L4_L7 = mosaicked_L4_L7.map(addValidCount);
var maxL4L7    = ee.Number(valid_L4_L7.aggregate_max('validCount'));
var filtered_L4_L7 = valid_L4_L7.filter(
  ee.Filter.gte('validCount', maxL4L7.multiply(fracKeep))
);

var valid_L8_L9 = mosaicked_L8_L9.map(addValidCount);
var maxL8L9    = ee.Number(valid_L8_L9.aggregate_max('validCount'));
var filtered_L8_L9 = valid_L8_L9.filter(
  ee.Filter.gte('validCount', maxL8L9.multiply(fracKeep))
);

var valid_ST = mosaicked_ST.map(addValidCountLST);
var maxST    = ee.Number(valid_ST.aggregate_max('validCount'));
var filtered_ST = valid_ST.filter(
  ee.Filter.gte('validCount', maxST.multiply(fracKeep))
);

  // ===== LST PROCESSING =====
  var convert_LST = function(img) {
    return img.multiply(0.00341802).add(149).subtract(273.15)
              .set('DATE_ACQUIRED', img.get('DATE_ACQUIRED'));
  };
  var maskZeroValuesAfterConversion_LST = function(image) {
    var mask = image.gt(10); // pixels >10°C retained
    return image.updateMask(mask).copyProperties(image, image.propertyNames());
  };
  var LST_converted = filtered_ST.map(convert_LST);
  var LST_ALL_DATES = LST_converted.map(maskZeroValuesAfterConversion_LST);

  
  // ===== DERIVE EVI, SAVI, AND NDVI ==========
  function addEVI(img) {
    var evi = img.expression(
      '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
        'NIR': img.select('NIR'),
        'RED': img.select('RED'),
        'BLUE': img.select('BLUE')
      }).rename('EVI')
      .set('DATE_ACQUIRED', img.get('DATE_ACQUIRED'));
    return img.addBands(evi);
  }
  
  function addSAVI(img) {
    var savi = img.expression(
      '((NIR - RED) / (NIR + RED + 0.5)) * (1.5)', {
        'NIR': img.select('NIR'),
        'RED': img.select('RED')
      }).rename('SAVI')
      .set('DATE_ACQUIRED', img.get('DATE_ACQUIRED'));
    return img.addBands(savi);
  }
  
  // Coefficients for harmonizing (NDVI/EVI/SAVI) bands (for Landsat 4-7)
  var coefficients = {
    itcps: ee.Image.constant([0.0003, 0.0061, 0.0412]),
    slopes: ee.Image.constant([0.8474, 0.9047, 0.8462])
  };
  
  function scaleToReflectance(img) {
    return img.select(['BLUE','RED','NIR']).multiply(0.0000275).subtract(0.2)
              .set('DATE_ACQUIRED', img.get('DATE_ACQUIRED'));
  }
  
  function Harmonize(img) {
    return img.select(['BLUE','RED','NIR'])
              .multiply(coefficients.slopes)
              .add(coefficients.itcps)
              .set('DATE_ACQUIRED', img.get('DATE_ACQUIRED'));
  }
  
  function addNDVI(img) {
    var ndvi = img.normalizedDifference(['NIR', 'RED']).rename('NDVI')
                 .set('DATE_ACQUIRED', img.get('DATE_ACQUIRED'));
    return img.addBands(ndvi);
  }
  
  // Apply scaling and harmonization to each sensor group.
  var scaled_L4_L7 = filtered_L4_L7.map(scaleToReflectance);
  var HARMONIZED = scaled_L4_L7.map(Harmonize);
  var scaled_L8_L9 = filtered_L8_L9.map(scaleToReflectance);
    // Map.addLayer(scaled_L4_L7.median(),{}, 'BOUND_city')
  print(scaled_L4_L7)
  
  // Compute indices for each group.
  var lsat_NDVI_EVI_SAVI = HARMONIZED.map(addNDVI).map(addEVI).map(addSAVI);
  var NDVI_EVI_SAVI_L8_L9 = scaled_L8_L9.map(addNDVI).map(addEVI).map(addSAVI);
  print(scaled_L8_L9,'scaled_L8_L9')

  // Merge into one collection and mask out invalid values. The merge requires data for both lsat_NDVI_EVI_SAVI and NDVI_EVI_SAVI_L8_L9.
  var FINAL_EVI_SAVI_NDVI = lsat_NDVI_EVI_SAVI.merge(NDVI_EVI_SAVI_L8_L9);
  
  function maskInvalidValues(image) {
    var validNDVI = image.select('NDVI').gt(-1).and(image.select('NDVI').lt(1));
    var validSAVI = image.select('SAVI').gt(-1).and(image.select('SAVI').lt(1));
    return image.updateMask(validNDVI)
                .updateMask(validSAVI)
                .set('DATE_ACQUIRED', image.get('DATE_ACQUIRED'));
  }
  
  var FINAL_MASKED_EVI_SAVI_NDVI = FINAL_EVI_SAVI_NDVI.map(maskInvalidValues);
  var finalCollection = FINAL_MASKED_EVI_SAVI_NDVI.select(['NDVI', 'EVI', 'SAVI']);
  print(finalCollection,'finalCollection')
  
  var LST_sorted = LST_ALL_DATES.sort('DATE_ACQUIRED');
  var finalCollection_sorted = finalCollection.sort('DATE_ACQUIRED');
  
  // IMPORTANT: do composites on FLOAT, not on scaled Ints
// Keep only NDVI band for composites
var NDVI_coll = finalCollection.select('NDVI');   // values in [-1, 1]
var LST_coll  = LST_ALL_DATES;                    // °C (you already subtracted 273.15)

  // ===== SCALING FOR EXPORT =====
  function scaleAndConvertToIntNDVI(image) {
    return image.multiply(10000).toInt16();
  }
  function scaleAndConvertToIntLST(image) {
    return image.multiply(100).toInt16();
  }
  
  var SVI_sorted_scaled = finalCollection_sorted.map(scaleAndConvertToIntNDVI);
  var LST_sorted_scaled = LST_sorted.map(scaleAndConvertToIntLST);
  
  var resampled_SVI = SVI_sorted_scaled.map(function(image) {
    return image.reproject({ crs: commonCrs, scale: RESAMPLE_SCALE });
  });
  var resampled_LST = LST_sorted_scaled.map(function(image) {
    return image.reproject({ crs: commonCrs, scale: RESAMPLE_SCALE });
  });
  
  var resampled_NDVI = resampled_SVI.select(['NDVI']);
  var resampled_SAVI = resampled_SVI.select(['SAVI']);
  var resampled_EVI = resampled_SVI.select(['EVI']);
  
  var NDVI_bands = resampled_NDVI.toBands();
  var SAVI_bands = resampled_SAVI.toBands();
  var EVI_bands = resampled_EVI.toBands();
  var LST_bands = resampled_LST.toBands();

  Export.image.toDrive({
  image: NDVI_bands,
  description: NDVIname,
  folder: OUTPUT_FOLDER,
  scale: 30,
  region: shp,
    fileFormat: 'GeoTIFF',
  formatOptions: {
    cloudOptimized: true
  }
});

Export.image.toDrive({
  image: LST_bands,
  description: LSTname,
  folder: OUTPUT_FOLDER,
  scale: 30,
  region: shp,
    fileFormat: 'GeoTIFF',
  formatOptions: {
    cloudOptimized: true
  }
});

   // ── BUILD DATE TABLES FROM system:index ──
  var ndviDates = resampled_NDVI.map(function(img) {
    // take last 8 characters of the system:index
    var dateStr = ee.String(img.get('system:index')).slice(-8);
    return ee.Feature(null, { 'DATE': dateStr });
  });

  var lstDates = resampled_LST.map(function(img) {
    var dateStr = ee.String(img.get('system:index')).slice(-8);
    return ee.Feature(null, { 'DATE': dateStr });
  });

  // 2) Export those tables to Drive as CSV
  Export.table.toDrive({
    collection: ndviDates,
    description: NDVIname + '_DATES_NDVI',
    folder: OUTPUT_FOLDER,
    fileNamePrefix: NDVIname + '_DATES_NDVI',
    fileFormat: 'CSV'
  });

  Export.table.toDrive({
    collection: lstDates,
    description: LSTname + '_DATES_LST',
    folder: OUTPUT_FOLDER,
    fileNamePrefix: LSTname + '_DATES_LST',
    fileFormat: 'CSV'
  });
  
}