// =============================================================================
// nrt-monitoring-demo.js
// Operational NRT forest change detection — demonstration script
//
// Paper:   Operational event-based forest change detection using
//          Sentinel-2 time series
// Authors: Viktor Myroniuk, Yevheniy Khan, Matthew J. Gregory,
//          Andrii Terentiev, Ihor Mischenko, Oleksandr Lesnik
// Version: 2026-09-01
//
// Description:
//   This script reproduces the near-real-time (NRT) forest change detection
//   workflow described in the paper using two pre-computed CCDC assets for a
//   demonstration tile in central Ukraine (50.73°N, 29.47°E).
//
//   The algorithm:
//   1. Builds a cloud-free Sentinel-2 collection within a monitoring window
//      ending at `currentDate`.
//   2. For each image, generates a phenologically-matched synthetic reference
//      by evaluating the CCDC harmonic model one year prior (`hindcastYears`).
//   3. Computes the Relativized differenced Normalized Burn Ratio (RdNBR)
//      between the synthetic reference and the observed image.
//   4. Applies an N-consecutive-flag confirmation rule (`nConsecutive`) to
//      separate confirmed changes from unconfirmed candidates.
//
//   Adjust `currentDate` (range: 2025-07-15 to 2025-08-28) to step through
//   the clearcut detection progression illustrated in the paper.
//
// Public assets used:
//   CCDC segments: projects/ee-victormyroniuk/assets/NRT/DEMO/
//                  N5050E29_2019-01-01_2025-11-01_ccdc
//   FNF mask:      projects/ee-victormyroniuk/assets/NRT/DEMO/
//                  N5050E29_2019-01-01_2025-11-01_cl_2024
//
// CCDC utilities (ccdc-utils.js) adapted from:
//   temporalSegmentation library by Daniel Wiell
//   https://code.earthengine.google.com/?accept_repo=users/wiell/temporalSegmentation
// =============================================================================


// ---------------------------------------------------------------------------
// 0. Dependencies
// ---------------------------------------------------------------------------
// Copy ccdc-utils.js to your GEE repository and update the path below.
var ccdc = require('users/victormyroniuk/<your-repo>:ccdc-utils');


// ---------------------------------------------------------------------------
// 1. Assets
// ---------------------------------------------------------------------------

// CCDC multi-segment array image (2019-01-01 to 2025-11-01, 20 m, S-2)
var ccdcImage    = ee.Image('projects/ee-victormyroniuk/assets/NRT/DEMO/N5050E29_2019-01-01_2025-11-01_ccdc');
var ccdcSegments = ccdc.Segments(ccdcImage);

// Forest / non-forest mask for 2024 (class 1 = forest)
var fnf = ee.Image('projects/ee-victormyroniuk/assets/NRT/DEMO/N5050E29_2019-01-01_2025-11-01_cl_2024')
            .eq(1).selfMask();


// ---------------------------------------------------------------------------
// 2. Area of interest and reference geometry
// ---------------------------------------------------------------------------

// Demonstration tile (approx. 1.5 × 1.0 km)
var aoi = ee.Geometry.Rectangle([29.4619, 50.7223, 29.4806, 50.7321]);

// Perimeter of the test clearcut used for visual validation in the paper
// (Adjust `currentDate` between 2025-07-15 and 2025-08-28 to track how
// the algorithm detects this event over time.)
var testClearcut = ee.Geometry.Polygon([
  [29.4686, 50.7278], [29.4685, 50.7267], [29.4698, 50.7267], [29.4698, 50.7265],
  [29.4722, 50.7264], [29.4722, 50.7268], [29.4715, 50.7271], [29.4700, 50.7271],
  [29.4697, 50.7272], [29.4697, 50.7276], [29.4694, 50.7276], [29.4694, 50.7278],
  [29.4686, 50.7278]
]);


// ---------------------------------------------------------------------------
// 3. Algorithm parameters
// ---------------------------------------------------------------------------

// End date of the monitoring window.  Change this to step through the demo.
var currentDate = '2025-08-02';

// Length of the monitoring window (days).  The collection spans
// [currentDate - monitoringWindow, currentDate].
// The paper uses a 60-day (≈ 2-month) window.
var monitoringWindow = 60;

// Number of years to look back when generating the CCDC synthetic reference.
// A 1-year offset provides a phenologically matched pre-disturbance baseline.
var hindcastYears = 1;

// Maximum scene-level cloud cover (%) used to pre-filter the S-2 collection.
var cloudCoverThreshold = 30;

// Number of consecutive positive RdNBR observations required to confirm a
// change (N in the paper).  Higher N → fewer false positives, longer lag.
var nConsecutive = 4;

// RdNBR threshold separating change from no-change pixels.
// Pixels with RdNBR ≥ rdnbrThreshold in N consecutive cloud-free images
// are flagged as confirmed forest disturbances.
var rdnbrThreshold = 655;


// ---------------------------------------------------------------------------
// 4. Build the Sentinel-2 collection for the monitoring window
// ---------------------------------------------------------------------------

var windowStart = ee.Date(currentDate)
                    .advance(ee.Number(monitoringWindow).multiply(-1), 'day');

var s2Collection = buildS2Collection({
  region:     aoi.bounds(),
  start:      windowStart,
  end:        currentDate,
  scale:      20,
  cloudCover: cloudCoverThreshold
});

// Collapse same-day acquisitions (overlapping S-2 swaths) to a single mosaic
// so that each date contributes at most one observation to the time series.
s2Collection = mosaicByDate(s2Collection);

// Diagnostic: number of usable images in the monitoring window
var nImages = s2Collection.size();
print('N images in monitoring window:', nImages);


// ---------------------------------------------------------------------------
// 5. Count clear (cloud-free) pixels per pixel across the window
//    Used as a data-availability diagnostic layer.
// ---------------------------------------------------------------------------

var clearPxCount = s2Collection.map(function (img) {
  return img.select('nir');
}).reduce(ee.Reducer.count());

Map.addLayer(
  clearPxCount,
  {min: 1, max: 10, palette: ['red', 'yellow', 'green']},
  'Clear pixels in monitoring window'
);


// ---------------------------------------------------------------------------
// 6. (Optional) Visualise individual images and the synthetic reference
// ---------------------------------------------------------------------------
// Uncomment this block to inspect all images used in change detection.
// -------------------------------------------------------------------------
// var refDate = ee.Date(currentDate).advance(ee.Number(hindcastYears).multiply(-1), 'year');
// var refSegment = ccdcSegments.findByDate(refDate, 'closest');
// var refImage   = ccdc.getSyntheticImage(ccdcSegments, refDate).clip(aoi);
//
// Map.addLayer(
//   refImage,
//   {bands: 'swir2,nir,red', min: [0, 500, 200], max: [1800, 5000, 3500]},
//   'CCDC synthetic reference ' + currentDate,
//   false
// );
//
// var icList   = s2Collection.toList(s2Collection.size());
// var listSize = icList.size().getInfo();
//
// for (var i = 0; i < listSize; i++) {
//   var img  = ee.Image(icList.get(i));
//   var date = ee.Date(img.get('system:time_start'))
//                .format('YYYY-MM-dd')
//                .getInfo();
//   Map.addLayer(
//     img.clip(aoi),
//     {bands: 'swir2,nir,red', min: [0, 500, 200], max: [1800, 5000, 3500]},
//     'S2 ' + date,
//     false
//   );
// }
// -------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// 7. Compute per-image RdNBR and stack into a 2-D array image
//
//    Array axis 0 = band index within one observation:
//      [0] binary flag   (1 = RdNBR ≥ threshold, 0 = no change, 2 = missing)
//      [1] day-of-year   (integer, 1–365)
//      [2] RdNBR value   (float)
//      [3] clear-pixel count (from Step 5)
//    Array axis 1 = time (one slice per image in the collection)
// ---------------------------------------------------------------------------

var rdnbrStack = ee.ImageCollection(
  s2Collection.map(function (currentImg) {

    // 1-year-prior reference date for CCDC model evaluation
    var imgDate  = ee.Date(currentImg.date());
    var refDate  = imgDate.advance(ee.Number(hindcastYears).multiply(-1), 'year');

    // Day of year for the current observation (1-based)
    var doy      = currentImg.date().getRelative('day', 'year').add(1);
    var doyImage = ee.Image(doy).toInt();

    // Synthetic pre-disturbance reference from CCDC
    var refImage = ccdc.getSyntheticImage(ccdcSegments, refDate);

    // RdNBR = (NBR_ref − NBR_obs) / sqrt(|NBR_ref / 1000|)
    var rdnbrImage = computeRdNBR(refImage.select('NBR'), currentImg.select('NBR'));

    // Binary change flag: 1 = change, 0 = no change
    var changeFlag = rdnbrImage.gte(rdnbrThreshold);

    // Encode missing observations as 2 so the array has a uniform footprint.
    // Pixel values:  0 = no change,  1 = change,  2 = missing observation
    return changeFlag
      .unmask({value: 2, sameFootprint: false})
      .addBands(doyImage)       // axis-0 index 1: day of year
      .addBands(rdnbrImage)     // axis-0 index 2: RdNBR magnitude
      .addBands(clearPxCount);  // axis-0 index 3: available clear pixels
  })
).toArray()           // stack into a 2-D array (axis 0 = band, axis 1 = time)
 .arrayTranspose()    // transpose: axis 0 = time, axis 1 = band
 .clip(aoi);

Map.addLayer(rdnbrStack, {}, 'RdNBR array (monitoring window)', false);


// ---------------------------------------------------------------------------
// 8. Remove missing observations (value = 2) from the per-pixel time series
// ---------------------------------------------------------------------------

// Mask applied along axis 0 (time): keep only valid observations (≠ 2)
var validObsMask  = rdnbrStack.arraySlice(0, 0, 1).neq(2);
var validObsArray = rdnbrStack.arrayMask(validObsMask);


// ---------------------------------------------------------------------------
// 9. Extract the N+1 most recent valid observations per pixel
//    (N+1 because the confirmation rule requires N consecutive positives
//     preceded by one observation to anchor the sequence start.)
// ---------------------------------------------------------------------------

var recentObsArray = validObsArray
  .arraySlice({axis: 1, start: ee.Number(nConsecutive).add(1).multiply(-1)});

Map.addLayer(recentObsArray, {}, 'Recent observations (N+1 window)', false);


// ---------------------------------------------------------------------------
// 10. Past-change mask
//     A pixel is labelled "past change" if all N+1 recent observations are
//     positive (sum of change flags = N+1), indicating it was already flagged
//     before the most recent confirmation window.
// ---------------------------------------------------------------------------

// Sum of change flags across all N+1 recent observations
var sumAllObs    = recentObsArray.arraySlice(0, 0, 1)
                    .arrayReduce({reducer: ee.Reducer.sum(), axes: [1]});

var pastChangeMask = sumAllObs
  .gt(nConsecutive)               // all N+1 positions were positive
  .arrayProject([0])
  .arrayFlatten([['past']])
  .updateMask(fnf)
  .selfMask();


// ---------------------------------------------------------------------------
// 11. Confirmed-change detection
//     A pixel is confirmed if the N most recent observations (excluding the
//     oldest anchor) are all positive:
//       • sumLastN  = N   (the N trailing flags are all 1)
//       • sumAllObs = N   (the anchor flag is 0, ruling out past changes)
//     This corresponds to the pattern  0-1-1-…-1  (N ones preceded by 0).
// ---------------------------------------------------------------------------

// Sum of change flags over the N most recent positions (drop the anchor)
var sumLastN = recentObsArray.arraySlice(0, 0, 1)
                .arraySlice(1, 1)               // skip the oldest anchor
                .arrayReduce({reducer: ee.Reducer.sum(), axes: [1]});

// A pixel is a confirmed new change when:
//   sumLastN  = nConsecutive  AND  sumAllObs = nConsecutive
var confirmedChangeMask = sumLastN.eq(nConsecutive).and(sumAllObs.eq(nConsecutive));

// Detection date: DOY of the first positive observation in the confirmed window
// (axis-0 index 1 = day-of-year band; axis-1 index 1 = second position = first
//  confirmed flag after the anchor)
var confirmedChangeDOY = recentObsArray
  .arraySlice(0, 1, 2)                     // DOY band
  .arrayMask(confirmedChangeMask)
  .arraySlice(1, 1, 2)                     // first confirmed position
  .arrayProject([0])
  .arrayFlatten([['start_DOY']])
  .updateMask(fnf).selfMask();

// Last RdNBR value in the confirmed window (intensity proxy)
var confirmedChangeRdNBR = recentObsArray
  .arraySlice(0, 2, 3)                     // RdNBR band
  .arrayMask(confirmedChangeMask)
  .arraySlice(1, -1)                       // most recent observation
  .arrayProject([0])
  .arrayFlatten([['last_rdnbr']])
  .updateMask(fnf).selfMask();

// Number of cloud-free images available at confirmed-change pixels
var confirmedChangeClearPx = recentObsArray
  .arraySlice(0, 3, 4)                     // clear-pixel-count band
  .arrayMask(confirmedChangeMask)
  .arraySlice(1, -1)
  .arrayProject([0])
  .arrayFlatten([['n_clearPx']])
  .updateMask(fnf).selfMask();


// ---------------------------------------------------------------------------
// 12. Visualisation
// ---------------------------------------------------------------------------

Map.addLayer(
  pastChangeMask,
  {palette: ['cyan']},
  'Past detected changes [within 1-year window]'
);

Map.addLayer(
  confirmedChangeDOY,
  {palette: ['magenta']},
  'Confirmed changes [detection DOY]'
);

Map.addLayer(
  confirmedChangeRdNBR,
  {palette: ['#f0a500', '#ff6f3c', '#d62828', '#800000'],
   min: rdnbrThreshold, max: 2000},
  'Confirmed changes [RdNBR intensity]'
);

// AOI and test clearcut outlines
var empty = ee.Image();
Map.addLayer(
  empty.paint({featureCollection: testClearcut, width: 1}),
  {palette: ['blue']},
  'Test clearcut perimeter'
);
Map.addLayer(
  empty.paint({featureCollection: aoi, width: 1}),
  {palette: ['black']},
  'Area of interest (AOI)'
);

Map.centerObject(aoi, 15);


// =============================================================================
// Helper functions
// =============================================================================

/**
 * Compute the Relativized differenced Normalized Burn Ratio (RdNBR).
 *
 *   RdNBR = (NBR_pre − NBR_post) / sqrt(|NBR_pre / 1000|)
 *
 * Reference:
 *   Miller, J.D. & Thode, A.E. (2007). Quantifying burn severity in a
 *   heterogeneous landscape with a relative version of the delta Normalized
 *   Burn Ratio (dNBR). RSE, 109(1), 66–80.
 *
 * @param {ee.Image} nbrPre   Pre-disturbance NBR (×10 000, single band)
 * @param {ee.Image} nbrPost  Post-disturbance NBR (×10 000, single band)
 * @returns {ee.Image}  RdNBR image (band name: 'rdnbr')
 */
function computeRdNBR(nbrPre, nbrPost) {
  return nbrPre.subtract(nbrPost)
    .divide(nbrPre.divide(1000).abs().sqrt())
    .rename('rdnbr');
}


/**
 * Compute the Normalized Burn Ratio (NBR) for a Sentinel-2 image.
 *   NBR = (NIR − SWIR2) / (NIR + SWIR2), scaled × 10 000
 *
 * @param {ee.Image} image  Image with bands 'nir' and 'swir2'
 * @returns {ee.Image}  Single-band image named 'NBR'
 */
function toNBR(image) {
  return image
    .normalizedDifference(['nir', 'swir2'])
    .multiply(10000.0)
    .select([0], ['NBR']);
}


/**
 * Build a cloud-masked Sentinel-2 SR collection over a given region and date
 * range, and compute NBR for each image.
 *
 * Cloud masking uses the Scene Classification Layer (SCL):
 *   SCL = 2  Dark-area pixels
 *   SCL = 4  Vegetation
 *   SCL = 5  Bare soil / not-vegetated
 *   SCL = 6  Water
 *
 * Only pixels in these classes are retained; cloud (SCL 8–10) and cloud
 * shadow (SCL 3) pixels are masked out.
 *
 * @param {Object} params
 * @param {ee.Geometry} params.region      Spatial filter geometry
 * @param {ee.Date|string} params.start    Start date
 * @param {ee.Date|string} params.end      End date
 * @param {number} params.scale            Spatial scale (m) for clipping
 * @param {number} params.cloudCover       Maximum scene cloud cover (%)
 * @returns {ee.ImageCollection}  Sorted (ascending time) collection with
 *   bands: blue, green, red, nir, swir1, swir2, NBR
 */
function buildS2Collection(params) {
  var filter = ee.Filter.and(
    ee.Filter.bounds(params.region),
    ee.Filter.date(params.start, params.end),
    ee.Filter.metadata('CLOUD_COVERAGE_ASSESSMENT', 'less_than', params.cloudCover)
  );

  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filter(filter)
    .map(function (img) {
      var scl = img.select('SCL');
      // Retain clear-sky land and water pixels; mask cloud and cloud shadow
      var clearPixels = scl.eq(2)   // dark area
        .or(scl.eq(4))              // vegetation
        .or(scl.eq(5))              // bare soil
        .or(scl.eq(6));             // water
      return img
        .updateMask(clearPixels)
        .select(
          ['B2',    'B3',     'B4',  'B8',  'B11',   'B12'],
          ['blue', 'green',  'red', 'nir', 'swir1', 'swir2']
        );
    });

  return s2
    .map(function (image) {
      return image.addBands(toNBR(image)).clip(params.region);
    })
    .sort('system:time_start');
}


/**
 * Collapse same-day acquisitions (from overlapping Sentinel-2 swaths) into
 * a single per-day mosaic (median composite).
 *
 * Each unique day is identified by a "year + day-of-year" string (format
 * 'yD', e.g. '2025183').  The timestamp of the first image on that day is
 * preserved so that downstream date arithmetic remains correct.
 *
 * @param {ee.ImageCollection} collection  Input collection (any bands)
 * @returns {ee.ImageCollection}  One image per calendar day
 */
function mosaicByDate(collection) {
  // Tag each image with its calendar date as a string key
  collection = collection.map(function (img) {
    var dateKey = ee.Date(img.get('system:time_start')).format('yD');
    return img.set('dateKey', dateKey);
  });

  var uniqueDates = collection.aggregate_array('dateKey').distinct();

  return ee.ImageCollection(
    uniqueDates.map(function (dateKey) {
      var sameDay   = collection.filter(ee.Filter.eq('dateKey', dateKey));
      var timestamp = sameDay.first().get('system:time_start');
      return sameDay.median().set('system:time_start', timestamp);
    })
  );
}
