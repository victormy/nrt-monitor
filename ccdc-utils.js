// =============================================================================
// ccdc-utils.js
// Minimal CCDC segment utilities for NRT forest change detection
//
// Paper:   Operational event-based forest change detection using
//          Sentinel-2 time series
// Author:  Viktor Myroniuk
// Version: 2026-09-01
//
// Attribution:
//   The harmonic evaluation and segment-selection logic in this file is
//   derived from the temporalSegmentation library by Daniel Wiell
//   (https://code.earthengine.google.com/?accept_repo=users/wiell/temporalSegmentation)
//   and has been simplified to the subset required by the NRT monitoring
//   workflow.  Only the functions actually used in nrt-monitoring-demo.js
//   are retained; all other functionality has been removed.
//
// Usage (GEE Code Editor):
//   var ccdc = require('users/<username>/<repo>:ccdc-utils');
//   var segments  = ccdc.Segments(ccdcImage);
//   var synImage  = ccdc.getSyntheticImage(segments, date);
// =============================================================================


// ---------------------------------------------------------------------------
// Date-format constants (matches GEE CCDC dateFormat parameter)
// ---------------------------------------------------------------------------
var J_DAYS          = 0;   // Julian days (default for GEE CCDC exports)
var FRACTIONAL_YEARS = 1;
var UNIX_TIME_MILLIS = 2;


// ---------------------------------------------------------------------------
// dateConversion — convert between ee.Date and the numeric t used in CCDC
// ---------------------------------------------------------------------------
var dateConversion = {

  /**
   * Convert an ee.Date (or JS Date string) to the numeric t value used
   * by CCDC coefficients.
   *
   * @param {ee.Date|string} date
   * @param {number} dateFormat  0 = Julian days, 1 = fractional years,
   *                             2 = Unix milliseconds
   * @returns {ee.Number|number}
   */
  toT: function (date, dateFormat) {
    if (date instanceof ee.Date) {
      date = ee.Date(date);
      switch (dateFormat) {
        case J_DAYS:
          var epochDay = 719529;
          return date.millis().divide(1000).divide(3600).divide(24).add(epochDay);
        case FRACTIONAL_YEARS:
          return date.get('year').add(date.getFraction('year'));
        case UNIX_TIME_MILLIS:
          return date.millis();
        default:
          throw Error('Unsupported dateFormat: ' + dateFormat +
            '. Use 0 (Julian days), 1 (fractional years), or 2 (Unix ms).');
      }
    } else {
      // JavaScript Date path (used by chartPoint — not needed for NRT demo
      // but kept for completeness)
      date = new Date(date);
      switch (dateFormat) {
        case 0: {
          var epochDay = 719529;
          return date.getTime() / 1000 / 3600 / 24 + epochDay;
        }
        case 1: {
          var firstOfYear     = new Date(Date.UTC(date.getFullYear(), 0, 1, 0, 0, 0));
          var firstOfNextYear = new Date(Date.UTC(date.getFullYear() + 1, 0, 1, 0, 0, 0));
          var fraction = (date - firstOfYear) / (firstOfNextYear - firstOfYear);
          return date.getFullYear() + fraction;
        }
        case 2:
          return date.getTime();
        default:
          throw Error('Unsupported dateFormat: ' + dateFormat);
      }
    }
  },

  /**
   * Compute the absolute difference in days between two t-images, given the
   * active dateFormat.
   */
  days: function (t1, t2, dateFormat) {
    var diff = t2.subtract(t1);
    switch (dateFormat) {
      case J_DAYS:          return diff;
      case FRACTIONAL_YEARS: return diff.multiply(365).round();
      case UNIX_TIME_MILLIS: return diff.divide(1000 * 3600 * 24).round();
      default:
        throw Error('Unsupported dateFormat: ' + dateFormat);
    }
  }
};


// ---------------------------------------------------------------------------
// harmonicSlice — evaluate the CCDC harmonic model at a given t
// ---------------------------------------------------------------------------
var harmonicSlice = (function () {

  /**
   * Angular frequency (omega) for the harmonic terms, keyed by dateFormat.
   */
  function getOmega(dateFormat) {
    switch (dateFormat) {
      case 0: return 2.0 * Math.PI / 365.25;           // Julian days
      case 1: return 2.0 * Math.PI;                    // fractional years
      case 2: return 2.0 * Math.PI / (60*60*24*365.25); // Unix ms
      default: throw Error('Unsupported dateFormat: ' + dateFormat);
    }
  }

  /**
   * Internal: build a GEE List of per-harmonic contribution images/numbers.
   *
   * @param {ee.Image|ee.Array} coefs  Band coefficients (8 values per band)
   * @param {ee.Image|ee.Number} t     Numeric time value
   * @param {number} dateFormat
   * @param {number} harmonics         Number of harmonic pairs to include (0–3)
   * @param {function} coefExtractor   Function(index) → scalar or image
   * @returns {ee.List}
   */
  function buildTerms(coefs, t, dateFormat, harmonics, coefExtractor) {
    dateFormat = (dateFormat === undefined) ? 0 : dateFormat;
    harmonics  = (harmonics  === undefined) ? 3 : harmonics;
    var w = getOmega(dateFormat);

    function c(i) { return coefExtractor(i); }

    return ee.List([
        // Intercept + linear trend
        c(0).add(c(1).multiply(t)),
        // 1st harmonic pair
        c(2).multiply(t.multiply(w    ).cos())
            .add(c(3).multiply(t.multiply(w    ).sin())),
        // 2nd harmonic pair
        c(4).multiply(t.multiply(w * 2).cos())
            .add(c(5).multiply(t.multiply(w * 2).sin())),
        // 3rd harmonic pair
        c(6).multiply(t.multiply(w * 3).cos())
            .add(c(7).multiply(t.multiply(w * 3).sin()))
      ])
      .slice(0, ee.Number(harmonics).add(1));
  }

  /**
   * Evaluate CCDC harmonic coefficients at time t and return a multi-band
   * ee.Image (one value per spectral band).
   *
   * @param {ee.Image} coefs    Array image — bands ending in '_coefs'
   * @param {ee.Image} t        Scalar image with the target t value
   * @param {number}   dateFormat
   * @param {number}   harmonics
   * @returns {ee.Image}
   */
  function sliceImage(coefs, t, dateFormat, harmonics) {
    return ee.ImageCollection(
        buildTerms(coefs, t, dateFormat, harmonics, function (index) {
          return coefs.arrayGet([index]);
        })
      )
      .reduce(ee.Reducer.sum())
      .regexpRename('(.*)_coefs_sum', '$1', false);
  }

  return {
    sliceImage: sliceImage,
    getOmega:   getOmega
  };
}());


// ---------------------------------------------------------------------------
// Internal helpers for Segments / Segment factories
// ---------------------------------------------------------------------------

/** Merge an array of plain objects into one (last writer wins). */
function mergeObjects(objects) {
  objects = objects || {};
  return objects.reduce(function (acc, o) {
    for (var a in o) { acc[a] = o[a]; }
    return acc;
  }, {});
}

/** JavaScript range helper: returns [start, start+1, …, end]. */
function sequence(start, end) {
  return Array.apply(start, Array(end)).map(function (_, i) { return i + start; });
}

/**
 * Ensure that 1-D and 2-D array bands are unmasked and that the overall
 * image mask reflects whether any segment exists at a pixel.
 */
function updateImageMask(segmentsImage) {
  var bands1D = segmentsImage.bandNames()
    .filter(ee.Filter.stringEndsWith('item', '_coefs').not());
  var bands2D = segmentsImage.bandNames()
    .filter(ee.Filter.stringEndsWith('item', '_coefs'));

  return segmentsImage
    .addBands(
      segmentsImage.select(bands1D)
        .unmask(ee.Array([], ee.PixelType.double())),
      null, true
    )
    .addBands(
      segmentsImage.select(bands2D)
        .unmask(ee.Array([[]], ee.PixelType.double())),
      null, true
    )
    .mask(segmentsImage.select(0).arrayLength(0).unmask(0));
}

/** Return a 0-based array of segment indexes for each pixel. */
function getSegmentIndexes(segmentsImage) {
  return segmentsImage
    .select(0).not().not()           // 1 per segment
    .arrayAccum(0, ee.Reducer.sum()) // cumulative count
    .subtract(1);                    // 0-based
}

/** Extract a single segment's bands as a flat (non-array) image. */
function getSegmentImage(segmentsImage, segmentIndex) {
  var bands1D = segmentsImage.bandNames()
    .filter(ee.Filter.stringEndsWith('item', '_coefs').not());
  var bands2D = segmentsImage.bandNames()
    .filter(ee.Filter.stringEndsWith('item', '_coefs'));

  var mask = getSegmentIndexes(segmentsImage).eq(segmentIndex.unmask(-1));

  var image1D = segmentsImage.select(bands1D).arrayMask(mask);
  image1D = image1D
    .mask(image1D.select(0).arrayLength(0).unmask(0))
    .arrayProject([0])
    .arrayGet([0]);

  var image2D = segmentsImage.select(bands2D)
    .arrayMask(mask.toArray(1).unmask(ee.Array([[]], ee.PixelType.double())));
  image2D = image2D
    .mask(image2D.select(0).arrayLength(0).unmask(0))
    .arrayProject([1]);

  return image1D.addBands(image2D);
}


// ---------------------------------------------------------------------------
// Segment — evaluate one CCDC segment at a given date
// ---------------------------------------------------------------------------

/**
 * Wrap a single-segment image (flat bands, not array) so that it can be
 * evaluated at an arbitrary date via slice().
 *
 * @param {ee.Image} segmentImage  Flat image with tStart, tEnd, *_coefs, …
 * @param {number}   dateFormat
 * @param {ee.Date}  [defaultDate] Optional: date to evaluate at by default
 * @returns {Object}  { slice, toImage }
 */
function Segment(segmentImage, dateFormat, defaultDate) {
  var defaultT = defaultDate
    ? ee.Image(dateConversion.toT(defaultDate, dateFormat))
    : segmentImage.expression('(i.tStart + i.tEnd) / 2', {i: segmentImage});

  /**
   * Evaluate the harmonic model at a given date and return a synthetic image.
   *
   * @param {Object} options
   * @param {number}  [options.harmonics=3]          Number of harmonic pairs
   * @param {string}  [options.strategy='mask']       'mask' | 'closest'
   * @param {number}  [options.extrapolateMaxDays=0]  Allow extrapolation (days)
   * @returns {ee.Image}  Synthetic reflectance image (one band per spectral band)
   */
  function slice(options) {
    var defaults = {
      t:                 defaultT,
      harmonics:         3,
      extrapolateMaxDays: 0,
      extrapolateMaxFraction: 0,
      strategy:          'mask'
    };
    options = mergeObjects([defaults, options]);

    var t      = ee.Image(options.date ? dateConversion.toT(options.date, dateFormat) : options.t);
    var tStart = segmentImage.select('tStart');
    var tEnd   = segmentImage.select('tEnd');
    var coefs  = segmentImage.select('.*_coefs');

    var daysFromStart = dateConversion.days(t, tStart, dateFormat);
    var daysFromEnd   = dateConversion.days(tEnd, t, dateFormat);

    if (options.strategy !== 'closest') {
      // Optionally allow extrapolation beyond segment bounds
      var extrapolateMaxDays = options.extrapolateMaxFraction
        ? dateConversion.days(tStart, tEnd, dateFormat)
            .multiply(options.extrapolateMaxFraction).round()
        : ee.Image(options.extrapolateMaxDays);
      extrapolateMaxDays = extrapolateMaxDays
        .where(extrapolateMaxDays.lt(0), ee.Image(Number.MAX_SAFE_INTEGER));

      var daysFromSegment = daysFromStart.max(daysFromEnd).max(0);

      return harmonicSlice.sliceImage(coefs, t, dateFormat, options.harmonics)
        .updateMask(extrapolateMaxDays.gte(daysFromSegment));

    } else {
      // 'closest' strategy: clamp t to segment boundaries instead of masking
      var tUsed = t
        .where(daysFromStart.gt(0), tStart)
        .where(daysFromEnd.gt(0),   tEnd);
      return harmonicSlice.sliceImage(coefs, tUsed, dateFormat, options.harmonics);
    }
  }

  /** Return the underlying flat segment image (or a band selection of it). */
  function toImage(selector) {
    return selector === undefined
      ? segmentImage
      : segmentImage.select(selector);
  }

  return {
    slice:   slice,
    toImage: toImage
  };
}


// ---------------------------------------------------------------------------
// Segments — wrap a full CCDC array image and locate segments by date
// ---------------------------------------------------------------------------

/**
 * Wrap a CCDC multi-segment array image so that individual segments can be
 * retrieved by date.
 *
 * @param {ee.Image} segmentsImage  CCDC output image (array bands)
 * @param {number}   [dateFormat=0] Date encoding used in the CCDC asset
 * @param {number}   [maxSegments=50]
 * @returns {Object}  { findByDate, toImage }
 */
function Segments(segmentsImage, dateFormat, maxSegments) {
  segmentsImage = updateImageMask(segmentsImage);
  // Pad coefficient arrays to 8 elements so arithmetic is always safe
  segmentsImage = segmentsImage
    .addBands(segmentsImage.select('.*_coefs').arrayPad([0, 8]), null, true);
  dateFormat  = (dateFormat  === undefined) ? 0  : dateFormat;
  maxSegments = maxSegments ? maxSegments : 50;

  /**
   * Compute a per-pixel segment index that best matches the given date under
   * the requested strategy.
   *
   * @param {ee.Date} date
   * @param {string}  [strategy='mask']  'mask' | 'closest' | 'previous' | 'next'
   * @returns {ee.Image}  int8 image — pixel value = segment index (0-based)
   */
  function segmentIndex(date, strategy) {
    strategy = strategy || 'mask';
    var t              = ee.Image(dateConversion.toT(date, dateFormat));
    var segmentIndexes = getSegmentIndexes(segmentsImage);
    var tStart         = segmentsImage.select('tStart');
    var tEnd           = segmentsImage.select('tEnd');
    var masked;

    function getPrevious() {
      return segmentIndexes
        .arrayMask(tStart.lte(t))
        .arrayReduce(ee.Reducer.lastNonNull(), [0]);
    }

    function getNext() {
      return segmentIndexes
        .arrayMask(tEnd.gt(t))
        .arrayReduce(ee.Reducer.first(), [0]);
    }

    if (strategy === 'mask') {
      // Only return a valid index if the date falls within a segment
      masked = segmentIndexes
        .arrayMask(tStart.lte(t).and(tEnd.gte(t)))
        .arrayReduce(ee.Reducer.first(), [0]);

    } else if (strategy === 'closest') {
      // Pick the segment whose boundary is nearest to the requested date
      var prevDist = tStart.subtract(t).abs()
        .arrayReduce(ee.Reducer.min(), [0]).arrayGet([0]);
      var nextDist = tEnd.subtract(t).abs()
        .arrayReduce(ee.Reducer.min(), [0]).arrayGet([0]);
      masked = getPrevious().where(nextDist.gt(prevDist), getNext());

    } else if (strategy === 'previous') {
      masked = getPrevious();

    } else if (strategy === 'next') {
      masked = getNext();

    } else {
      throw new Error('Unsupported strategy: "' + strategy + '". ' +
        'Supported: mask (default), closest, previous, next.');
    }

    return masked
      .updateMask(masked.arrayLength(0).gt(0))
      .arrayFlatten([['segmentIndex']])
      .int8();
  }

  /**
   * Retrieve the CCDC segment that best matches the given date.
   *
   * @param {ee.Date} date
   * @param {string}  [strategy='mask']  See segmentIndex() for options
   * @returns {Segment}
   */
  function findByDate(date, strategy) {
    return Segment(
      getSegmentImage(segmentsImage, segmentIndex(date, strategy)),
      dateFormat,
      date
    );
  }

  /** Return the raw array image (or a band selection). */
  function toImage(selector) {
    return selector === undefined
      ? segmentsImage
      : segmentsImage.select(selector);
  }

  return {
    findByDate: findByDate,
    toImage:    toImage
  };
}


// ---------------------------------------------------------------------------
// getSyntheticImage — high-level helper used by the NRT monitoring script
// ---------------------------------------------------------------------------

/**
 * Generate a phenologically-matched synthetic reference image from a CCDC
 * segment at the requested date.
 *
 * Internally this calls Segment.slice() with strategy='closest' and a
 * 30-day extrapolation allowance so that dates just outside a segment's
 * bounds still return a valid image rather than a masked result.
 *
 * @param {Object}  segments  A Segments() wrapper object
 * @param {ee.Date} date      Target date for the synthetic image
 * @returns {ee.Image}  Synthetic reflectance image (bands: blue, green, red,
 *                      nir, swir1, swir2, NBR, …)
 */
function getSyntheticImage(segments, date) {
  var segment = segments.findByDate(date, 'closest');
  return ee.Image(
    segment.slice({
      harmonics:         3,
      strategy:          'closest',
      extrapolateMaxDays: 30
    })
  );
}


// ---------------------------------------------------------------------------
// Module exports
// ---------------------------------------------------------------------------
exports = {
  Segments:          Segments,
  Segment:           Segment,
  getSyntheticImage: getSyntheticImage,
  dateConversion:    dateConversion,
  harmonicSlice:     harmonicSlice
};
