/*
Remote Sensing Applications: Group Glacial Retreat
Part 1.2 (GEE): Glacial Analysis - get Data
therefore we use Landsat Collection 2 Level-2 (LT04/LT05/LE07/LC08, T1+T2), NASADEM (30m), ROI_Fox (Part 1.1)

Output: Glacial_TieSeries (CSV), SnowMask (TIFF), Snowline (TIFF)

Course: Remote Sensing Applications
Group: Glacial Retreat
By Justin Lingg-Laham, Liev Paustian, Lisa Wohlgemuth

Workflow:
1. Per-scene preprocessing For each Landsat scene:
-Cloud/shadow masking using QA_PIXEL (remove clouds + shadows)
-Band scaling (Surface Reflectance) → green and swir1
-Land Surface Temperature (LST) from the thermal band (°C)
-Clip to roi

2. Two time windows per year

A) Glacier area / snow-area (Aug–May)
Goal: “How large is the snowcovered area in a given year?”
therefore:
-Collect images within the Aug–May window
-Create a median composite (one image per year, robust against outliers)
-Compute NDSI
-Define the mask: NDSI > threshold and swir1 < swir1Max
→ results in a binary mask (1 = snow, 0 = not snow)

Calculate area: mask × pixelArea → sum → km² (glacier_km2)
clearFraction: fraction of valid (unmasked) pixels within the ROI
→ Output per year: glacier_km2, LST_median_C, imgCount, clearFraction, and quality flags.

B) Snowline (Jan–Mar, summer window)
Goal: “Where is the end-of-summer boundary between snow and melted-out terrain?”
therefore:
-Collect only images within Jan–Mar (summer)
-Median composite
-NDSI → summer snow/ice mask (gMaskSnow)
-Quality filter (“doSnowline”): only if Snow ≥ 3
-clearFraction_snow ≥ 0.20 (Threshold)
-glacier_km2 ≥ 5 km² (to avoid tiny patches)
-Remove small patches (connectedPixelCount) → keep only large contiguous areas
-Compute the mask edge
→ this is the “boundary line” between snow and non-snow

then Overlay the DEM and sample elevation along the edge
→ snowline_edge_m = median elevation of edge pixels

also: snowElev_p10/p50/p90_m = elevation distribution of snow pixels 
→ Output per year: a snowline elevation 
*/

// Define ROI for available Data
var bufferAnalysis_m = 300;
var bufferSearch_m   = 8000; // Search for Available Images
var simplify_m       = 20;

var roiBase = ee.FeatureCollection(fox).geometry();
var roiAnalysis = roiBase.buffer(bufferAnalysis_m);
if (simplify_m > 0) roiAnalysis = roiAnalysis.simplify(simplify_m);
var roiSearch = roiBase.buffer(bufferSearch_m);

Map.centerObject(roiBase, 10);

// load DEM
var dem = ee.Image("NASA/NASADEM_HGT/001").select("elevation").clip(roiAnalysis);

// Define Seasons (Glacial Analysis)
var seasonStartMonth = 8, seasonStartDay = 1;   // 01 Aug
var seasonEndMonth   = 5, seasonEndDay   = 31;  // 31 May

// Snowline summer season
var snowStartMonth = 1, snowStartDay = 1;   // 01 Jan
var snowEndMonth   = 3, snowEndDay   = 31;  // 31 Mar

// Parameters based on Try and Error
var ndsiThr   = 0.40; 
var swir1Max  = 0.20; 

var minScenes = 2;               // for the season quality flag
var minClearFraction = 0.10;     // for the season quality flag
var scaleArea = 30;              // 30 m resolution

// Snowline gating (strict)
var minScenes_snowline = 3;
var minClearFraction_snowline = 0.20;
var minGlacierKm2_forSnowline = 5;

// Cloudmask (Wright et al. 2024)
function maskC02L2(img) {
  var qa = img.select("QA_PIXEL");
  var cloud  = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow = qa.bitwiseAnd(1 << 4).neq(0);
  return img.updateMask(cloud.or(shadow).not());
}
function clipAnalysis(img){ return img.clip(roiAnalysis); }

// Landsat 8 processing
function prepL8L2(img){
  var green = img.select("SR_B3").multiply(0.0000275).add(-0.2).rename("green");
  var swir1 = img.select("SR_B6").multiply(0.0000275).add(-0.2).rename("swir1");
  var lstC  = img.select("ST_B10").multiply(0.00341802).add(149.0).subtract(273.15).rename("LST");
  return ee.Image.cat([green, swir1, lstC]).copyProperties(img, img.propertyNames());
}
// Landsat 4-7 processing
function prepL457L2(img){
  var green = img.select("SR_B2").multiply(0.0000275).add(-0.2).rename("green");
  var swir1 = img.select("SR_B5").multiply(0.0000275).add(-0.2).rename("swir1");
  var lstC  = img.select("ST_B6").multiply(0.00341802).add(149.0).subtract(273.15).rename("LST");
  return ee.Image.cat([green, swir1, lstC]).copyProperties(img, img.propertyNames());
}

function colT1T2(idT1, idT2, prepFn){
  var t1 = ee.ImageCollection(idT1).filterBounds(roiSearch).map(maskC02L2).map(prepFn).map(clipAnalysis);
  var t2 = ee.ImageCollection(idT2).filterBounds(roiSearch).map(maskC02L2).map(prepFn).map(clipAnalysis);
  return t1.merge(t2);
}

// Load Landsat Data
var lt04 = colT1T2("LANDSAT/LT04/C02/T1_L2", "LANDSAT/LT04/C02/T2_L2", prepL457L2);
var lt05 = colT1T2("LANDSAT/LT05/C02/T1_L2", "LANDSAT/LT05/C02/T2_L2", prepL457L2);
var le07 = colT1T2("LANDSAT/LE07/C02/T1_L2", "LANDSAT/LE07/C02/T2_L2", prepL457L2);
var lc08 = colT1T2("LANDSAT/LC08/C02/T1_L2", "LANDSAT/LC08/C02/T2_L2", prepL8L2);

// merge Landsat together
var landsat = lt04.merge(lt05).merge(le07).merge(lc08);

// startyear
var firstDate = ee.Date(landsat.sort("system:time_start").first().get("system:time_start"));
var startYear = ee.Number.parse(firstDate.format("YYYY")); 
var endYear   = 2025;

print("First available date:", firstDate);
print("Start year used:", startYear);

var years = ee.List.sequence(startYear, endYear);

// Glacial Statistics (Aug-May)
function areaKm2(binaryMask){
  var m2 = ee.Number(binaryMask.rename("m")
    .multiply(ee.Image.pixelArea())
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: roiAnalysis,
      scale: scaleArea,
      tileScale: 8,
      bestEffort: true
    }).get("m"));
  return m2.divide(1000000);
}

function clearFraction(comp){
  var valid = comp.select("green").mask();
  var frac = ee.Number(valid.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roiAnalysis,
    scale: 60,
    tileScale: 8,
    bestEffort: true
  }).get("green"));
  return frac;
}

function fmtDate(ms){ return ee.Date(ms).format("YYYY-MM-dd"); }

// Statistics per Year
function statsForYear(y){
  y = ee.Number(y);
  //  Define Season (Aug–May) 
  var startS = ee.Date.fromYMD(y, seasonStartMonth, seasonStartDay);
  var endS   = ee.Date.fromYMD(y.add(1), seasonEndMonth, seasonEndDay).advance(1,"day");

  var col = landsat.filterDate(startS, endS);
  var n = col.size();
  var windowUsed = "SEASON";

  var first = ee.Algorithms.If(n.gt(0), fmtDate(col.sort("system:time_start").first().get("system:time_start")), null);
  var last  = ee.Algorithms.If(n.gt(0), fmtDate(col.sort("system:time_start", false).first().get("system:time_start")), null);

  var comp = ee.Image(ee.Algorithms.If(n.gt(0), col.median(), ee.Image(0)));

  // clear fraction 
  var cf_num = ee.Number(ee.Algorithms.If(n.gt(0), clearFraction(comp), 0));
  var cf = ee.Algorithms.If(n.gt(0), cf_num, null);

  var ndsi = ee.Image(ee.Algorithms.If(
    n.gt(0),
    comp.normalizedDifference(["green","swir1"]).rename("NDSI"),
    ee.Image(0).rename("NDSI")
  ));

  var gMask = ndsi.gt(ndsiThr).and(comp.select("swir1").lt(swir1Max)).rename("m");

  // glacier area 
  var glacier_km2_num = ee.Number(ee.Algorithms.If(n.gt(0), areaKm2(gMask), 0));
  var glacier_km2 = ee.Algorithms.If(n.gt(0), glacier_km2_num, null);

  var lstMed = ee.Algorithms.If(
    n.gt(0),
    comp.select("LST").reduceRegion({
      reducer: ee.Reducer.median(),
      geometry: roiAnalysis,
      scale: 120,
      tileScale: 8,
      bestEffort: true
    }).get("LST"),
    null
  );

  // Season quality flaggs (for better understanding in R)
  var okScenes  = n.gte(minScenes);
  var okClear   = ee.Algorithms.If(n.gt(0), cf_num.gte(minClearFraction), false);
  var okOverall = ee.Algorithms.If(n.gt(0), n.gte(minScenes).and(cf_num.gte(minClearFraction)), false);

  //  Define Snowline summer (Jan–Mar)
  var snowStart = ee.Date.fromYMD(y, snowStartMonth, snowStartDay);
  var snowEnd   = ee.Date.fromYMD(y, snowEndMonth, snowEndDay).advance(1,"day");

  var colSnow = landsat.filterDate(snowStart, snowEnd);
  var nSnow = colSnow.size();

  var compSnow = ee.Image(ee.Algorithms.If(nSnow.gt(0), colSnow.median(), ee.Image(0)));

  // clearFraction_snow robust
  var cfSnow_num = ee.Number(ee.Algorithms.If(nSnow.gt(0), clearFraction(compSnow), 0));
  var cfSnow = ee.Algorithms.If(nSnow.gt(0), cfSnow_num, null);

  var ndsiSnow = ee.Image(ee.Algorithms.If(
    nSnow.gt(0),
    compSnow.normalizedDifference(["green","swir1"]).rename("NDSI"),
    ee.Image(0).rename("NDSI")
  ));

  var gMaskSnow = ndsiSnow.gt(ndsiThr).and(compSnow.select("swir1").lt(swir1Max)).rename("m");

  // clean snowline (using Gating)
  var hasBigIce = glacier_km2_num.gte(minGlacierKm2_forSnowline);
  var hasGoodSnow = nSnow.gte(minScenes_snowline).and(cfSnow_num.gte(minClearFraction_snowline));
  var doSnowline = hasGoodSnow.and(hasBigIce);

  // largest-blob cleaning on summer mask
  var g0 = gMaskSnow.selfMask();
  var patchSize = g0.connectedPixelCount(1024, true);
  var gClean = g0.updateMask(patchSize.gte(200)); // 200 px ca. 0.18 km2 = 30m
  
  var gC0 = gClean.unmask(0);
  var edge = gC0.focal_max(2).neq(gC0.focal_min(2)).selfMask();

  var snowline_edge_m = ee.Algorithms.If(
    doSnowline,
    dem.updateMask(edge).reduceRegion({
      reducer: ee.Reducer.median(),
      geometry: roiAnalysis,
      scale: 30,
      tileScale: 8,
      bestEffort: true
    }).get("elevation"),
    null
  );

  var elevDict = ee.Dictionary(ee.Algorithms.If(
    doSnowline,
    dem.updateMask(gClean).reduceRegion({
      reducer: ee.Reducer.percentile([10, 50, 90]),
      geometry: roiAnalysis,
      scale: 30,
      tileScale: 8,
      bestEffort: true
    }),
    ee.Dictionary({})
  ));

// Define Elevation percentiles
  var elev_p10 = ee.Algorithms.If(elevDict.contains("elevation_p10"), elevDict.get("elevation_p10"), null);
  var elev_p50 = ee.Algorithms.If(elevDict.contains("elevation_p50"), elevDict.get("elevation_p50"), null);
  var elev_p90 = ee.Algorithms.If(elevDict.contains("elevation_p90"), elevDict.get("elevation_p90"), null);

//Output per year: imgCount, glacier_km2, LST_median_C, clearFraction + snowline_edge_m, snowElev_p10/p50/p90, nSnow, clearFraction_snow, doSnowline
  return ee.Feature(null, {
    year: y,
    // Glacial
    windowUsed: windowUsed,
    imgCount: n,
    firstDate: first,
    lastDate: last,
    clearFraction: cf,

    glacier_km2: glacier_km2,
    LST_median_C: lstMed,

    ok_scenes: okScenes,
    ok_clear: okClear,
    ok_overall: okOverall,

    // snowline 
    nSnow: nSnow,
    clearFraction_snow: cfSnow,
    doSnowline: doSnowline,

    snowline_edge_m: snowline_edge_m,
    snowElev_p10_m: elev_p10,
    snowElev_p50_m: elev_p50,
    snowElev_p90_m: elev_p90,

    // parameters
    season: seasonStartMonth + "/" + seasonStartDay + "–" + seasonEndMonth + "/" + seasonEndDay,
    snowWindow: snowStartMonth + "/" + snowStartDay + "–" + snowEndMonth + "/" + snowEndDay,
    ndsiThr: ndsiThr,
    swir1Max: swir1Max,
    minScenes_quality: minScenes,
    minClearFraction: minClearFraction,
    minScenes_snowline: minScenes_snowline,
    minClearFraction_snowline: minClearFraction_snowline,
    minGlacierKm2_forSnowline: minGlacierKm2_forSnowline,
    scale_m: scaleArea,
    dem: "NASADEM_HGT_001"
  });
}

var statsALL = ee.FeatureCollection(years.map(statsForYear));

function seasonMaskImageForYear(y){
  y = ee.Number(y);

  var startS = ee.Date.fromYMD(y, seasonStartMonth, seasonStartDay);
  var endS   = ee.Date.fromYMD(y.add(1), seasonEndMonth, seasonEndDay).advance(1,"day");

  var col = landsat.filterDate(startS, endS);
  var n = col.size();

  // median composite 
  var comp = ee.Image(ee.Algorithms.If(n.gt(0), col.median(), ee.Image(0)));

  var ndsi = comp.normalizedDifference(["green","swir1"]).rename("NDSI");
  var gMask = ndsi.gt(ndsiThr).and(comp.select("swir1").lt(swir1Max)).rename("snowice");

  // make Binary Output
  var out = gMask.unmask(0).toByte()
    .set({
      year: y,
      imgCount: n,
      start: startS.format("YYYY-MM-dd"),
      end: endS.format("YYYY-MM-dd")
    });

  return out;
}
// Export every Stats as CSV
Export.table.toDrive({
  collection: statsALL,
  description: "FoxGlacier_LandsatC02L2_PHYS_CLEAN_" + startYear + "_" + endYear,
  folder: "GEE_Exports",
  fileNamePrefix: "FoxGlacier_LandsatC02L2_PHYS_CLEAN_" + startYear + "_" + endYear,
  fileFormat: "CSV"
});
var exportYears = years; 

exportYears.getInfo().forEach(function(y){
  var img = seasonMaskImageForYear(y);

  Export.image.toDrive({
    image: img,
    description: "FOX_snowiceMask_SEASON_" + y,
    folder: "GEE_Exports",
    fileNamePrefix: "FOX_snowiceMask_SEASON_" + y,
    region: roiAnalysis,
    scale: 30,
    crs: "EPSG:2193",
  });
});

// Snowline function to extract for every Year
function snowlineEdgeImageForYear(y){
  y = ee.Number(y);

  var snowStart = ee.Date.fromYMD(y, snowStartMonth, snowStartDay);
  var snowEnd   = ee.Date.fromYMD(y, snowEndMonth, snowEndDay).advance(1,"day");

  var colSnow = landsat.filterDate(snowStart, snowEnd);
  var nSnow = colSnow.size();

  var compSnow = ee.Image(ee.Algorithms.If(nSnow.gt(0), colSnow.median(), ee.Image(0)));
  var ndsiSnow = compSnow.normalizedDifference(["green","swir1"]).rename("NDSI");
  var gMaskSnow = ndsiSnow.gt(ndsiThr).and(compSnow.select("swir1").lt(swir1Max)).rename("m");

  // cleaning
  var g0 = gMaskSnow.selfMask();
  var patchSize = g0.connectedPixelCount(1024, true);
  var gClean = g0.updateMask(patchSize.gte(200)).rename("snow");

  var gC0 = gClean.unmask(0);
  var edge = gC0.focal_max(2).neq(gC0.focal_min(2)).rename("edge").selfMask();

  // make it binary
  var snow01 = gClean.unmask(0).toByte().rename("snow");
  var edge01 = edge.unmask(0).toByte().rename("edge");

  return ee.Image.cat([snow01, edge01]).set({year: y, nSnow: nSnow});
}
exportYears.getInfo().forEach(function(y){
  var img = snowlineEdgeImageForYear(y);

  Export.image.toDrive({
    image: img,
    description: "FOX_snowline_snow_edge_" + y,
    folder: "GEE_Exports",
    fileNamePrefix: "FOX_snowline_snow_edge_" + y,
    region: roiAnalysis,
    scale: 30,
    crs: "EPSG:2193",
  });
});


