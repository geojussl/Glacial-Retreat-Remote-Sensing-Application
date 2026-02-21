/*
Remote Sensing Applications: Group Glacial Retreat
Part 1.4 (GEE): Hypsometry
Therefore we use Landsat Collection 2 Level-2 (LT04/LT05/LE07/LC08, T1+T2), NASADEM (30m), ROI_FOX

Outputs:
- Hypsometry bands/classes (CSV), valid Years (Hypsometry after quality filter)  

Course: Remote Sensing Applications
Group: Glacial Retreat
By Justin Lingg-Laham, Liev Paustian, Lisa Wohlgemuth

Workflow:

1. Define Elevation-Bands (Classes)= called elev_bin
- Elevation bands: binSize (50 m) → elev_bin = floor(elevation / binSize) * binSize

2. Per-scene preprocessing For each Landsat scene (Same workflow as in Part 1.2 Analysis):
-Cloud/shadow masking using QA_PIXEL (remove clouds + shadows)
-Band scaling (Surface Reflectance) → green and swir1
-Land Surface Temperature (LST) from the thermal band (°C)
-Clip to roi

3. Calculate Hypsometry per year (area per elevation band)
For each valid year:
- Combine: elevation bands (elev_bin) + glacier mask
- Compute area per pixel (km²) and keep only pixels within the mask
- reduceRegion with group reducer:
  - sum(area_km2) grouped by elev_bin
→ Output per year: a table of (elev_bin, area_km2) for that year

4. Period aggregation (Early vs Late) = (pre2006 vs post2006)
Define two periods:
- Early: 1989–2005
- Late:  2006–2022
For each period:
- Collect all yearly hypsometry rows
- Compute mean area per elevation bin across years
- Track n_years_in_bin (how many years contributed to that bin)
→ Output: hypsometry curves for Early and Late, comparable by elevation bin

Outputs (files):
- Fox_hypsometry_bins.csv 
contains: elev_bin (m), area_km2 (mean snowarea in that band), year_group, n_years_in_bin
- Fox_hypsometry_validYears.csv
*/

// Define ROI
var bufferAnalysis_m = 300;
var simplify_m = 20;

var roiBase = ee.FeatureCollection(fox).geometry();
var roi = roiBase.buffer(bufferAnalysis_m);
if (simplify_m > 0) roi = roi.simplify(simplify_m);

Map.centerObject(roiBase, 10);

// Use DEM
var dem = ee.Image("NASA/NASADEM_HGT/001").select("elevation").clip(roiAnalysis);

// Define Groups (Have to be equal)
var groupA = {name: "Early", y0: 1989, y1: 2005}; //pre2006
var groupB = {name: "Late",  y0: 2006, y1: 2022}; //post2006

// Define Elevation Bands/Classes (bins)
var binSize = 50; // a band ebvery 50 m
var elevBin = dem.divide(binSize).floor().toInt().rename("elev_bin");

// Define Seasons
var seasonStartMonth = 8, seasonStartDay = 1;
var seasonEndMonth   = 5, seasonEndDay   = 31;

//Same Instructions as in Part 1.2 (Landsat processing)
var ndsiThr  = 0.40;
var swir1Max = 0.20;

var minScenes_hypso = 4;
var minClearFrac_hypso = 0.20;

function maskC02L2(img) {
  var qa = img.select("QA_PIXEL");
  var cloud  = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow = qa.bitwiseAnd(1 << 4).neq(0);
  return img.updateMask(cloud.or(shadow).not());
}

function prepL8L2(img){
  var green = img.select("SR_B3").multiply(0.0000275).add(-0.2).rename("green");
  var swir1 = img.select("SR_B6").multiply(0.0000275).add(-0.2).rename("swir1");
  return ee.Image.cat([green, swir1]).copyProperties(img, img.propertyNames());
}
function prepL457L2(img){
  var green = img.select("SR_B2").multiply(0.0000275).add(-0.2).rename("green");
  var swir1 = img.select("SR_B5").multiply(0.0000275).add(-0.2).rename("swir1");
  return ee.Image.cat([green, swir1]).copyProperties(img, img.propertyNames());
}
function colT1T2(idT1, idT2, prepFn){
  var t1 = ee.ImageCollection(idT1).filterBounds(roi).map(maskC02L2).map(prepFn).map(function(i){return i.clip(roi);});
  var t2 = ee.ImageCollection(idT2).filterBounds(roi).map(maskC02L2).map(prepFn).map(function(i){return i.clip(roi);});
  return t1.merge(t2);
}

var lt04 = colT1T2("LANDSAT/LT04/C02/T1_L2","LANDSAT/LT04/C02/T2_L2", prepL457L2);
var lt05 = colT1T2("LANDSAT/LT05/C02/T1_L2","LANDSAT/LT05/C02/T2_L2", prepL457L2);
var le07 = colT1T2("LANDSAT/LE07/C02/T1_L2","LANDSAT/LE07/C02/T2_L2", prepL457L2);
var lc08 = colT1T2("LANDSAT/LC08/C02/T1_L2","LANDSAT/LC08/C02/T2_L2", prepL8L2);
var landsat = lt04.merge(lt05).merge(le07).merge(lc08);

function clearFraction(comp){
  var valid = comp.select("green").mask();
  return ee.Number(valid.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roi,
    scale: 120,
    bestEffort: true,
    tileScale: 16
  }).get("green"));
}
// Connect the Pixels together
function despeckle(mask){
  var cc = mask.selfMask().connectedPixelCount(1024, true);
  return mask.updateMask(cc.gte(200));
}
// Define the Snowmask for every year (Part 1.2)
function glacierMaskForYear(year){
  year = ee.Number(year);

  var startS = ee.Date.fromYMD(year, seasonStartMonth, seasonStartDay);
  var endS   = ee.Date.fromYMD(year.add(1), seasonEndMonth, seasonEndDay).advance(1,"day");

  var col = landsat.filterDate(startS, endS);
  var n = col.size();

  var comp = ee.Image(ee.Algorithms.If(n.gt(0), col.median(), ee.Image(0)));

  var cf = ee.Number(ee.Algorithms.If(n.gt(0), clearFraction(comp), 0));

  var ndsi = ee.Image(ee.Algorithms.If(
    n.gt(0),
    comp.normalizedDifference(["green","swir1"]).rename("NDSI"),
    ee.Image(0).rename("NDSI")
  ));

  var rawMask = ndsi.gt(ndsiThr).and(comp.select("swir1").lt(swir1Max));
  var ok = n.gte(minScenes_hypso).and(cf.gte(minClearFrac_hypso));

  var m = despeckle(rawMask).rename("mask");

  return ee.Image(ee.Algorithms.If(
    ok,
    m.unmask(0).clip(roi),
    ee.Image(0).rename("mask")
  ));
}

//Hypsometry: Define Mean per bin per year
function hypsometryForYear(year, year_group){
  year = ee.Number(year);
  var mask = glacierMaskForYear(year);

  var maskedBins = elevBin.updateMask(mask);
  var areaKm2 = ee.Image.pixelArea().divide(1000000).rename("area_km2");
  var areaInMask = areaKm2.updateMask(mask);

  var dict = areaInMask.addBands(maskedBins).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: "elev_bin"}),
    geometry: roi,
    scale: 60,          
    bestEffort: true,
    tileScale: 16
  });

  var groups = ee.List(dict.get("groups"));

  return ee.FeatureCollection(groups.map(function(g){
    g = ee.Dictionary(g);
    var b = ee.Number(g.get("elev_bin")).multiply(binSize);
    var a = ee.Number(g.get("sum"));
    return ee.Feature(null, {
      year: year,
      year_group: year_group,
      elev_bin: b,
      area_km2: a
    });
  }));
}
// Function for Hypsometry per Period
function hypsometryForPeriod(y0, y1, year_group){
  var years = ee.List.sequence(y0, y1);

  var all = ee.FeatureCollection(years.map(function(y){
    return hypsometryForYear(y, year_group);
  })).flatten();

  // mean area per bin across years (only where bin exists)
  var bins = ee.List(all.aggregate_array("elev_bin")).distinct();

  return ee.FeatureCollection(bins.map(function(b){
    b = ee.Number(b);
    var sub = all.filter(ee.Filter.eq("elev_bin", b));
    var meanA = ee.Number(sub.aggregate_mean("area_km2"));
    var nObs  = ee.Number(sub.size());
    return ee.Feature(null, {
      elev_bin: b,
      area_km2: meanA,
      year_group: year_group,
      n_years_in_bin: nObs
    });
  }));
}
// Clip Hypsometry on defined Groups
var hypsoEarly = hypsometryForPeriod(groupA.y0, groupA.y1, groupA.name); //pre2006
var hypsoLate  = hypsometryForPeriod(groupB.y0, groupB.y1, groupB.name); //post2006

var hypsoAll = hypsoEarly.merge(hypsoLate)
  .map(function(f){
    var g = ee.String(f.get("year_group"));
    var b = ee.Number(f.get("elev_bin")).format("%05d"); // make a int (%d) out of a string with a predefined length (05)
    return f.set("sort_key", g.cat("_").cat(b));
  })
  .sort("sort_key");

print("Hypsometry sample", hypsoAll.limit(10));
Map.addLayer(dem, {min: 0, max: 3000}, "DEM");

// export the Hypsometry layers (csv)
Export.table.toDrive({
  collection: hypsoAll,
  description: "Fox_hypsometry_bins_" + groupA.name + "_vs_" + groupB.name,
  folder: "GEE_Exports",
  fileNamePrefix: "Fox_hypsometry_bins",
  fileFormat: "CSV"
});

// validYears summary
function isValidYear(year){
  year = ee.Number(year);
  var m = glacierMaskForYear(year);
  var sum = ee.Number(m.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 120,
    bestEffort: true,
    tileScale: 16
  }).get("mask"));
  return sum.gt(0);
}

function validYearsList(y0, y1){
  var years = ee.List.sequence(y0, y1);
  var maybe = years.map(function(y){
    y = ee.Number(y);
    return ee.Algorithms.If(isValidYear(y), y, null);
  });
  return ee.List(maybe).removeAll([null]);
}

function groupValidStats(groupObj){
  var y0 = ee.Number(groupObj.y0);
  var y1 = ee.Number(groupObj.y1);
  var name = ee.String(groupObj.name);

  var valids = validYearsList(y0, y1);
  var nValid = valids.size();

  var yearsStr = ee.String(valids.map(function(y){
    return ee.Number(y).format();
  }).join(","));

  return ee.Feature(null, {
    year_group: name,
    y0: y0,
    y1: y1,
    n_valid_years: nValid,
    valid_years: yearsStr
  });
}

var validStats = ee.FeatureCollection([
  groupValidStats(groupA),
  groupValidStats(groupB)
]);

print("Valid year stats", validStats);

Export.table.toDrive({
  collection: validStats,
  description: "Fox_hypsometry_validYears_" + groupA.name + "_vs_" + groupB.name,
  folder: "GEE_Exports",
  fileNamePrefix: "Fox_hypsometry_validYears",
  fileFormat: "CSV"
});
