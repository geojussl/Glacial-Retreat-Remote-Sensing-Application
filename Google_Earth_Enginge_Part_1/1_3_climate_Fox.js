/*
Remote Sensing Applications: Group Glacial Retreat
Part 1.3 (GEE): Climate Data Extraction 
therefore we use ECMWF ERA5-Land monthly (temperature_2m, total_precipitation), ROI_Fox

Output: climate variables (Preciptation, temperature) (CSV)

Course: Remote Sensing Applications
Group: Glacial Retreat
By Justin Lingg-Laham, Liev Paustian, Lisa Wohlgemuth

Workflow:

Two time windows per year (See Part 1.2):
A) Season (Aug–May, hydrological “glacier season”)
B) Summer (Jan–Mar, melt season)

1. Load ERA5-Land monthly and convert units
- temperature_2m (Kelvin)  → converted to °C  (t2m_C)
- total_precipitation (meters) → converted to mm (tp_mm)

2. Reduce climate signals over the ROI
For each year:
A) Season window (Aug y → May y+1)
- t2m_C_mean_season: mean 2 m temperature (°C) over ROI
- tp_mm_sum_season: sum of monthly precipitation totals (mm) over ROI
- nMonths_season: number of monthly images used (For Data Availability)

B) Summer window (Jan–Mar of the same year)
- t2m_C_mean_summer: mean 2 m temperature (°C) over ROI
- tp_mm_sum_summer: sum of monthly precipitation totals (mm) over ROI
- t2m_C_pos_sum_summer: positive temperature sum (proxy)
  (monthly mean temperatures clamped at 0°C and summed; units ≈ °C-months)
- nMonths_summer: number of monthly images used

Outputs (files):
- Fox_Climate_ERA5Land_1989_2025.csv
  Columns include:
  year, nMonths_season, nMonths_summer,
  t2m_C_mean_season, tp_mm_sum_season,
  t2m_C_mean_summer, tp_mm_sum_summer,
  t2m_C_pos_sum_summer,
  season, summer, dataset
*/

// define ROI
var bufferAnalysis_m = 300;
var bufferSearch_m   = 8000;
var simplify_m       = 20;

var roiBase = ee.FeatureCollection(fox).geometry();
var roiAnalysis = roiBase.buffer(bufferAnalysis_m);
if (simplify_m > 0) roiAnalysis = roiAnalysis.simplify(simplify_m);
var roiSearch = roiBase.buffer(bufferSearch_m);

Map.centerObject(roiBase, 10);

// define Windows
var seasonStartMonth = 8, seasonStartDay = 1;   // 01 Aug
var seasonEndMonth   = 5, seasonEndDay   = 31;  // 31 May 

var summerStartMonth = 1, summerStartDay = 1;   // 01 Jan
var summerEndMonth   = 3, summerEndDay   = 31;  // 31 Mar

// define Timeline
var startYear = 1989;
var endYear   = 2025;
var years = ee.List.sequence(startYear, endYear);

// Get Climate Data (ERA5-Land)
var era = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY")
  .filterBounds(roiSearch);

// convert Units
function toT2mC(img){
  // kelvin to celsius
  return img.select("temperature_2m").subtract(273.15).rename("t2m_C")
    .copyProperties(img, img.propertyNames());
}
function toTPmm(img){
  // total_preciptation in mm
  return img.select("total_precipitation").multiply(1000).rename("tp_mm")
    .copyProperties(img, img.propertyNames());
}

// Mean over ROI 
function meanOverROI(img, bandName){
  var v = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roiAnalysis,
    scale: 1000,        
    bestEffort: true,
    tileScale: 4
  }).get(bandName);
  return v;
}

function sumOverROIMulti(col, bandName){
  var imgSum = col.select(bandName).sum();
  return meanOverROI(imgSum, bandName);
}

function meanOverROIMulti(col, bandName){
  var imgMean = col.select(bandName).mean();
  return meanOverROI(imgMean, bandName);
}

// positive temperature sum (proxy) (monthly means > 0)
function posTempSum(colT){
  // colT is monthly temp in °C
  var pos = colT.map(function(img){
    var t = img.select("t2m_C");
    var tpos = t.max(0); 
    return tpos.rename("tpos").copyProperties(img, img.propertyNames());
  });
  var sumImg = pos.select("tpos").sum();
  return meanOverROI(sumImg, "tpos"); // units: °C-months (proxy)
}

// Statistics per Year
function climateForYear(y){
  y = ee.Number(y);

  // Season Aug-May
  var startS = ee.Date.fromYMD(y, seasonStartMonth, seasonStartDay);
  var endS   = ee.Date.fromYMD(y.add(1), seasonEndMonth, seasonEndDay).advance(1,"day");

  // Season jan-Mar
  var startU = ee.Date.fromYMD(y, summerStartMonth, summerStartDay);
  var endU   = ee.Date.fromYMD(y, summerEndMonth, summerEndDay).advance(1,"day");

  // build collections
  var colS = era.filterDate(startS, endS);
  var colU = era.filterDate(startU, endU);

  var nS = colS.size();
  var nU = colU.size();

  // Temperature (°C mean)
  var tS = colS.map(toT2mC);
  var tU = colU.map(toT2mC);

  var t2m_C_mean_season = ee.Algorithms.If(nS.gt(0), meanOverROIMulti(tS, "t2m_C"), null);
  var t2m_C_mean_summer = ee.Algorithms.If(nU.gt(0), meanOverROIMulti(tU, "t2m_C"), null);

  // Preciptation (mm)
  var pS = colS.map(toTPmm);
  var pU = colU.map(toTPmm);

  var tp_mm_sum_season = ee.Algorithms.If(nS.gt(0), sumOverROIMulti(pS, "tp_mm"), null);
  var tp_mm_sum_summer = ee.Algorithms.If(nU.gt(0), sumOverROIMulti(pU, "tp_mm"), null);

  // Proxy (monthly mean positive temp sum)
  var t2m_C_pos_sum_summer = ee.Algorithms.If(nU.gt(0), posTempSum(tU), null);

  return ee.Feature(null, {
    year: y,

    nMonths_season: nS,
    nMonths_summer: nU,

    t2m_C_mean_season: t2m_C_mean_season,
    tp_mm_sum_season: tp_mm_sum_season,

    t2m_C_mean_summer: t2m_C_mean_summer,
    tp_mm_sum_summer: tp_mm_sum_summer,

    t2m_C_pos_sum_summer: t2m_C_pos_sum_summer,

    season: seasonStartMonth + "/" + seasonStartDay + "–" + seasonEndMonth + "/" + seasonEndDay,
    summer: summerStartMonth + "/" + summerStartDay + "–" + summerEndMonth + "/" + summerEndDay,
    dataset: "ERA5-Land MONTHLY"
  });
}

var climateFC = ee.FeatureCollection(years.map(climateForYear));
print("Climate sample", climateFC.limit(5));

// Export
Export.table.toDrive({
  collection: climateFC,
  description: "Fox_Climate_ERA5Land_" + startYear + "_" + endYear,
  folder: "GEE_Exports",
  fileNamePrefix: "Fox_Climate_ERA5Land_" + startYear + "_" + endYear,
  fileFormat: "CSV"
});
