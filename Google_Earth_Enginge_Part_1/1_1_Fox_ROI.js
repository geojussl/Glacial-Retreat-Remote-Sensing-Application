/*
Remote Sensing Applications: Group Glacial Retreat
Part 1.1 (GEE): Define Region of Interest (ROI): Fox Glacier, New Zealand
therefore we use Global Land Ice Measurements from Space (GLIMS)

Output: ROI as Shapefile

Course: Remote Sensing Applications
Group: Glacial Retreat
By Justin Lingg-Laham, Liev Paustian, Lisa Wohlgemuth
*/

//find data (Search Radius)
var bufferMeters = 500;     
var searchRadius = 15000;   
var minYear = 1970;       

var foxPt = ee.Geometry.Point([169.96, -43.47]); 
var nzBox = ee.Geometry.Rectangle([166, -46, 173, -41]);

Map.setOptions("TERRAIN");
Map.centerObject(foxPt, 10);
Map.addLayer(foxPt, {color:"yellow"}, "Fox point", false);
Map.addLayer(nzBox, {color:"cyan"}, "nzBox", false);

// Load GLIMS (NASA)
var glims = ee.FeatureCollection("GLIMS/current").filterBounds(nzBox);

// Function for geometry area in km²
function addGeomAreaKm2(f) {
  return f.set("geom_km2", ee.Number(f.geometry().area()).divide(1000000));
}

// Find candidates near to ROI (Data-availability)
var cand = glims
  .filterBounds(foxPt.buffer(searchRadius))
  .filter(ee.Filter.eq("line_type", "glac_bound"))
  .map(addGeomAreaKm2);

print('Candidates near Fox', cand.size());

// Show Top Candidates (Overview)
var top = cand.sort("geom_km2", false).limit(30);
Map.addLayer(top.style({color:"cyan", fillColor:"00000000", width: 2}), {}, "Top near Fox (geom_km2)", true);
print("Top near Fox: ", top);

var biggestCand = ee.Feature(cand.sort("geom_km2", false).first());
var foxId = ee.String(biggestCand.get("glac_id"));
print("Selected glac_id (Fox target)", foxId);
print("Biggest candidate properties", biggestCand);

var foxAll = glims
  .filter(ee.Filter.eq("glac_id", foxId))
  .filter(ee.Filter.eq("line_type", "glac_bound"))
  .map(addGeomAreaKm2);
// Stats for Top candidates
print("foxAll size", foxAll.size());
print("foxAll dates", foxAll.aggregate_array("src_date"));
print("foxAll area property (dirty)", foxAll.aggregate_array("area"));
print("foxAll geom_km2", foxAll.aggregate_array("geom_km2"));

Map.addLayer(foxAll.style({color:"gold", fillColor:"00000000", width: 2}), {}, "FOX outlines (all)", true);

var cleaned = foxAll
  .filter(ee.Filter.gt("area", 0)) 
  .filter(ee.Filter.gte("src_date", minYear + "-01-01T00:00:00"));

print("cleaned size", cleaned.size());
print("cleaned dates", cleaned.aggregate_array("src_date"));
print("cleaned geom_km2", cleaned.aggregate_array("geom_km2"));

var cleanedSafe = ee.FeatureCollection(ee.Algorithms.If(cleaned.size().gt(0), cleaned, foxAll));
//ROI out of cleaned Data (Roi Big just for testing the data)
var biggestOutline = ee.Feature(cleanedSafe.sort("geom_km2", false).first());

var roiUnion = cleanedSafe.union().geometry().buffer(bufferMeters); 
var roiBig   = biggestOutline.geometry().buffer(bufferMeters);      
// Final ROI
var roi = roiUnion;

Map.centerObject(roi, 11);
Map.addLayer(roiUnion, {color:"00BFFF"}, "ROI UNION + buffer", true);
Map.addLayer(roiBig,   {color:"FF0000"}, "ROI BIGGEST + buffer", false);

// Area ROI
var roiAreaKm2 = ee.Number(roi.area()).divide(1000000);
print("ROI area (km2)", roiAreaKm2);

// Export ROI as Shapefile
// make a FeatureCollection
var roiFC = ee.FeatureCollection([
  ee.Feature(roi, {
    name: 'Fox_ROI_union_buffer',
    buffer_m: bufferMeters,
    minYear: minYear,
    searchRadius_m: searchRadius
  })
]);

// Export 
Export.table.toDrive({
  collection: roiFC,
  description: 'FOX_ROI_union_buffer_SHP',
  folder: 'GEE_Exports',
  fileNamePrefix: 'FOX_ROI_union_buffer',
  fileFormat: 'SHP'
});
