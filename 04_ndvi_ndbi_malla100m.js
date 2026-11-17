/**************************************************************************
 * 04_ndvi_ndbi_malla100m.js
 * NDVI y NDBI (Sentinel-2 L2A) por período estacional, promediados en cada
 * celda de la malla de 100 m (grid_lst_periodos). Insumo de la correlación
 * a nivel de celda LST–NDVI–NDBI (actividad 14) y del OLS/GWR (actividad 15).
 *
 * Si ya tenés NDVI y NDBI por período dentro de la malla, no necesitás este script.
 **************************************************************************/

// ======================= PARÁMETROS =======================
// Subí grid_lst_periodos.shp como asset y pegá aquí la ruta
var malla = ee.FeatureCollection('projects/icu-guate-ciudad/assets/grid_lst_periodos');
var CAMPO_ID = 'id';                 // identificador único de cada celda (revisalo en QGIS)
var CARPETA_DRIVE = 'GEE_ICU';

var aoi = malla.geometry().bounds();

// Mismos períodos que la LST
var periodos = [
  {nombre: 'Ver24', ini: '2024-01-01', fin: '2024-05-01'},
  {nombre: 'Inv24', ini: '2024-05-01', fin: '2024-11-01'},
  {nombre: 'Ver25', ini: '2025-01-01', fin: '2025-05-01'},
  {nombre: 'Inv25', ini: '2025-05-01', fin: '2025-11-01'}
];

// ======================= SENTINEL-2 =======================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate('2024-01-01', '2025-11-01')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60));

var s2nubes = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
  .filterBounds(aoi)
  .filterDate('2024-01-01', '2025-11-01');

var s2unida = ee.ImageCollection(ee.Join.saveFirst('prob_nubes').apply({
  primary: s2,
  secondary: s2nubes,
  condition: ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'})
}));

function indicesS2(img) {
  var scl = img.select('SCL');
  var okScl = scl.neq(1).and(scl.neq(3)).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  var prob = ee.Image(img.get('prob_nubes')).select('probability');
  var sr = img.select(['B4', 'B8', 'B11']).multiply(0.0001);
  var ndvi = sr.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndbi = sr.normalizedDifference(['B11', 'B8']).rename('NDBI');
  return ndvi.addBands(ndbi).updateMask(okScl.and(prob.lt(40)))
    .copyProperties(img, ['system:time_start']);
}
var s2indices = s2unida.map(indicesS2);

// ======================= COMPUESTOS POR PERÍODO =======================
var pila = ee.Image.cat(periodos.map(function (p) {
  var med = s2indices.filterDate(p.ini, p.fin).median();
  return med.rename(['NDVI_' + p.nombre, 'NDBI_' + p.nombre]);
}));

// Media por celda de 100 m, calculada a 10 m (resolución nativa de B4/B8)
var resultado = pila.reduceRegions({
  collection: malla.select([CAMPO_ID]),
  reducer: ee.Reducer.mean(),
  scale: 10,
  tileScale: 4
}).map(function (f) { return f.setGeometry(null); });

print('Vista previa:', resultado.limit(5));

var columnas = [CAMPO_ID];
periodos.forEach(function (p) { columnas.push('NDVI_' + p.nombre, 'NDBI_' + p.nombre); });

Export.table.toDrive({
  collection: resultado,
  description: 'ndvi_ndbi_malla100m',
  folder: CARPETA_DRIVE,
  fileFormat: 'CSV',
  selectors: columnas
});

Map.centerObject(aoi, 13);
Map.addLayer(pila.select('NDVI_Ver24'), {min: 0, max: 0.8, palette: ['#d73027', '#f7f7f7', '#1a9850']}, 'NDVI Ver24');
Map.addLayer(malla.style({color: '555555', fillColor: '00000000', width: 0.3}), {}, 'Malla 100 m');
