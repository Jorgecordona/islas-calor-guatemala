// ========================================================================
// ISLAS DE CALOR URBANO - CIUDAD DE GUATEMALA (Zonas 1, 4, 9, 10)
// Fase 5.1.3 - Variables ambientales: NDVI y NDBI (Sentinel-2)
// Para correlación posterior con LST (Pearson/Spearman, OLS/GWR)
// ========================================================================

// --------------------------------------------------------------------
// 1. ÁREA DE ESTUDIO
// --------------------------------------------------------------------
var zonasEstudio = ee.FeatureCollection('projects/icu-guate-ciudad/assets/Zonas14910_Guate_city');
Map.centerObject(zonasEstudio, 12);
Map.addLayer(zonasEstudio, {color: 'red'}, 'Zonas de estudio');

// --------------------------------------------------------------------
// 2. PERIODOS ESTACIONALES (mismos que en los scripts de LST)
// --------------------------------------------------------------------
var periodos = {
  'Verano2024':   ['2024-01-01', '2024-04-30'],
  'Invierno2024': ['2024-05-01', '2024-10-31'],
  'Verano2025':   ['2025-01-01', '2025-04-30'],
  'Invierno2025': ['2025-05-01', '2025-10-31']
};

// --------------------------------------------------------------------
// 3. MÁSCARA DE NUBES (QA60, Sentinel-2 SR)
// --------------------------------------------------------------------
function enmascararNubesS2(imagen) {
  var qa = imagen.select('QA60');
  var nubes = 1 << 10;
  var cirros = 1 << 11;
  var mask = qa.bitwiseAnd(nubes).eq(0).and(qa.bitwiseAnd(cirros).eq(0));
  return imagen.updateMask(mask).divide(10000); // reflectancia a escala 0-1
}

// --------------------------------------------------------------------
// 4. CÁLCULO DE NDVI Y NDBI POR ESCENA
// NDVI = (NIR - Red) / (NIR + Red)       -> bandas B8, B4
// NDBI = (SWIR1 - NIR) / (SWIR1 + NIR)   -> bandas B11, B8
// --------------------------------------------------------------------
function calcularIndices(imagen) {
  var ndvi = imagen.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndbi = imagen.normalizedDifference(['B11', 'B8']).rename('NDBI');
  return imagen.addBands(ndvi).addBands(ndbi);
}

// --------------------------------------------------------------------
// 5. PROCESAMIENTO POR PERIODO: composición, diagnóstico y exportación
// --------------------------------------------------------------------
Object.keys(periodos).forEach(function(nombrePeriodo) {
  var fechas = periodos[nombrePeriodo];
  var inicio = fechas[0];
  var fin = fechas[1];

  var coleccion = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(zonasEstudio)
    .filterDate(inicio, fin)
    .map(enmascararNubesS2)
    .map(calcularIndices);

  print('Escenas disponibles ' + nombrePeriodo + ':', coleccion.size());

  var compuesto = coleccion.select(['NDVI', 'NDBI']).median().clip(zonasEstudio);

  var stats = compuesto.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: zonasEstudio.geometry(),
    scale: 10,
    maxPixels: 1e9
  });
  print('NDVI/NDBI min/max ' + nombrePeriodo + ':', stats);

  Map.addLayer(compuesto.select('NDVI'), {min: -0.2, max: 0.8, palette: ['brown', 'yellow', 'green']},
    'NDVI_' + nombrePeriodo);
  Map.addLayer(compuesto.select('NDBI'), {min: -0.5, max: 0.5, palette: ['green', 'white', 'red']},
    'NDBI_' + nombrePeriodo);

  Export.image.toDrive({
    image: compuesto,
    description: 'NDVI_NDBI_' + nombrePeriodo,
    folder: 'GEE_IslasCalor',
    fileNamePrefix: 'NDVI_NDBI_' + nombrePeriodo,
    region: zonasEstudio.geometry(),
    scale: 10,
    crs: 'EPSG:32615',
    maxPixels: 1e9
  });
});
