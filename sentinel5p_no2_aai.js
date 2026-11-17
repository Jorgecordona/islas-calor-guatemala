var zonasEstudio = ee.FeatureCollection("projects/icu-guate-ciudad/assets/Zonas14910_Guate_city");

Map.centerObject(zonasEstudio, 12);
Map.addLayer(zonasEstudio.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'Zonas 1-4-9-10');
print('Número de zonas cargadas:', zonasEstudio.size());
// ==============================
// 1. Área de estudio (tu asset)
// ==============================
var zonasEstudio = ee.FeatureCollection("projects/icu-guate-ciudad/assets/Zonas14910_Guate_city");

Map.centerObject(zonasEstudio, 12);
Map.addLayer(zonasEstudio.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'Zonas 1-4-9-10');

// ==============================
// 2. Periodo de análisis (2024-2025)
// ==============================
var fechaInicio = '2024-01-01';
var fechaFin = '2025-12-31';

// ==============================
// 3. Colección Sentinel-2 nivel superficie
// ==============================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(zonasEstudio)
  .filterDate(fechaInicio, fechaFin)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

print('Escenas Sentinel-2 disponibles 2024-2025:', s2.size());

// ==============================
// 4. Máscara de nubes (banda QA60)
// ==============================
function enmascararNubesS2(imagen) {
  var qa = imagen.select('QA60');
  var bitNubes = 1 << 10;
  var bitCirros = 1 << 11;
  var mascara = qa.bitwiseAnd(bitNubes).eq(0)
    .and(qa.bitwiseAnd(bitCirros).eq(0));
  return imagen.updateMask(mascara).divide(10000)
    .copyProperties(imagen, ['system:time_start']);
}

var s2Limpio = s2.map(enmascararNubesS2);

// ==============================
// 5. Composite de una fecha objetivo (empezamos con temporada seca)
// ==============================
var composite = s2Limpio
  .filterDate('2024-03-01', '2024-04-01')
  .median()
  .clip(zonasEstudio);

Map.addLayer(composite, {bands: ['B4','B3','B2'], min: 0, max: 0.3}, 'Sentinel-2 Marzo 2024 (RGB)');

// ==============================
// 6. NDVI
// ==============================
var ndvi = composite.normalizedDifference(['B8', 'B4']).rename('NDVI');
Map.addLayer(ndvi, {min: -1, max: 1, palette: ['brown','white','green']}, 'NDVI Marzo 2024');

// ==============================
// 7. Exportar a Drive
// ==============================
Export.image.toDrive({
  image: composite.select(['B2','B3','B4','B8']),
  description: 'S2_Marzo2024_ZonasEstudio',
  folder: 'GEE_IslasCalor',
  region: zonasEstudio,
  scale: 10,
  crs: 'EPSG:32615',
  maxPixels: 1e13
});
// ==============================
// 1. Área de estudio (tu asset)
// ==============================
var zonasEstudio = ee.FeatureCollection("projects/icu-guate-ciudad/assets/Zonas14910_Guate_city");

Map.centerObject(zonasEstudio, 12);
Map.addLayer(zonasEstudio.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'Zonas 1-4-9-10');

// ==============================
// 2. Periodo de análisis (2024-2025)
// ==============================
var fechaInicio = '2024-01-01';
var fechaFin = '2025-12-31';

// ==============================
// 3. Colección Sentinel-2 nivel superficie
// ==============================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(zonasEstudio)
  .filterDate(fechaInicio, fechaFin)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

print('Escenas Sentinel-2 disponibles 2024-2025:', s2.size());

// ==============================
// 4. Máscara de nubes (banda QA60)
// ==============================
function enmascararNubesS2(imagen) {
  var qa = imagen.select('QA60');
  var bitNubes = 1 << 10;
  var bitCirros = 1 << 11;
  var mascara = qa.bitwiseAnd(bitNubes).eq(0)
    .and(qa.bitwiseAnd(bitCirros).eq(0));
  return imagen.updateMask(mascara).divide(10000)
    .copyProperties(imagen, ['system:time_start']);
}

var s2Limpio = s2.map(enmascararNubesS2);

// ==============================
// 5. Función reutilizable: genera composite + NDVI + tarea de exportación
// ==============================
function procesarComposite(fechaIni, fechaFinComposite, etiqueta) {
  var composite = s2Limpio
    .filterDate(fechaIni, fechaFinComposite)
    .median()
    .clip(zonasEstudio);

  Map.addLayer(composite, {bands: ['B4','B3','B2'], min: 0, max: 0.3}, 'RGB ' + etiqueta);

  var ndvi = composite.normalizedDifference(['B8', 'B4']).rename('NDVI');
  Map.addLayer(ndvi, {min: -1, max: 1, palette: ['brown','white','green']}, 'NDVI ' + etiqueta);

  Export.image.toDrive({
    image: composite.select(['B2','B3','B4','B8']),
    description: 'S2_' + etiqueta,
    folder: 'GEE_IslasCalor',
    region: zonasEstudio,
    scale: 10,
    crs: 'EPSG:32615',
    maxPixels: 1e13
  });
}

// ==============================
// 6. Fechas representativas: verano/invierno x 2024/2025
// ==============================
procesarComposite('2024-03-01', '2024-04-01', 'Verano2024');
procesarComposite('2024-08-01', '2024-09-01', 'Invierno2024');
procesarComposite('2025-03-01', '2025-04-01', 'Verano2025');
procesarComposite('2025-08-01', '2025-09-01', 'Invierno2025');
// ==============================
// Sentinel-5P: NO2 troposférico y Aerosol Index (AAI)
// ==============================
// (Requiere que "zonasEstudio" ya esté definido, igual que en tu script de Sentinel-2)
// var zonasEstudio = ee.FeatureCollection("projects/icu-guate-ciudad/assets/Zonas14910_Guate_city");

// --- Colección NO2 ---
var s5pNO2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2')
  .select(['tropospheric_NO2_column_number_density', 'cloud_fraction'])
  .filterBounds(zonasEstudio);

function enmascararNO2(imagen) {
  var mascaraCalidad = imagen.select('cloud_fraction').lt(0.3);
  return imagen.select('tropospheric_NO2_column_number_density')
    .updateMask(mascaraCalidad)
    .copyProperties(imagen, ['system:time_start']);
}

var no2Limpio = s5pNO2.map(enmascararNO2);

// --- Colección Aerosol Index ---
var s5pAAI = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_AER_AI')
  .select('absorbing_aerosol_index')
  .filterBounds(zonasEstudio);

// ==============================
// Función: composite mensual (promedio) + diagnóstico + exportación
// ==============================
function procesarPolucion(coleccion, fechaIni, fechaFinRango, etiqueta, paleta, rangoVis) {
  var coleccionFecha = coleccion.filterDate(fechaIni, fechaFinRango);

  print('Escenas para ' + etiqueta + ':', coleccionFecha.size());

  var composite = coleccionFecha.mean().clip(zonasEstudio);

  var valorPromedio = composite.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: zonasEstudio,
    scale: 1113,
    maxPixels: 1e9
  });
  print('Valor promedio ' + etiqueta + ':', valorPromedio);

  Map.addLayer(composite, {min: rangoVis[0], max: rangoVis[1], palette: paleta}, etiqueta);

  Export.image.toDrive({
    image: composite,
    description: etiqueta,
    folder: 'GEE_IslasCalor',
    region: zonasEstudio,
    scale: 1113, // resolución nativa aproximada de Sentinel-5P
    crs: 'EPSG:32615',
    maxPixels: 1e13
  });
}

var paletaNO2 = ['black','blue','purple','cyan','green','yellow','red'];
var paletaAAI = ['blue','white','red'];

// ==============================
// NO2: las mismas 4 fechas representativas
// ==============================
procesarPolucion(no2Limpio, '2024-03-01','2024-04-01','NO2_Verano2024', paletaNO2, [0, 0.0002]);
procesarPolucion(no2Limpio, '2024-08-01','2024-09-01','NO2_Invierno2024', paletaNO2, [0, 0.0002]);
procesarPolucion(no2Limpio, '2025-03-01','2025-04-01','NO2_Verano2025', paletaNO2, [0, 0.0002]);
procesarPolucion(no2Limpio, '2025-08-01','2025-09-01','NO2_Invierno2025', paletaNO2, [0, 0.0002]);

// ==============================
// Aerosol Index (AAI): las mismas 4 fechas
// ==============================
procesarPolucion(s5pAAI, '2024-03-01','2024-04-01','AAI_Verano2024', paletaAAI, [-1, 2]);
procesarPolucion(s5pAAI, '2024-08-01','2024-09-01','AAI_Invierno2024', paletaAAI, [-1, 2]);
procesarPolucion(s5pAAI, '2025-03-01','2025-04-01','AAI_Verano2025', paletaAAI, [-1, 2]);
procesarPolucion(s5pAAI, '2025-08-01','2025-09-01','AAI_Invierno2025', paletaAAI, [-1, 2]);
