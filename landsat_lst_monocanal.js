// ==============================
// 1. Área de estudio
// ==============================
var zonasEstudio = ee.FeatureCollection("projects/icu-guate-ciudad/assets/Zonas14910_Guate_city");

Map.centerObject(zonasEstudio, 12);
Map.addLayer(zonasEstudio.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'Zonas 1-4-9-10');

// ==============================
// 2. Dos colecciones: TOA (temperatura de brillo) y Nivel-2 (reflectancia corregida)
// ==============================
var l8Toa = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA').filterBounds(zonasEstudio);
var l9Toa = ee.ImageCollection('LANDSAT/LC09/C02/T1_TOA').filterBounds(zonasEstudio);
var toaColeccion = l8Toa.merge(l9Toa).filterDate('2024-01-01', '2025-12-31');

var l8Sr = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(zonasEstudio);
var l9Sr = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(zonasEstudio);
var srColeccion = l8Sr.merge(l9Sr).filterDate('2024-01-01', '2025-12-31');

print('Escenas TOA disponibles:', toaColeccion.size());
print('Escenas Nivel-2 (SR) disponibles:', srColeccion.size());

// Unimos cada escena TOA con su par de Reflectancia de Superficie (mismo ID de escena)
var filtroUnion = ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'});
var unidas = ee.Join.inner().apply(toaColeccion, srColeccion, filtroUnion);

var landsat = ee.ImageCollection(unidas.map(function(par) {
  var toa = ee.Image(par.get('primary'));
  var sr = ee.Image(par.get('secondary'));
  return toa.addBands(sr.select(['SR_B4', 'SR_B5']));
}));

print('Escenas unidas (TOA + Nivel-2):', landsat.size());

// ==============================
// 3. Máscara de nubes (QA_PIXEL de la colección TOA)
// ==============================
function enmascararNubesLandsat(imagen) {
  var qa = imagen.select('QA_PIXEL');
  var bitNubeDilatada = 1 << 1;
  var bitCirros = 1 << 2;
  var bitNube = 1 << 3;
  var bitSombraNube = 1 << 4;
  var mascara = qa.bitwiseAnd(bitNubeDilatada).eq(0)
    .and(qa.bitwiseAnd(bitCirros).eq(0))
    .and(qa.bitwiseAnd(bitNube).eq(0))
    .and(qa.bitwiseAnd(bitSombraNube).eq(0));
  return imagen.updateMask(mascara);
}

// ==============================
// 4. Cálculo de LST por imagen (NDVI desde Nivel-2 + emisividad + mono-window)
// ==============================
function calcularLST(imagen) {
  imagen = enmascararNubesLandsat(imagen);

  // Temperatura de brillo desde la banda térmica TOA (ya en Kelvin)
  var bt = imagen.select('B10').rename('BT');

  // NDVI desde reflectancia de superficie Nivel-2 (con factor de escala oficial)
  var srEscalado = imagen.select(['SR_B4', 'SR_B5']).multiply(0.0000275).add(-0.2);
  var ndvi = srEscalado.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');

  var pv = ndvi.subtract(0.2).divide(0.3).pow(2).rename('PV');
  var emisividad = ndvi.expression(
    '(NDVI < 0.2) ? 0.97 : (NDVI > 0.5) ? 0.99 : (0.004 * PV + 0.986)', {
      'NDVI': ndvi, 'PV': pv
    }).rename('EM');

  var lambda = 10.895;
  var rho = 14380;

  var lst = bt.expression(
    '(BT / (1 + (lambda * (BT / rho)) * log(EM))) - 273.15', {
      'BT': bt, 'EM': emisividad, 'lambda': lambda, 'rho': rho
    }).rename('LST');

  return lst.copyProperties(imagen, ['system:time_start']);
}

var lstColeccion = landsat.map(calcularLST);

// ==============================
// 5. Función: composite (mediana) por temporada + diagnóstico + exportación
// ==============================
function procesarLST(fechaIni, fechaFinRango, etiqueta) {
  var coleccionFecha = lstColeccion.filterDate(fechaIni, fechaFinRango);

  print('Escenas para ' + etiqueta + ':', coleccionFecha.size());

  var composite = coleccionFecha.median().clip(zonasEstudio);

  var minMax = composite.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: zonasEstudio,
    scale: 30,
    maxPixels: 1e9
  });
  print('Rango LST ' + etiqueta + ' (°C):', minMax);

  Map.addLayer(composite, {min: 19, max: 37, palette: ['blue','green','yellow','orange','red']}, 'LST ' + etiqueta);

  Export.image.toDrive({
    image: composite,
    description: 'LST_' + etiqueta,
    folder: 'GEE_IslasCalor',
    region: zonasEstudio,
    scale: 30,
    crs: 'EPSG:32615',
    maxPixels: 1e13
  });
}

// ==============================
// 6. Mismas 4 temporadas
// ==============================
procesarLST('2024-01-01', '2024-04-30', 'Verano2024');
procesarLST('2024-05-01', '2024-10-31', 'Invierno2024');
procesarLST('2025-01-01', '2025-04-30', 'Verano2025');
procesarLST('2025-05-01', '2025-10-31', 'Invierno2025');
