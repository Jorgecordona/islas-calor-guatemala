// ==============================
// 1. Área de estudio (tu asset)
// ==============================
var zonasEstudio = ee.FeatureCollection("projects/icu-guate-ciudad/assets/Zonas14910_Guate_city");

Map.centerObject(zonasEstudio, 12);
Map.addLayer(zonasEstudio.style({color: 'red', fillColor: '00000000', width: 2}), {}, 'Zonas 1-4-9-10');

// ==============================
// 2. Colección Sentinel-2 nivel superficie (sin filtro de nubosidad de escena completa;
// confiamos en el enmascarado por píxel QA60 más abajo)
// ==============================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(zonasEstudio)
  .filterDate('2024-01-01', '2025-12-31');

print('Escenas Sentinel-2 disponibles 2024-2025:', s2.size());

// ==============================
// 3. Máscara de nubes (banda QA60)
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
// 4. Función reutilizable: composite + NDVI + exportación
// ==============================
function procesarComposite(fechaIni, fechaFinComposite, etiqueta) {
  var coleccionFecha = s2Limpio.filterDate(fechaIni, fechaFinComposite);

  print('Escenas para ' + etiqueta + ':', coleccionFecha.size());

  var composite = coleccionFecha.median().clip(zonasEstudio);

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
// 5. Temporadas completas en vez de un solo mes
// ==============================
// Guatemala: "verano" = temporada seca (ene-abril, dentro del periodo de estudio) 
//            "invierno" = temporada lluviosa (mayo-octubre)
procesarComposite('2024-01-01', '2024-04-30', 'Verano2024');
procesarComposite('2024-05-01', '2024-10-31', 'Invierno2024');
procesarComposite('2025-01-01', '2025-04-30', 'Verano2025');
procesarComposite('2025-05-01', '2025-10-31', 'Invierno2025');
