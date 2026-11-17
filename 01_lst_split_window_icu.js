// ========================================================================
// ISLAS DE CALOR URBANO - CIUDAD DE GUATEMALA (Zonas 1, 4, 9, 10)
// LST - Algoritmo SPLIT-WINDOW (Du et al., 2015)
// Landsat 8/9 Collection 2 - Bandas termales B10 y B11 (TOA, ya en Kelvin)
// Emisividad: método NDVI threshold (Sobrino / Skokovic et al., 2014),
//             usando reflectancia de superficie Nivel-2 (corregida atmosféricamente)
// Vapor de agua atmosférico (CWV): dinámico, NCEP/NCAR Reanalysis (pr_wtr)
// + Umbralización estadística para delimitación de zonas ICU (Fase 5.1.2)
// ========================================================================

// --------------------------------------------------------------------
// 1. ÁREA DE ESTUDIO
// --------------------------------------------------------------------
var zonasEstudio = ee.FeatureCollection('projects/icu-guate-ciudad/assets/Zonas14910_Guate_city');
Map.centerObject(zonasEstudio, 12);
Map.addLayer(zonasEstudio, {color: 'red'}, 'Zonas de estudio');

// --------------------------------------------------------------------
// 2. PERIODOS ESTACIONALES
// --------------------------------------------------------------------
var periodos = {
  'Verano2024':   ['2024-01-01', '2024-04-30'],
  'Invierno2024': ['2024-05-01', '2024-10-31'],
  'Verano2025':   ['2025-01-01', '2025-04-30'],
  'Invierno2025': ['2025-05-01', '2025-10-31']
};

// --------------------------------------------------------------------
// 3. MÁSCARA DE NUBES (QA_PIXEL, Landsat Collection 2)
// --------------------------------------------------------------------
function enmascararNubes(imagen) {
  var qa = imagen.select('QA_PIXEL');
  var dilatedCloud = 1 << 1;
  var cirrus = 1 << 2;
  var cloud = 1 << 3;
  var cloudShadow = 1 << 4;
  var mask = qa.bitwiseAnd(dilatedCloud).eq(0)
    .and(qa.bitwiseAnd(cirrus).eq(0))
    .and(qa.bitwiseAnd(cloud).eq(0))
    .and(qa.bitwiseAnd(cloudShadow).eq(0));
  return imagen.updateMask(mask);
}

// --------------------------------------------------------------------
// 4. COEFICIENTES DEL ALGORITMO SPLIT-WINDOW (Du et al., 2015, Tabla 3)
// --------------------------------------------------------------------
function obtenerCoeficientes(cwv) {
  var b0 = ee.Image.constant(-0.41165);
  var b1 = ee.Image.constant(1.00522);
  var b2 = ee.Image.constant(0.14543);
  var b3 = ee.Image.constant(-0.27297);
  var b4 = ee.Image.constant(4.06655);
  var b5 = ee.Image.constant(-6.92512);
  var b6 = ee.Image.constant(-18.27461);
  var b7 = ee.Image.constant(0.24468);

  var r1 = cwv.gte(0).and(cwv.lt(2.00));
  var r2 = cwv.gte(2.00).and(cwv.lt(3.25));
  var r3 = cwv.gte(3.25).and(cwv.lt(4.25));
  var r4 = cwv.gte(4.25).and(cwv.lt(5.20));
  var r5 = cwv.gte(5.20).and(cwv.lte(6.30));

  b0 = b0.where(r1, -2.78009).where(r2, 11.00824).where(r3, 9.62610).where(r4, 0.61258).where(r5, -0.34808);
  b1 = b1.where(r1, 1.01408).where(r2, 0.95995).where(r3, 0.96202).where(r4, 0.99124).where(r5, 0.98123);
  b2 = b2.where(r1, 0.15833).where(r2, 0.17243).where(r3, 0.13834).where(r4, 0.10051).where(r5, 0.05599);
  b3 = b3.where(r1, -0.34991).where(r2, -0.28852).where(r3, -0.17262).where(r4, -0.09664).where(r5, -0.03518);
  b4 = b4.where(r1, 4.04487).where(r2, 7.11492).where(r3, 7.87883).where(r4, 7.85758).where(r5, 11.96444);
  b5 = b5.where(r1, 3.55414).where(r2, 0.42684).where(r3, 5.17910).where(r4, 6.86626).where(r5, 9.06710);
  b6 = b6.where(r1, -8.88394).where(r2, -6.62025).where(r3, -13.26611).where(r4, -15.00742).where(r5, -14.74085);
  b7 = b7.where(r1, 0.09152).where(r2, -0.06381).where(r3, -0.07603).where(r4, -0.01185).where(r5, -0.20471);

  return {b0: b0, b1: b1, b2: b2, b3: b3, b4: b4, b5: b5, b6: b6, b7: b7};
}

// --------------------------------------------------------------------
// 5. VAPOR DE AGUA ATMOSFÉRICO (CWV) - NCEP/NCAR Reanalysis
// --------------------------------------------------------------------
function obtenerCWV(fechaImagen) {
  var fecha = ee.Date(fechaImagen);
  var coleccionCWV = ee.ImageCollection('NCEP_RE/surface_wv')
    .filterDate(fecha, fecha.advance(1, 'day'));
  var cwv_kgm2 = coleccionCWV.mean().select('pr_wtr').unmask(35);
  return cwv_kgm2.divide(10).rename('CWV'); // g/cm2
}

// --------------------------------------------------------------------
// 6. EMISIVIDAD POR BANDA (método del umbral de NDVI)
// --------------------------------------------------------------------
function calcularEmisividad(ndvi, epsSuelo, epsVeg) {
  var Pv = ndvi.subtract(0.2).divide(0.3).pow(2).clamp(0, 1);
  var F = 0.55;

  var mezcla = ee.Image.constant(epsVeg).multiply(Pv)
    .add(ee.Image.constant(epsSuelo).multiply(ee.Image.constant(1).subtract(Pv)))
    .add(ee.Image.constant(1).subtract(epsSuelo).multiply(epsVeg).multiply(F)
      .multiply(ee.Image.constant(1).subtract(Pv)));

  var emisividad = mezcla
    .where(ndvi.lt(0.2), epsSuelo)
    .where(ndvi.gt(0.5), epsVeg);

  return emisividad.rename('emisividad');
}

// --------------------------------------------------------------------
// 7. UNIR TOA (térmico) CON NIVEL-2 (reflectancia de superficie, para NDVI)
// --------------------------------------------------------------------
function unirColecciones(toa, l2) {
  var filtro = ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'});
  var unidas = ee.Join.inner().apply(toa, l2, filtro);
  return ee.ImageCollection(unidas.map(function(par) {
    return ee.Image.cat(par.get('primary'), par.get('secondary'));
  }));
}

// --------------------------------------------------------------------
// 8. CÁLCULO DE LST POR ESCENA (algoritmo Split-Window completo)
// --------------------------------------------------------------------
function calcularLST_SW(imagen) {
  var BTi = imagen.select('B10');
  var BTj = imagen.select('B11');

  var sr = imagen.select(['SR_B4', 'SR_B5']).multiply(0.0000275).add(-0.2);
  var ndvi = sr.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');

  var eps10 = calcularEmisividad(ndvi, 0.971, 0.987);
  var eps11 = calcularEmisividad(ndvi, 0.977, 0.989);
  var epsMedia = eps10.add(eps11).divide(2);
  var deltaEps = eps10.subtract(eps11);

  var cwv = obtenerCWV(imagen.date());
  var coef = obtenerCoeficientes(cwv);

  var promedioT = BTi.add(BTj).divide(2);
  var diferenciaT = BTi.subtract(BTj);

  var term1 = coef.b1
    .add(coef.b2.multiply(ee.Image.constant(1).subtract(epsMedia)).divide(epsMedia))
    .add(coef.b3.multiply(deltaEps).divide(epsMedia.pow(2)))
    .multiply(promedioT);

  var term2 = coef.b4
    .add(coef.b5.multiply(ee.Image.constant(1).subtract(epsMedia)).divide(epsMedia))
    .add(coef.b6.multiply(deltaEps).divide(epsMedia.pow(2)))
    .multiply(diferenciaT.divide(2));

  var term3 = coef.b7.multiply(diferenciaT.pow(2));

  var lst_K = coef.b0.add(term1).add(term2).add(term3);
  var lst_C = lst_K.subtract(273.15).rename('LST_SW');

  return imagen.addBands(lst_C).addBands(ndvi).addBands(epsMedia.rename('emisividad_media'));
}

// --------------------------------------------------------------------
// 9. PROCESAMIENTO POR PERIODO + 10. UMBRALIZACIÓN ESTADÍSTICA ICU
// --------------------------------------------------------------------
Object.keys(periodos).forEach(function(nombrePeriodo) {
  var fechas = periodos[nombrePeriodo];
  var inicio = fechas[0];
  var fin = fechas[1];

  var toa = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA')
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_TOA'))
    .filterBounds(zonasEstudio)
    .filterDate(inicio, fin)
    .map(enmascararNubes);

  var l2 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(zonasEstudio)
    .filterDate(inicio, fin)
    .select(['SR_B4', 'SR_B5']);

  var unidas = unirColecciones(toa, l2);
  var conLST = unidas.map(calcularLST_SW);

  print('Escenas disponibles ' + nombrePeriodo + ':', conLST.size());

  var compuestoLST = conLST.select('LST_SW').median().clip(zonasEstudio);

  var statsLST = compuestoLST.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: zonasEstudio.geometry(),
    scale: 30,
    maxPixels: 1e9
  });
  print('LST Split-Window min/max ' + nombrePeriodo + ' (°C):', statsLST);

  var cwvRef = obtenerCWV(ee.Date(inicio).advance(
    ee.Date(fin).difference(ee.Date(inicio), 'day').divide(2), 'day'));
  print('CWV aproximado a mitad de periodo ' + nombrePeriodo + ' (g/cm2):',
    cwvRef.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: zonasEstudio.geometry(),
      scale: 27830,
      maxPixels: 1e9
    }));

  Map.addLayer(compuestoLST, {min: 20, max: 40, palette: ['blue', 'yellow', 'red']},
    'LST_SW_' + nombrePeriodo);

  Export.image.toDrive({
    image: compuestoLST,
    description: 'LST_SplitWindow_' + nombrePeriodo,
    folder: 'GEE_IslasCalor',
    fileNamePrefix: 'LST_SplitWindow_' + nombrePeriodo,
    region: zonasEstudio.geometry(),
    scale: 30,
    crs: 'EPSG:32615',
    maxPixels: 1e9
  });

  // ------------------------------------------------------------------
  // 10. UMBRALIZACIÓN ESTADÍSTICA - DELIMITACIÓN DE ZONAS ICU (Fase 5.1.2)
  // Clasificación en 4 clases según desviaciones estándar respecto a la
  // media de LST del área de estudio (propia de cada periodo):
  //   1 = Fría        : LST < media - 1 std
  //   2 = Templada     : media - 1 std <= LST <= media + 1 std
  //   3 = Cálida (ICU) : media + 1 std < LST <= media + 2 std
  //   4 = Muy cálida   : LST > media + 2 std
  // (Tu programa define ICU como "LST > media + 1 std", es decir,
  //  las clases 3 y 4 combinadas)
  // ------------------------------------------------------------------
  var statsZonales = compuestoLST.reduceRegion({
    reducer: ee.Reducer.mean().combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true}),
    geometry: zonasEstudio.geometry(),
    scale: 30,
    maxPixels: 1e9
  });
  var mediaLST = ee.Number(statsZonales.get('LST_SW_mean'));
  var stdLST = ee.Number(statsZonales.get('LST_SW_stdDev'));

  var umbral1 = mediaLST.subtract(stdLST);          // fría / templada
  var umbral2 = mediaLST.add(stdLST);               // templada / cálida (> 1 std = ICU)
  var umbral3 = mediaLST.add(stdLST.multiply(2));   // cálida / muy cálida

  var zonasICU = ee.Image.constant(1)
    .where(compuestoLST.gt(umbral1), 2)
    .where(compuestoLST.gt(umbral2), 3)
    .where(compuestoLST.gt(umbral3), 4)
    .updateMask(compuestoLST.mask())
    .rename('ICU_clase')
    .clip(zonasEstudio);

  print('Umbrales ICU ' + nombrePeriodo + ' (°C) - media:', mediaLST,
    ' | std:', stdLST, ' | fría<', umbral1, ' | templada<=', umbral2, ' | cálida<=', umbral3);

  Map.addLayer(zonasICU, {min: 1, max: 4, palette: ['2166ac', 'f7f7f7', 'fdae61', 'd73027']},
    'ICU_clases_' + nombrePeriodo);

  Export.image.toDrive({
    image: zonasICU,
    description: 'ICU_clases_' + nombrePeriodo,
    folder: 'GEE_IslasCalor',
    fileNamePrefix: 'ICU_clases_' + nombrePeriodo,
    region: zonasEstudio.geometry(),
    scale: 30,
    crs: 'EPSG:32615',
    maxPixels: 1e9
  });
    // ------------------------------------------------------------------
  // 11. ÁREA (ha) Y % DE CADA CLASE ICU DENTRO DEL ÁREA DE ESTUDIO
  // ------------------------------------------------------------------
  var areaHa = ee.Image.pixelArea().divide(10000); // m2 -> hectáreas

  var agrupado = areaHa.addBands(zonasICU)
    .reduceRegion({
      reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'clase'}),
      geometry: zonasEstudio.geometry(),
      scale: 30,
      maxPixels: 1e9
    });

  var listaGrupos = ee.List(agrupado.get('groups'));
  var areaTotalHa = listaGrupos.map(function(g) {
    return ee.Dictionary(g).getNumber('sum');
  }).reduce(ee.Reducer.sum());

  var nombresClase = ee.Dictionary({'1': 'Fría', '2': 'Templada', '3': 'Cálida (ICU)', '4': 'Muy cálida (núcleo ICU)'});

  var resumenClases = listaGrupos.map(function(g) {
    g = ee.Dictionary(g);
    var clase = g.getNumber('clase');
    var areaClaseHa = g.getNumber('sum');
    var porcentaje = areaClaseHa.divide(ee.Number(areaTotalHa)).multiply(100);
    return ee.Dictionary({
      clase: clase,
      nombre: nombresClase.get(clase.format('%d')),
      area_ha: areaClaseHa,
      porcentaje: porcentaje
    });
  });

  print('Área (ha) y % por clase ICU - ' + nombrePeriodo + ':', resumenClases);
});
