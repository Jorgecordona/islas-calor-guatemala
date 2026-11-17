/**************************************************************************
 * 03_serie_mensual_zonas.js
 * Serie mensual 2024–2025 por zona (1, 4, 9, 10) y para el área total:
 *   - LST (Landsat 8/9, Split-Window, Du et al., 2015)
 *   - NDVI y NDBI (Sentinel-2 L2A, máscara SCL + s2cloudless)
 *   - NO2 troposférico y CO (Sentinel-5P/TROPOMI L3)
 *   - PM2.5 (MERRA-2 M2T1NXAER, colección NASA/GSFC/MERRA/aer/2)
 *   - Contexto ERA5-Land: temperatura del aire, viento y humedad relativa
 * Salida: una tabla CSV (zona × mes) a Google Drive.
 * Actividades 12 y 13 del cronograma; insumo de la actividad 14.
 **************************************************************************/

// ======================= 1. PARÁMETROS =======================
// Subí Zonas14910_Guate_city.shp como asset y pegá aquí la ruta
var zonas = ee.FeatureCollection('projects/icu-guate-ciudad/assets/zonas14910_guate_city');  // reemplazar por tu propio asset
var CAMPO_ZONA = 'ZONA';            // nombre del campo con el número de zona (revisalo en la tabla de atributos)
var INICIO = '2024-01-01';
var FIN    = '2026-01-01';          // exclusivo
var ESCALA = 30;                    // m; las variables gruesas se promedian por área de zona
var CARPETA_DRIVE = 'GEE_ICU';

var aoi = zonas.geometry().dissolve();

// Regiones de salida: cada zona + el área total
var regiones = zonas.map(function (f) {
  return ee.Feature(f.geometry(), {zona: ee.Algorithms.String(f.get(CAMPO_ZONA))});
}).merge(ee.FeatureCollection([ee.Feature(aoi, {zona: 'Total'})]));

// ======================= 2. LST: Split-Window (Du et al., 2015) =======================
// Mismo algoritmo que gee/01_lst_split_window_icu.js (compuestos estacionales):
// bandas TOA B10 y B11, emisividad por banda con umbral de NDVI (Skoković et al., 2014)
// y vapor de agua dinámico de NCEP/NCAR (pr_wtr).
function enmascararNubes(imagen) {
  var qa = imagen.select('QA_PIXEL');
  var mascara = qa.bitwiseAnd(1 << 1).eq(0)   // nube dilatada
    .and(qa.bitwiseAnd(1 << 2).eq(0))         // cirros
    .and(qa.bitwiseAnd(1 << 3).eq(0))         // nube
    .and(qa.bitwiseAnd(1 << 4).eq(0));        // sombra de nube
  return imagen.updateMask(mascara);
}

function obtenerCoeficientes(cwv) {
  var b0 = ee.Image.constant(-0.41165), b1 = ee.Image.constant(1.00522);
  var b2 = ee.Image.constant(0.14543), b3 = ee.Image.constant(-0.27297);
  var b4 = ee.Image.constant(4.06655), b5 = ee.Image.constant(-6.92512);
  var b6 = ee.Image.constant(-18.27461), b7 = ee.Image.constant(0.24468);
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

function obtenerCWV(fechaImagen) {
  var fecha = ee.Date(fechaImagen);
  var cwv = ee.ImageCollection('NCEP_RE/surface_wv')
    .filterDate(fecha, fecha.advance(1, 'day')).mean().select('pr_wtr').unmask(35);
  return cwv.divide(10).rename('CWV');  // kg/m² → g/cm²
}

function calcularEmisividad(ndvi, epsSuelo, epsVeg) {
  var Pv = ndvi.subtract(0.2).divide(0.3).pow(2).clamp(0, 1);
  var F = 0.55;
  var mezcla = ee.Image.constant(epsVeg).multiply(Pv)
    .add(ee.Image.constant(epsSuelo).multiply(ee.Image.constant(1).subtract(Pv)))
    .add(ee.Image.constant(1).subtract(epsSuelo).multiply(epsVeg).multiply(F)
      .multiply(ee.Image.constant(1).subtract(Pv)));
  return mezcla.where(ndvi.lt(0.2), epsSuelo).where(ndvi.gt(0.5), epsVeg);
}

function calcularLST_SW(imagen) {
  var BTi = imagen.select('B10');
  var BTj = imagen.select('B11');
  var sr = imagen.select(['SR_B4', 'SR_B5']).multiply(0.0000275).add(-0.2);
  var ndvi = sr.normalizedDifference(['SR_B5', 'SR_B4']);
  var eps10 = calcularEmisividad(ndvi, 0.971, 0.987);
  var eps11 = calcularEmisividad(ndvi, 0.977, 0.989);
  var epsMedia = eps10.add(eps11).divide(2);
  var deltaEps = eps10.subtract(eps11);
  var coef = obtenerCoeficientes(obtenerCWV(imagen.date()));
  var promedioT = BTi.add(BTj).divide(2);
  var diferenciaT = BTi.subtract(BTj);
  var term1 = coef.b1.add(coef.b2.multiply(ee.Image.constant(1).subtract(epsMedia)).divide(epsMedia))
    .add(coef.b3.multiply(deltaEps).divide(epsMedia.pow(2))).multiply(promedioT);
  var term2 = coef.b4.add(coef.b5.multiply(ee.Image.constant(1).subtract(epsMedia)).divide(epsMedia))
    .add(coef.b6.multiply(deltaEps).divide(epsMedia.pow(2))).multiply(diferenciaT.divide(2));
  var term3 = coef.b7.multiply(diferenciaT.pow(2));
  return coef.b0.add(term1).add(term2).add(term3).subtract(273.15).rename('LST')
    .copyProperties(imagen, ['system:time_start']);
}

var toa = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA')
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_TOA'))
  .filterBounds(aoi).filterDate(INICIO, FIN).map(enmascararNubes);
var l2 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
  .filterBounds(aoi).filterDate(INICIO, FIN).select(['SR_B4', 'SR_B5']);
var unidasL = ee.Join.inner().apply(toa, l2,
  ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'}));
var coleccionLST = ee.ImageCollection(unidasL.map(function (par) {
  return ee.Image.cat(par.get('primary'), par.get('secondary'))
    .copyProperties(ee.Image(par.get('primary')), ['system:time_start']);
})).map(calcularLST_SW);

// ======================= 3. SENTINEL-2: NDVI y NDBI =======================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(INICIO, FIN)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 60));

var s2nubes = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
  .filterBounds(aoi)
  .filterDate(INICIO, FIN);

var s2unida = ee.ImageCollection(ee.Join.saveFirst('prob_nubes').apply({
  primary: s2,
  secondary: s2nubes,
  condition: ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'})
}));

function indicesS2(img) {
  var scl = img.select('SCL');
  // Excluye: 1 saturado, 3 sombra de nube, 8-9 nube media/alta, 10 cirros
  var okScl = scl.neq(1).and(scl.neq(3)).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  var prob = ee.Image(img.get('prob_nubes')).select('probability');
  var mascara = okScl.and(prob.lt(40));
  var sr = img.select(['B4', 'B8', 'B11']).multiply(0.0001);
  var ndvi = sr.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndbi = sr.normalizedDifference(['B11', 'B8']).rename('NDBI');
  return ndvi.addBands(ndbi).updateMask(mascara)
    .copyProperties(img, ['system:time_start']);
}
var s2indices = s2unida.map(indicesS2);

// ======================= 4. SENTINEL-5P: NO2 y CO =======================
// Unidades de salida: NO2 en µmol/m², CO en mmol/m²
var no2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2')
  .filterBounds(aoi).filterDate(INICIO, FIN)
  .select('tropospheric_NO2_column_number_density')
  .map(function (i) {
    return i.multiply(1e6).rename('NO2').copyProperties(i, ['system:time_start']);
  });

var co = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_CO')
  .filterBounds(aoi).filterDate(INICIO, FIN)
  .select('CO_column_number_density')
  .map(function (i) {
    return i.multiply(1e3).rename('CO').copyProperties(i, ['system:time_start']);
  });

// ======================= 5. MERRA-2: PM2.5 =======================
// PM2.5 = DUST2.5 + SS2.5 + BC + 1.4·OC + 1.375·SO4  (Buchard et al., 2016)
// kg/m³ → µg/m³ (× 1e9). Datos horarios; se promedian por mes.
var pm25 = ee.ImageCollection('NASA/GSFC/MERRA/aer/2')
  .filterDate(INICIO, FIN)
  .map(function (i) {
    return i.expression('(DU + SS + BC + 1.4 * OC + 1.375 * SO4) * 1e9', {
      DU: i.select('DUSMASS25'),
      SS: i.select('SSSMASS25'),
      BC: i.select('BCSMASS'),
      OC: i.select('OCSMASS'),
      SO4: i.select('SO4SMASS')
    }).rename('PM25').copyProperties(i, ['system:time_start']);
  });

// ======================= 6. ERA5-Land (contexto) =======================
var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
  .filterDate(INICIO, FIN)
  .map(function (i) {
    var t = i.select('temperature_2m').subtract(273.15);
    var td = i.select('dewpoint_temperature_2m').subtract(273.15);
    // Humedad relativa por la fórmula de Magnus (%)
    var hr = td.expression('100 * exp(17.625*Td/(243.04+Td)) / exp(17.625*T/(243.04+T))',
      {T: t, Td: td});
    var u = i.select('u_component_of_wind_10m');
    var v = i.select('v_component_of_wind_10m');
    var viento = u.pow(2).add(v.pow(2)).sqrt();
    return t.rename('T_aire').addBands(hr.rename('HR')).addBands(viento.rename('Viento'))
      .copyProperties(i, ['system:time_start']);
  });

// ======================= 7. COMPUESTOS MENSUALES =======================
// Si un mes no tiene imágenes, se devuelve una banda vacía (sale como celda vacía en el CSV)
function compuesto(col, banda, ini, fin, reductor) {
  var c = col.filterDate(ini, fin).select(banda);
  var img = ee.Algorithms.If(
    c.size().gt(0),
    c.reduce(reductor).rename(banda),
    ee.Image.constant(0).rename(banda).updateMask(0)
  );
  return ee.Image(img).toFloat();
}

function periodo(d) {
  var m = d.get('month');
  var anio = d.get('year').format('%d');
  var est = ee.String(ee.Algorithms.If(m.lte(4), 'Verano',
            ee.Algorithms.If(m.lte(10), 'Invierno', 'Transicion')));
  return est.cat(anio);
}

var nMeses = ee.Date(FIN).difference(ee.Date(INICIO), 'month').round();
var meses = ee.List.sequence(0, nMeses.subtract(1)).map(function (k) {
  return ee.Date(INICIO).advance(k, 'month');
});

var hayLST = (coleccionLST !== null);

var tabla = ee.FeatureCollection(meses.map(function (d) {
  d = ee.Date(d);
  var ini = d;
  var fin = d.advance(1, 'month');

  var bandas = [
    compuesto(s2indices, 'NDVI', ini, fin, ee.Reducer.median()),
    compuesto(s2indices, 'NDBI', ini, fin, ee.Reducer.median()),
    compuesto(no2, 'NO2', ini, fin, ee.Reducer.mean()),
    compuesto(co, 'CO', ini, fin, ee.Reducer.mean()),
    compuesto(pm25, 'PM25', ini, fin, ee.Reducer.mean()),
    compuesto(era5, 'T_aire', ini, fin, ee.Reducer.mean()),
    compuesto(era5, 'HR', ini, fin, ee.Reducer.mean()),
    compuesto(era5, 'Viento', ini, fin, ee.Reducer.mean())
  ];
  var nLandsat = 0;
  if (hayLST) {
    var lstMes = compuesto(coleccionLST, 'LST', ini, fin, ee.Reducer.median());
    bandas.unshift(lstMes.mask().unmask(0).rename('pct_LST_valido'));
    bandas.unshift(lstMes);
    nLandsat = coleccionLST.filterDate(ini, fin).size();
  }
  var img = ee.Image.cat(bandas);

  var red = img.reduceRegions({
    collection: regiones,
    reducer: ee.Reducer.mean(),
    scale: ESCALA,
    tileScale: 4
  });

  var props = {
    anio: d.get('year'),
    mes: d.get('month'),
    fecha: d.format('YYYY-MM'),
    periodo: periodo(d),
    n_escenas_landsat: nLandsat,
    n_escenas_s2: s2indices.filterDate(ini, fin).size()
  };
  return red.map(function (f) { return f.set(props).setGeometry(null); });
})).flatten();

print('Vista previa (primeras 10 filas):', tabla.limit(10));

var columnas = ['zona', 'fecha', 'anio', 'mes', 'periodo', 'n_escenas_landsat', 'n_escenas_s2',
  'LST', 'pct_LST_valido', 'NDVI', 'NDBI', 'NO2', 'CO', 'PM25', 'T_aire', 'HR', 'Viento'];

Export.table.toDrive({
  collection: tabla,
  description: 'serie_mensual_zonas_2024_2025',
  folder: CARPETA_DRIVE,
  fileFormat: 'CSV',
  selectors: columnas
});

// Control visual rápido
Map.centerObject(aoi, 13);
Map.addLayer(zonas, {color: 'black'}, 'Zonas 1, 4, 9 y 10');
