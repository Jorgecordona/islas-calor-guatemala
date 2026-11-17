// ========================================================================
// SERIE POR ESCENA 2019-2025: LST (Split-Window) + contaminación + ERA5-Land
// a la hora exacta del paso de Landsat.
// Salida (Drive/GEE_ICU): una fila por escena Landsat 8/9, con la cobertura
// de LST para que puedas filtrar después (0.5, 0.8, ...) sin volver a GEE.
// Archivos: Escenas_LST_ERA5_<año>.csv  (7 tareas) o un solo archivo (ver 0.)
// ========================================================================

// --------------------------------------------------------------------
// 0. PARÁMETROS
// --------------------------------------------------------------------
var zonasEstudio = ee.FeatureCollection('projects/icu-guate-ciudad/assets/Zonas14910_Guate_city');
var geom = zonasEstudio.geometry();
var centro = geom.centroid(1);
Map.centerObject(zonasEstudio, 12);

var ANIO_INI = 2019;
var ANIO_FIN = 2025;               // inclusivo
var DIVIDIR_POR_ANIO = true;       // true: 7 tareas (más seguro); false: 1 sola tarea
var COBERTURA_EXPORT = 0.3;        // se exportan escenas con cobertura >= 0.3 (filtra más fino en Python)
var SIN_DATO = -9999;
var CARPETA = 'GEE_ICU';

// --------------------------------------------------------------------
// 1. TUS FUNCIONES DE LST (sin cambios, Split-Window Du et al. 2015)
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

function obtenerCWV(fechaImagen) {
  var fecha = ee.Date(fechaImagen);
  var coleccionCWV = ee.ImageCollection('NCEP_RE/surface_wv')
    .filterDate(fecha, fecha.advance(1, 'day'));
  var cwv_kgm2 = coleccionCWV.mean().select('pr_wtr').unmask(35);
  return cwv_kgm2.divide(10).rename('CWV'); // g/cm2
}

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

function unirColecciones(toa, l2) {
  var filtro = ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'});
  var unidas = ee.Join.inner().apply(toa, l2, filtro);
  return ee.ImageCollection(unidas.map(function(par) {
    return ee.Image.cat(par.get('primary'), par.get('secondary'));
  }));
}

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
// 2. CONTAMINACIÓN (mismas funciones del script de 24 meses)
// --------------------------------------------------------------------
function valorOSinDato(v) {
  return ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(v, null), SIN_DATO, v));
}

var BANDAS_PM = ['DUSMASS25', 'SSSMASS25', 'BCSMASS', 'OCSMASS', 'SO4SMASS'];
var merra = ee.ImageCollection('NASA/GSFC/MERRA/aer/2').select(BANDAS_PM);

// PM2.5 = DU25 + SS25 + BC + 1.4*OC + 1.375*SO4 (µg/m3). Celda MERRA-2 sobre el centro del área.
function pm25Medio(col) {
  var m = col.mean();
  var pm = m.select('DUSMASS25').add(m.select('SSSMASS25')).add(m.select('BCSMASS'))
    .add(m.select('OCSMASS').multiply(1.4))
    .add(m.select('SO4SMASS').multiply(1.375))
    .multiply(1e9).rename('PM25');
  var v = pm.reduceRegion({reducer: ee.Reducer.mean(), geometry: centro, scale: 55000, maxPixels: 1e9}).get('PM25');
  return ee.Number(ee.Algorithms.If(col.size().gt(0), valorOSinDato(v), SIN_DATO));
}

function no2Col(ini, fin) {
  return ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2')
    .filterBounds(geom).filterDate(ini, fin)
    .map(function (i) { return i.updateMask(i.select('cloud_fraction').lt(0.3)); })
    .select('tropospheric_NO2_column_number_density');
}
function aaiCol(ini, fin) {
  return ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_AER_AI')
    .filterBounds(geom).filterDate(ini, fin)
    .select('absorbing_aerosol_index');
}
function s5pMedio(col, banda, factor) {
  var v = col.mean().reduceRegion({reducer: ee.Reducer.mean(), geometry: geom, scale: 1113.2, maxPixels: 1e9}).get(banda);
  var vf = ee.Algorithms.If(ee.Algorithms.IsEqual(v, null), null, ee.Number(v).multiply(factor));
  return ee.Number(ee.Algorithms.If(col.size().gt(0), valorOSinDato(vf), SIN_DATO));
}

// --------------------------------------------------------------------
// 3. ERA5-Land A LA HORA DEL PASO (interpolación lineal entre horas)
// --------------------------------------------------------------------
var era5h = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY');
var era5d = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR');
var BANDAS_H = ['temperature_2m', 'dewpoint_temperature_2m', 'skin_temperature',
                'u_component_of_wind_10m', 'v_component_of_wind_10m',
                'volumetric_soil_water_layer_1', 'surface_solar_radiation_downwards_hourly'];

// t: ee.Date del paso del satélite. Devuelve un ee.Dictionary con las variables del clima.
function era5Paso(t) {
  var ms = ee.Number(t.millis());
  var h0 = ee.Date(ms.divide(3600000).floor().multiply(3600000));
  var h1 = h0.advance(1, 'hour');
  var w = ms.subtract(h0.millis()).divide(3600000);             // fracción de hora transcurrida

  var i0 = ee.Image(era5h.filterDate(h0, h1).select(BANDAS_H).first());
  var i1 = ee.Image(era5h.filterDate(h1, h1.advance(1, 'hour')).select(BANDAS_H).first());
  var img = i0.multiply(ee.Number(1).subtract(w)).add(i1.multiply(w));
  var d = ee.Dictionary(img.reduceRegion({reducer: ee.Reducer.mean(), geometry: centro, scale: 11132, maxPixels: 1e9}));

  var T = ee.Number(d.get('temperature_2m')).subtract(273.15);
  var Td = ee.Number(d.get('dewpoint_temperature_2m')).subtract(273.15);
  var HR = Td.multiply(17.625).divide(Td.add(243.04)).exp()
    .divide(T.multiply(17.625).divide(T.add(243.04)).exp()).multiply(100);
  var u = ee.Number(d.get('u_component_of_wind_10m'));
  var v = ee.Number(d.get('v_component_of_wind_10m'));
  var viento = u.pow(2).add(v.pow(2)).sqrt();

  // Lluvia previa: 24 h antes del paso (horaria) y 7 días UTC completos previos (diaria)
  var pp24 = era5h.filterDate(t.advance(-24, 'hour'), t).select('total_precipitation_hourly').sum()
    .reduceRegion({reducer: ee.Reducer.mean(), geometry: centro, scale: 11132, maxPixels: 1e9})
    .get('total_precipitation_hourly');
  var dia = ee.Date(ms.divide(86400000).floor().multiply(86400000));
  var pp7 = era5d.filterDate(dia.advance(-7, 'day'), dia).select('total_precipitation_sum').sum()
    .reduceRegion({reducer: ee.Reducer.mean(), geometry: centro, scale: 11132, maxPixels: 1e9})
    .get('total_precipitation_sum');

  return ee.Dictionary({
    T_aire_C: T,
    T_rocio_C: Td,
    HR_pct: HR,
    VIENTO_ms: viento,
    RAD_Wm2: ee.Number(d.get('surface_solar_radiation_downwards_hourly')).divide(3600),
    T_skin_ERA5_C: ee.Number(d.get('skin_temperature')).subtract(273.15),
    HUM_SUELO: d.get('volumetric_soil_water_layer_1'),
    PP_24h_mm: ee.Number(pp24).multiply(1000),
    PP_7d_mm: ee.Number(pp7).multiply(1000)
  });
}

// --------------------------------------------------------------------
// 4. PIPELINE POR PERIODO
// --------------------------------------------------------------------
var pixTotal = ee.Number(ee.Image.constant(1).clip(geom)
  .reduceRegion({reducer: ee.Reducer.count(), geometry: geom, scale: 30, maxPixels: 1e9})
  .get('constant'));

var SELECTORES = ['fecha', 'mes', 'anio', 'satelite', 'hora_UTC', 'LST_media', 'LST_sd', 'cobertura',
  'PM25_dia', 'NO2_dia_umol_m2', 'AAI_dia',
  'T_aire_C', 'T_rocio_C', 'HR_pct', 'VIENTO_ms', 'RAD_Wm2', 'T_skin_ERA5_C',
  'HUM_SUELO', 'PP_24h_mm', 'PP_7d_mm'];

function procesarPeriodo(ini, fin, etiqueta) {
  var toa = ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA')
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_TOA'))
    .filterBounds(zonasEstudio).filterDate(ini, fin)
    .map(enmascararNubes);
  var l2 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(zonasEstudio).filterDate(ini, fin)
    .select(['SR_B4', 'SR_B5']);
  var conLST = unirColecciones(toa, l2).map(calcularLST_SW);

  // Etapa 1: solo LST y cobertura (barato). Se descartan escenas casi totalmente nubladas.
  var etapa1 = ee.FeatureCollection(conLST.map(function (img) {
    img = ee.Image(img);
    var st = ee.Dictionary(img.select('LST_SW').reduceRegion({
      reducer: ee.Reducer.mean()
        .combine({reducer2: ee.Reducer.count(), sharedInputs: true})
        .combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true}),
      geometry: geom, scale: 30, maxPixels: 1e9
    }));
    return ee.Feature(null, {
      t_ms: img.get('system:time_start'),
      satelite: img.get('SPACECRAFT_ID'),
      LST_media: st.get('LST_SW_mean'),
      LST_sd: st.get('LST_SW_stdDev'),
      cobertura: ee.Number(st.get('LST_SW_count', 0)).divide(pixTotal)
    });
  })).filter(ee.Filter.gte('cobertura', COBERTURA_EXPORT));

  // Etapa 2: contaminación y ERA5-Land solo para las escenas que quedaron
  var etapa2 = etapa1.map(function (f) {
    var t = ee.Date(f.get('t_ms'));
    var dia = ee.Date(t.format('YYYY-MM-dd'));
    var colDia = merra.filterDate(dia, dia.advance(1, 'day'));
    var clima = era5Paso(t);
    return f.set({
      fecha: t.format('YYYY-MM-dd'),
      mes: t.format('YYYY-MM'),
      anio: t.get('year'),
      hora_UTC: ee.Number(t.get('hour')).add(ee.Number(t.get('minute')).divide(60)),
      PM25_dia: pm25Medio(colDia),
      NO2_dia_umol_m2: s5pMedio(no2Col(dia, dia.advance(1, 'day')), 'tropospheric_NO2_column_number_density', 1e6),
      AAI_dia: s5pMedio(aaiCol(dia, dia.advance(1, 'day')), 'absorbing_aerosol_index', 1)
    }).set(clima);
  });

  print('Escenas Landsat (TOA) ' + etiqueta + ':', toa.size());
  Export.table.toDrive({
    collection: etapa2,
    description: 'Escenas_LST_ERA5_' + etiqueta,
    folder: CARPETA, fileNamePrefix: 'Escenas_LST_ERA5_' + etiqueta, fileFormat: 'CSV',
    selectors: SELECTORES
  });
}

if (DIVIDIR_POR_ANIO) {
  for (var y = ANIO_INI; y <= ANIO_FIN; y++) {
    procesarPeriodo(y + '-01-01', (y + 1) + '-01-01', String(y));
  }
} else {
  procesarPeriodo(ANIO_INI + '-01-01', (ANIO_FIN + 1) + '-01-01', ANIO_INI + '_' + ANIO_FIN);
}
