// load feature collection
var area = ee.FeatureCollection('users/facu_ruarte/tesis/area_estudio_v1SJ_181225');
Map.centerObject(area, 10);

// --------------------------------------------------------------------------- //
// 1. Funciones (Máscara y Escalado)
// --------------------------------------------------------------------------- //

// Máscara basada en SCL (Scene Classification map)
// Clases SCL: 3 (Sombra), 8 (Nube media), 9 (Nube alta), 10 (Cirrus), 11 (Nieve)
// Mantenemos: 4 (Veg), 5 (Suelo), 6 (Agua), 1 (Saturado), 2 (Oscuro), 7 (No clasificado)
function maskS2clouds(image) {
  var scl = image.select('SCL');
  // Seleccionamos lo que NO queremos (Nubes, Sombras, Nieve si aplica)
  var unwanted = scl.eq(3).or(scl.eq(8)).or(scl.eq(9)).or(scl.eq(10)).or(scl.eq(11)); 
  return image.updateMask(unwanted.not());
}

// Función de escalado
function scaling_SR(image){
  // Seleccionamos bandas ópticas y escalamos. Copiamos propiedades de tiempo.
  return image.select('B.*').multiply(0.0001)
              .copyProperties(image, ['system:time_start']);
}

// Funciones de Índices
var addIndicesS2 = function(img) {
    // --- Agua ---
    var ndwi = img.normalizedDifference(['B3','B8']).rename('NDWI');
    var ndmi = img.normalizedDifference(['B8','B11']).rename('NDMI');
    var msi = img.expression('SWIR1 / NIR', {
        'NIR' : img.select('B8'),
        'SWIR1' : img.select('B11')
    }).rename('MSI');
   
    // --- Vegetación ---
    var ndvi = img.normalizedDifference(['B8','B4']).rename('NDVI');
    
    var evi = img.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
        'NIR':img.select('B8'),
        'RED':img.select('B4'),
        'BLUE': img.select('B2')
    }).rename('EVI');

    var evi2 = img.expression(
        '2.4 * (NIR - RED) / (NIR + RED + 1.0)', {
        'NIR': img.select('B8'),
        'RED': img.select('B4')
    }).rename('EVI2');
     
    var ndre1 = img.normalizedDifference(['B8','B5']).rename('NDRE1');

    // --- Disturbios / Suelo ---
    var bnri = img.normalizedDifference(['B8','B12']).rename('BNRI');

    var bsi = img.expression(
        '((Red + SWIR1) - (NIR + Blue)) / ((Red + SWIR1) + (NIR + Blue))', {
        'NIR' : img.select('B8'),
        'SWIR1' : img.select('B11'),
        'Red' : img.select('B4'),
        'Blue' : img.select('B2')
    }).rename('BSI');

    return img.addBands([ndwi, ndmi, msi, ndvi, evi, evi2, ndre1, bnri, bsi]);
};  

// --------------------------------------------------------------------------- //
// 2. Colecciones Sentinel-2 
// --------------------------------------------------------------------------- //

var s2_base = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                  .filterBounds(area)
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30)) // Subí un poco el umbral para permitir que la máscara SCL trabaje
                  .map(maskS2clouds) // 1. Enmascarar
                  .map(scaling_SR)   // 2. Escalar
                  .map(addIndicesS2); // 3. Índices

// Definir temporadas
var S2_summer = s2_base.filterDate('2022-01-01', '2022-03-31');
var S2_winter = s2_base.filterDate('2021-07-01', '2021-09-30');

print("# Escenas Verano: ", S2_summer.size());
print("# Escenas Invierno: ", S2_winter.size());

// --------------------------------------------------------------------------- //
// 3. Reducción (Composite) y Clipping
// --------------------------------------------------------------------------- //

// Seleccionamos bandas de interés usando Regex (Bandas B y los índices calculados)
// Nota: B[2-4,8] selecciona B2, B3, B4, B8.
var bandsToSelect = ['B[2-4,8]', 'NDVI', 'NDWI']; 

var S2_sum_med = S2_summer.select(bandsToSelect).median().clip(area);
var S2_win_med = S2_winter.select(bandsToSelect).median().clip(area);

// Visualización rápida
var vizNDVI = {min: 0, max: 0.8, palette: ['red', 'yellow', 'green']};
Map.addLayer(S2_sum_med.select('NDVI'), vizNDVI, "Summer NDVI", false);
Map.addLayer(S2_win_med.select('NDVI'), vizNDVI, "Winter NDVI", false);
Map.addLayer(area, {'color': 'blue'}, "Study Area", false); 

// --------------------------------------------------------------------------- //
// 4. Diferencias y Stack Final
// --------------------------------------------------------------------------- //

// Diferencias (Verano - Invierno)
// Importante: Al hacer median() las bandas no cambian de nombre automáticamente en GEE
// a menos que uses reduce(ee.Reducer.median()).
// Como usé .median() directo arriba, los nombres siguen siendo 'NDVI', 'NDWI', etc.
var dif_NDVI = S2_sum_med.select('NDVI').subtract(S2_win_med.select('NDVI')).rename('Dif_NDVI');
var dif_NDWI = S2_sum_med.select('NDWI').subtract(S2_win_med.select('NDWI')).rename('Dif_NDWI');

// Renombrado Dinámico (Más seguro que hacerlo manual)
function addSuffix(image, suffix) {
  var names = image.bandNames();
  var newNames = names.map(function(n) { return ee.String(n).cat(suffix) });
  return image.rename(newNames);
}

var summer_final = addSuffix(S2_sum_med, '_sum');
var winter_final = addSuffix(S2_win_med, '_win');

print('Bandas Verano final:', summer_final.bandNames());

// Stack Final
var S2_stack = summer_final.addBands(winter_final)
                           .addBands(dif_NDVI)
                           .addBands(dif_NDWI);

print("Stack Final para Clasificación: ", S2_stack);

// Visualización de Diferencia
Map.addLayer(dif_NDVI, {min: -0.3, max: 0.3, palette: ['red', 'white', 'green']}, 'Diferencia NDVI (Sum-Win)');

// --------------------------------------------------------------------------- //
// 5. Exportar a Google Drive
// --------------------------------------------------------------------------- //

Export.image.toDrive({
  image: S2_stack.toFloat(),        
  
  description: 'Export_Sentinel2_Stack_Tesis',
  folder: 'TESIS_GEE',
  fileNamePrefix: 'S2_Seasonal_Composite_SanJuan',
  region: area,
  scale: 10,
  maxPixels: 1e13,
  crs: 'EPSG:4326',
  fileFormat: 'GeoTIFF'
});
