// 1. CARGA DE ÁREA DE ESTUDIO
var area = ee.FeatureCollection('users/facu_ruarte/tesis/area_estudio_v1SJ_181225');
Map.centerObject(area, 10);
Map.addLayer(area, {color: 'red'}, 'Área de Estudio', false);

// 2. MODELO DIGITAL DE ELEVACIÓN (DEM)
// Usamos SRTM de 30 metros
var dem = ee.Image("USGS/SRTMGL1_003")
  .clip(area);

// 3. CÁLCULO DE LA PENDIENTE (SLOPE)
// El resultado está en grados por defecto
var slope = ee.Terrain.slope(dem);

// 4. CÁLCULO DE DISTANCIA A RÍOS
// Opción A: Usar red de drenaje global (HydroSHEDS)
// Nota: Puedes reemplazar esto por tu propio Asset de ríos si lo prefieres
var rivers = ee.FeatureCollection("WWF/HydroSHEDS/v1/FreeFlowingRivers")
  .filterBounds(area);

// Creamos una máscara: 1 donde hay río, 0 donde no
// Convertimos los vectores a una imagen ráster de 30m
var riversRaster = rivers
  .style({color: '0000FF', width: 1})
  .mask()
  .select('vis-red'); // Seleccionamos la máscara de dibujo

// Calculamos la distancia euclidiana en metros
// fastDistanceTransform calcula distancia en píxeles, luego multiplicamos por la escala
var distanceToRivers = riversRaster.fastDistanceTransform()
  .sqrt()
  .multiply(ee.Image.pixelArea().sqrt())
  .rename('distancia_rios');

// 5. VISUALIZACIÓN EN EL MAPA
var demVis = {min: 500, max: 3000, palette: ['000000', '444444', 'FFFFFF']};
var slopeVis = {min: 0, max: 30, palette: ['white', 'orange', 'red']};
var distVis = {min: 0, max: 2000, palette: ['white', 'blue', 'darkblue']};

Map.addLayer(dem, demVis, '1. DEM (SRTM 30m)');
Map.addLayer(slope, slopeVis, '2. Pendiente (Grados)');
Map.addLayer(distanceToRivers.clip(area), distVis, '3. Distancia a Ríos (m)');
Map.addLayer(rivers, {color: '0000FF'}, 'Vectores de Ríos (HydroSHEDS)');


// 6. EXPORTACIÓN DE RESULTADOS
Export.image.toDrive({
  image: dem,
  description: 'DEM_SanJuan_GEE',
  scale: 30,
  region: area,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: slope,
  description: 'Slope_SanJuan_GEE',
  scale: 30,
  region: area,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: distanceToRivers.clip(area),
  description: 'Distancia_Rios_GEE',
  scale: 30,
  region: area,
  maxPixels: 1e13
}); 
