// load feature collection (= vector) from personal Asset 
var area = ee.FeatureCollection('users/facu_ruarte/Sur_SanJuan_Export');


Map.centerObject(area, 10);



// throw('Stop in line 22')


// --------------------------------------------------------------------------- //
// Colección Sentinel-2 
// --------------------------------------------------------------------------- //

// Filter and bound S2 collection
var S2_summer = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                    .filterDate('2022-01-01', '2022-03-31')
                    .filterBounds(area)
                    // Pre-filter to get less cloudy granules.
                    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));
                  
print("Sentinel Summer: ", S2_summer);
print("# scenes: ", S2_summer.size());



var S2_winter = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                    .filterDate('2021-07-01', '2021-09-30')
                    .filterBounds(area)
                    // Pre-filter to get less cloudy granules.
                    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));
                  
print("Sentinel Winter: ", S2_winter);
print("# scenes: ", S2_winter.size());

// Write a function for Cloud masking
function maskCloudAndShadowsSR(image) {
    var cloudProb = image.select('MSK_CLDPRB');
    var snowProb = image.select('MSK_SNWPRB');
    var cloud = cloudProb.lt(5);
    var snow = snowProb.lt(5);
    var scl = image.select('SCL'); 
    var shadow = scl.eq(3); // 3 = cloud shadow
    var cirrus = scl.eq(10); // 10 = cirrus
    // Cloud probability less than 5% or cloud shadow classification
    var mask = (cloud.and(snow)).and(cirrus.neq(1)).and(shadow.neq(1));
    return image.updateMask(mask);
  }
  

var S2_summer = S2_summer.map(maskCloudAndShadowsSR);
var S2_winter = S2_winter.map(maskCloudAndShadowsSR);  

// throw('stop in line 55') 
  
  
var viz = {
    bands: ['B4','B3','B2'],
    min: 0,
    max: 0.5
}

Map.addLayer(S2_summer.median(), viz, "summer True Color (432)", false);
Map.addLayer(S2_winter.median(), viz, "winter True Color (432)", false);

//throw('stop in line')


function scaling_SR(image){
    return image.select('B.*').multiply(0.0001)
                .copyProperties(image, ['system:time_start']);
  }
  
var S2summ_scaled = S2_summer.map(scaling_SR); 
print('Summer Scaled S2', S2summ_scaled);

var S2win_scaled = S2_winter.map(scaling_SR); 
print('Winter Scaled S2', S2win_scaled);

  
var RGBvis2 = {
    min: 0.0,
    max: 0.6,
    bands: ['B4', 'B3', 'B2'],
  };

Map.addLayer(S2summ_scaled, RGBvis2, 'Final Summer RGB', false);

// throw('stop')









// ----------------------------------------------------------------------------------- //
//     Spectral Indices                                                                //
// ----------------------------------------------------------------------------------- //


var addIndicesS2 = function(img) {
    /// Spectral indices related to water stress
     
     // NDWI (Normalized Difference Water Index)
     var ndwi = img.normalizedDifference(['B3','B8']).toDouble().rename('NDWI');
 
     // NDMI (Normalized Difference Moisture Index)
     var ndmi = img.normalizedDifference(['B8','B11']).toDouble().rename('NDMI');
 
     // MSI (Moisture Stress Index)
     var msi = img.expression('SWIR1 / NIR', {
     'NIR' : img.select('B8'),
     'SWIR1' : img.select('B11')
     }).rename('MSI');
    
    
    /// Spectral indices related to photosynthesis
    
     // NDVI (Normalized Difference Vegetation Index)
     var ndvi = img.normalizedDifference(['B8','B4']).toDouble().rename('NDVI');
     
     // EVI (Enhance Vegetation Index)
     var evi = img.expression(
         '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
         'NIR':img.select('B8'),
         'RED':img.select('B4'),
         'BLUE': img.select('B2')
       }).rename('EVI');

     var evi2 = img.expression(
      '2.4 * (B08 - B04) / (B08 + B04 + 1.0)', {
        'B08': img.select('B8'),
        'B04': img.select('B4')
      }).rename('EVI2');
     
     // NDRE-1 (Normalized Difference Red Edge 1)
     var ndre1 = img.normalizedDifference(['B8','B5']).toDouble().rename('NDRE1');
 
 
    /// Spectral indices related to forest disturbance 
 
     // NBRI (Normalized Burnt Ratio Index)
     var bnri = img.normalizedDifference(['B8','B12']).toDouble().rename('BNRI');
 
     // BSI (Bare Soil Index)
     var bsi = img.expression('(Red + SWIR1) - (NIR + Blue) / (Red + SWIR1) + (NIR + Blue)', {
         'NIR' : img.select('B8'),
         'SWIR1' : img.select('B11'),
         'Red' : img.select('B4'),
         'Blue' : img.select('B2')
       }).rename('BSI');
 
 
       return img
         .addBands(ndwi)
         .addBands(ndmi)
         .addBands(msi)
         .addBands(ndvi)
         .addBands(evi)
         .addBands(evi2)
         .addBands(ndre1)
         .addBands(bnri)
         .addBands(bsi);
     };  

     
var S2summ_indices = S2summ_scaled.map(addIndicesS2);
print(S2summ_indices);

var S2win_indices = S2win_scaled.map(addIndicesS2);
print(S2win_indices);


             
// Generating a image composite
var bands = ['B[2-4,8]', 'NDVI', 'NDWI'];

// Reducing collection by meadian
var S2_sum_med = S2summ_indices.select(bands).reduce(ee.Reducer.median());
print("Summer composite: ", S2_sum_med.bandNames()); 

var S2_win_med = S2win_indices.select(bands).reduce(ee.Reducer.median());
print("Winter composite: ", S2_win_med.bandNames()); 

// throw('stop in line 146')
// clipping images
var S2_sum_clipped = S2_sum_med.clip(area);

var S2_win_clipped = S2_win_med.clip(area);



Map.addLayer(S2_sum_clipped.select('NDVI_median'), {}, "Summer NDVI", false);
Map.addLayer(S2_win_clipped.select('NDVI_median'), {}, "Winter NDVI", false);

Map.addLayer(S2_sum_clipped.select('NDWI_median'), {}, "Summer NDWI", false);
Map.addLayer(S2_win_clipped.select('NDWI_median'), {}, "Winter NDWI", false);

Map.addLayer(area, {'color': 'blue'}, "Study Area",false); 


// throw ('stop')

// ----------------------------------------------------------------------------------- //
//     Diferencias Summer - Winter                                                     //
// ----------------------------------------------------------------------------------- //


var diferencia_NDVI = S2_sum_clipped.select('NDVI_median').subtract(S2_win_clipped.select('NDVI_median')).rename('Dif_NDVI');
print(diferencia_NDVI.bandNames());                        
                        
// 4. Mostrar el resultado (Opcional: configura una visualización adecuada)
Map.addLayer(diferencia_NDVI, {min: -0.5, max: 0.5, palette: ['IndianRed', 'white', 'ForestGreen']}, 'Diferencia NDVI');



var diferencia_NDWI = S2_sum_clipped.select('NDWI_median').subtract(S2_win_clipped.select('NDWI_median')).rename('Dif_NDWI');
print(diferencia_NDWI.bandNames());
                        
// 4. Mostrar el resultado (Opcional: configura una visualización adecuada)
Map.addLayer(diferencia_NDWI, {min: -0.5, max: 0.5, palette: ['Lavender', 'white', 'MediumSlateBlue']}, 'Diferencia NDWI');

// Throw ('stop')
// Add sufijo de winter y summer
var bandas = ['B2_median', 'B3_median', 'B4_median', 
              'B8_median', 'NDVI_median', 'NDWI_median'];

var summer_comp = S2_sum_clipped.rename(bandas.map(function(b){return ee.String(b).cat('_sum')}));
print('Summer Properties', summer_comp);

var winter_comp = S2_win_clipped.rename(bandas.map(function(b){return ee.String(b).cat('_win')}));
print('Winter Properties', winter_comp);


var S2_stack_summer = summer_comp.addBands(winter_comp)
                                 .addBands(diferencia_NDVI)
                                 .addBands(diferencia_NDWI);
print("Stack Bands: ", S2_stack_summer.bandNames());
/*
// 2. Muestrear la Imagen Apilada (Stack)
// Se toman muestras aleatorias de píxeles para entrenar el clusterer.
var entrenamiento = S2_stack_summer.sample({
  region: area,
  scale: 10, // Utiliza la escala de resolución de Sentinel-2 (o ajusta según tu necesidad)
  numPixels: 5000, // Número de píxeles a muestrear
  seed: 0, // Semilla para la reproducibilidad
  tileScale: 16 // Ajuste para grandes áreas si es necesario
});

// 3. Instanciar y Entrenar el Clusterer K-Means
// Se crea una instancia del clusterer wekaKMeans con el número deseado de clusters (ej. 10).
var n_clusters = 8; 
var clusters = ee.Clusterer.wekaKMeans(n_clusters)
  .train(entrenamiento);

// 4. Aplicar el Clusterer a la Imagen Completa
var resultado_cluster = S2_stack_summer.cluster(clusters);

var clusters = ee.Image('users/facu_ruarte/Classification_RAE/KMeans_Clustering_k8_Asset_Export');
print(clusters); 
// throw ('stop')
Map.addLayer(clusters.randomVisualizer(), {}, 
             //'K-means classified 8 - random colors');



// 1. Asumimos que tu ráster clasificado es 'Kclassified_k8' 
//    y que se generó con 8 clústeres (de 0 a 7).
var K = 8;


// 2. Definir los parámetros de visualización (visParams)
var visParams_custom = {
  // El valor mínimo siempre debe ser 0.
  min: 0,
  // El valor máximo debe ser el número de clústeres menos uno (K - 1).
  max: K - 1, // En este caso, 7 (para clústeres 0, 1, 2, 3, 4, 5, 6, 7)
  
  // La paleta debe ser un array con K colores, en el orden de los clústeres (0 al 7).
  palette: [
    '0000FF',  // Clúster 0: Azul (p. ej., Agua)
    '00FF00',  // Clúster 1: Verde (p. ej., Vegetación Densa)
    'ADFF2F',  // Clúster 2: Verde-amarillo (p. ej., Vegetación Escasa)
    'FFFF00',  // Clúster 3: Amarillo (p. ej., Suelo Desnudo)
    'FF0000',  // Clúster 4: Rojo (p. ej., Urbano Denso)
    '8B4513',  // Clúster 5: Marrón (p. ej., Suelo Agrícola)
    'FFFFFF',  // Clúster 6: Blanco (p. ej., Nubes/Nieve)
    'A9A9A9'   // Clúster 7: Gris (p. ej., Sombra)
  ]
};

// 3. Aplicar la paleta de colores al mapa
Map.addLayer(
  clusters, 
  visParams_custom, 
  'K-means clasificado (8 Clústeres Personalizados)'
);

// Nota: Si quieres mantener tu título original, simplemente reemplaza el randomVisualizer
// Map.addLayer(Kclassified_k8.randomVisualizer(), {}, 'K-means classified 8 - random colors');
// por esto:
// Map.addLayer(Kclassified_k8, visParams_custom, 'K-means classified 8 - colores personalizados');


// throw ('stop')


// Definir parámetros de visualización para las diferentes escenas
var visRGB = {bands: ['B4', 'B3', 'B2'], min: 0.0, max: 0.6, gamma: 1.2};
var visNIR = {bands: ['B8'], min: 0, max: 15, palette: ['black', 'red', 'yellow']};

// Calcular un índice (NDVI) para otra escena

var visNDVI = {min: 0, max: 0.8, palette: ['blue', 'white', 'green']};


var left = ui.Map();
var center = ui.Map();
var right = ui.Map();
ui.root.clear();
ui.root.add(left);
ui.root.add(center);
ui.root.add(right);

// Link maps, so when you drag one map, the other will be moved in sync.
ui.Map.Linker([left, center, right], 'change-bounds');

// Center map to the area of interest
right.centerObject(geometry, 12); //Number indicates the zoom level
center.centerObject(geometry, 12)
left.centerObject(geometry, 12); //Number indicates the zoom level

left.addLayer(diferencia_NDVI, 
   {min: -0.5, max: 0.5, palette: ['red', 'gray', 'green']}, 
    "Diferencia NDVI");
center.addLayer(diferencia_NDWI, 
   {min: -0.5, max: 0.5, palette: ['lightblue', 'gray', 'blue']}, 
    "Combinación: 543");
right.addLayer(clusters.randomVisualizer(), {}, 
             'K-means classified 8 - random colors');


var legend1 = ui.Panel({style: {position: 'bottom-left',padding: '8px 15px'}});
legend1.add(ui.Label({
  value: "Diferencia Estacional NDVI",
  style: {fontWeight: 'bold',fontSize: '18px',margin: '0 0 4px 0',padding: '0px'}
}));


var legend2 = ui.Panel({style: {position: 'bottom-left',padding: '8px 15px'}});
legend2.add(ui.Label({
  value: "Diferencia Estacional NDWI",
  style: {fontWeight: 'bold',fontSize: '18px',margin: '0 0 4px 0',padding: '0px'}
}));


var legend3 = ui.Panel({
  style: {
    position: 'bottom-right',
    padding: '8px 15px'
  }
});

var legend3a = ui.Label({
  value: 'K-means',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});


legend3.add(legend3a);

// Creates the content of the legend
var content = function(color, label) {
      // Create the color boxes
      var box = ui.Label({
        style: {
          backgroundColor: '#' + color,
          // Set box height and width
          padding: '9px',
          margin: '0 0 4px 0'
        }
      });
      // Create the labels
      var labels = ui.Label({
        value: label,
        style: {margin: '0 0 4px 6px'}
      });
      return ui.Panel({
        widgets: [box, labels],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};
//  Set legend colors
var classcolor = ['654321','FFA500','FFFF00', '00FF00', '008000'];
// Set legend labels
var labelName = ['<=0','0 - 0.2','0.2 - 0.4','0.4 - 0.6', '>0.6'];
// Combine legend colou and labels
for (var i = 0; i < 5; i++) {
  legend3.add(content(classcolor[i], labelName[i]));
  }  

left.add(legend1);
center.add(legend2);
right.add(legend3); 
*/

// 1. Cargar el Asset Exportado
var asset_cluster = ee.Image('users/facu_ruarte/Classification_RAE/kmeans_cluster_k8')
  .select('KMeans_Cluster_k8')
  .rename('cluster'); // <--- CORRECCIÓN ADICIONAL: Se renombra la banda para asegurar compatibilidad con el visor.
print('Asset Cargado para Visualización (Bandas):', asset_cluster.bandNames()); // Ahora debería decir ['cluster']

var K = 8;

// 2. Definición de Parámetros de Visualización

// Parámetros para RGB (Bandas B4, B3, B2)
// CORRECCIÓN: Usar los nombres de banda disponibles después del stack/median reduction.
var visRGB = {bands: ['B4_median_sum', 'B3_median_sum', 'B2_median_sum'], min: 0.0, max: 0.3, gamma: 1.2};

// Parámetros para Diferencia NDVI
var visDiffNDVI = {min: -0.5, max: 0.5, palette: ['red', 'Silver', 'green']};

// Parámetros para Diferencia NDWI
var visDiffNDWI = {min: -0.5, max: 0.5, palette: ['lightblue', 'Silver', 'blue']};

// Parámetros para K-Means Custom
var visParams_custom = {
  min: 0,
  max: K - 1, // 7
  palette: [
    'C0C0C0',  // Clúster 0: Silver (Plata)
    'A52A2A',  // Clúster 1: Brown (Marrón)
    '808000',  // Clúster 2: Olive (Oliva)
    'EEE8AA',  // Clúster 3: PaleGoldenrod (Dorado Pálido)
    '228B22',  // Clúster 4: ForestGreen (Verde Bosque)
    '3CB371',  // Clúster 5: MediumSeaGreen (Verde Marino Medio)
    'A9A9A9',  // Clúster 6: DarkGray (Gris Oscuro)
    'A9A9A9'   // Clúster 7: Gris (Sombra)
  ]
};

// --- 3. Configuración de Mapas Vinculados (UI) ---

// Definición de los 4 mapas
var map1 = ui.Map(); // RGB
var map2 = ui.Map(); // Diferencia NDVI
var map3 = ui.Map(); // Diferencia NDWI (Nuevo)
var map4 = ui.Map(); // Cluster K-Means


// Limpiar la raíz y añadir los 4 mapas
ui.root.clear();
ui.root.add(map1);
ui.root.add(map2);
ui.root.add(map3);
ui.root.add(map4);

// --- Añadir Flecha Cardinal (North Arrow) ---
// La barra de escala (scale bar) se incluye por defecto en el mapa.

// Crear la Flecha Cardinal (North Arrow) para map4
var northArrow = ui.Panel({
  widgets: [
    // Triángulo Unicode (Flecha)
    ui.Label('\u25B2', {fontSize: '20px', fontWeight: 'bold', color: 'black'}), 
    // Letra N
    ui.Label('N', {fontSize: '14px', fontWeight: 'bold', color: 'black', margin: '-8px 0 0 0'})
  ],
  layout: ui.Panel.Layout.flow('vertical'),
  style: {
    position: 'top-left', 
    padding: '4px 8px', 
    backgroundColor: 'rgba(255, 255, 255, 0.7)', 
    border: '1px solid #ccc',
    borderRadius: '4px'
  }
});

// Añadir la Flecha Cardinal al mapa del Cluster (map4)
map4.add(northArrow);

// Sincronizar los 4 mapas
ui.Map.Linker([map1, map2, map3, map4], 'change-bounds');

// Centrar los 4 mapas en el área de interés
map1.centerObject(geometry, 12); 
map2.centerObject(geometry, 12);
map3.centerObject(geometry, 12);
map4.centerObject(geometry, 12); 

// **ASIGNACIÓN DE CAPAS (Orden Cuádruple):**

// Map 1 (Izquierda): RGB
map1.addLayer(S2_stack_summer, 
  visRGB, 
  "1. Imagen Satelital RGB"); 
  
// Map 2: Diferencia NDVI
map2.addLayer(diferencia_NDVI, 
  visDiffNDVI, 
  "2. Diferencia NDVI");
  
// Map 3: Diferencia NDWI
map3.addLayer(diferencia_NDWI, 
  visDiffNDWI, 
  "3. Diferencia NDWI");

// Map 4 (Derecha): Cluster K-Means
map4.addLayer(
  asset_cluster, 
  visParams_custom, 
  '4. K-means clasificado (Asset - Custom)'
); 

// --- 4. Configuración de Leyendas (UI) ---

// Función auxiliar para crear ítems de leyenda
var createLegendItem = function(color, label) {
       var box = ui.Label({
         style: {
           backgroundColor: '#' + color,
           padding: '8px',
           margin: '0 0 4px 0'
         }
       });
       var labels = ui.Label({
         value: label,
         style: {margin: '0 0 4px 6px'}
       });
       return ui.Panel({
         widgets: [box, labels],
         layout: ui.Panel.Layout.Flow('horizontal')
       });
};

// --- Leyenda 1: RGB (Solo Título) ---
var legend1 = ui.Panel({style: {position: 'bottom-left',padding: '8px 15px'}});
legend1.add(ui.Label({
  value: "1. Imagen Satelital RGB", 
  style: {fontWeight: 'bold',fontSize: '18px',margin: '0 0 4px 0',padding: '0px'}
}));


// --- Leyenda 2: Diferencia NDVI ---
var legend2 = ui.Panel({style: {position: 'bottom-left',padding: '8px 15px'}});
legend2.add(ui.Label({
  value: "2. Diferencia Estacional NDVI", 
  style: {fontWeight: 'bold',fontSize: '18px',margin: '0 0 4px 0',padding: '0px'}
}));
legend2.add(createLegendItem('FF0000', 'Menor NDVI (Pérdida de vegetación)'));
legend2.add(createLegendItem('808080', 'Cambio Nulo (Gris)'));
legend2.add(createLegendItem('008000', 'Mayor NDVI (Aumento de vegetación)'));


// --- Leyenda 3: Diferencia NDWI ---
var legend3 = ui.Panel({style: {position: 'bottom-left',padding: '8px 15px'}});
legend3.add(ui.Label({
  value: "3. Diferencia Estacional NDWI", 
  style: {fontWeight: 'bold',fontSize: '18px',margin: '0 0 4px 0',padding: '0px'}
}));
legend3.add(createLegendItem('ADD8E6', 'Menor NDWI (Pérdida de agua)'));
legend3.add(createLegendItem('808080', 'Cambio Nulo (Gris)'));
legend3.add(createLegendItem('0000FF', 'Mayor NDWI (Aumento de agua)'));


// --- Leyenda 4: K-Means Cluster (Simplificada) ---
var legend4 = ui.Panel({
  style: {
    position: 'bottom-right',
    padding: '8px 15px'
  }
});

legend4.add(ui.Label({
  value: '4. K-means Acople de Vegetación',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
}));

// Definición de las CLASES DE INTERÉS (Solo 4 de 8 clústeres)
var simplifiedClasses = [
  // Clúster 5 (Marrón: Suelo Agrícola)
  { cluster: 1, color: 'A52A2A', label: 'Vegetación No asociada' }, 
  // Clúster 2 (Verde-amarillo: Veg. Escasa)
  { cluster: 5, color: '3CB371', label: 'Vegetación Parcialmente Acoplada' },
  // Clúster 1 (Verde: Vegetación Densa)
  { cluster: 4, color: '228B22', label: 'Vegetación Acoplada' },
  // Clúster 0 (Azul: Agua) o Clúster 3 (Amarillo: Suelo Desnudo)
  { cluster: 3, color: 'EEE8AA', label: 'Vegetación No Acoplada' } // Asumo el clúster 3 (Suelo Desnudo) como "No Acoplada"
];

// Revisa el índice del color en visParams_custom.palette para verificar la asociación:
// 0: 'C0C0C0' (Silver)
// 1: 'A52A2A' (Brown)
// 2: '808000' (Olive)
// 3: 'EEE8AA' (PaleGoldenrod)
// 4: '228B22' (ForestGreen)
// 5: '3CB371' (MediumSeaGreen)
// 6: 'A9A9A9' (DarkGray)
// 7: 'A9A9A9' (Gris)
// Bucle simplificado para mostrar solo las 4 clases definidas
simplifiedClasses.forEach(function(item) {
  // Nota: Creamos el item usando el color y la etiqueta definidos en el array simplifiedClasses
  legend4.add(createLegendItem(item.color, item.label));
});

// Añadir leyendas a los mapas correspondientes
map1.add(legend1);
map2.add(legend2);
map3.add(legend3);
map4.add(legend4);
