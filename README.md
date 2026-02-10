# Mapeo de Bosques Nativos y Freatófitas en San Juan: Integración Bio-Topográfica 🌵💧
Este repositorio contiene el ecosistema de procesamiento y modelado desarrollado para mi tesis de grado en la Universidad Nacional de San Juan (UNSJ). El proyecto se centra en la identificación de vegetación freatófita (*Neltuma flexuosa*) mediante la fusión de datos multiespectrales y variables geomorfológicas.

🔬 Contexto de la Investigación

La clasificación de bosques nativos en zonas áridas requiere ir más allá de la reflectancia satelital. Este estudio integra la dinámica hídrica del terreno para diferenciar bosques con acceso a napas freáticas. Se utiliza un enfoque de Machine Learning robusto para generar mapas de probabilidad de presencia con alta resolución espacial.

🛠️ Stack Tecnológico e Idiomas

- R (Lenguaje principal): Procesamiento de datos espaciales y modelado estadístico avanzado.
- Librerías Clave: terra (Rasters), sf (Vectores), ranger (Random Forest), fastshap & shapviz (Interpretabilidad).
- Google Earth Engine (GEE): Preprocesamiento y obtención de métricas multitemporales.
- QGIS / SAGA GIS: Cálculo de proxies topográficos avanzados.

📊 Metodología y Feature Engineering

El modelo se alimenta de un stack de 11 variables predictoras alineadas bajo el sistema de referencia EPSG:32719 (UTM 19S):
1. Variables de Vegetación (Sentinel-2)
   - Dinámica Estacional: Diferencias de $NDVI$ y $NDWI$ entre estaciones.
   - Composiciones Multitemporales: Medianas estacionales (invierno/verano) para capturar la fenología estable del bosque nativo.
2. Variables Topográficas y Geomorfológicas
   Se incorporan proxies hídricos para capturar la aptitud del terreno:
   - TWI (Topographic Wetness Index): Índice de humedad topográfica.
   - HAND (Height Above Nearest Drainage): Altura sobre el drenaje más cercano, clave para detectar proximidad a acuíferos.
   - Distancia a Ríos y Pendiente (Slope): Análisis de proximidad y morfología del terreno.

🧠 Modelado e Interpretabilidad

A diferencia de los modelos de "caja negra", este flujo de trabajo prioriza la explicabilidad:
- Algoritmo: Random Forest (Ranger) con 500 árboles y entrenamiento supervisado de presencia/ausencia.
- Valores SHAP (Shapley Additive Explanations): Implementación de fastshap para cuantificar el impacto real de cada variable (como el TWI o el HAND) en la predicción final.
- Análisis Waterfall: Visualización detallada de la contribución de cada predictor en puntos específicos del territorio.

🚀 Estructura del Repositorio

- `analysis_freatofitas.R`: Script principal de alineación, entrenamiento y generación de mapas de probabilidad.
- `data/`: Estructura para rasters de entrada (Sentinel, TWI, HAND, DEM).
- `outputs/`:
  - `Mapa_Probabilidad_Freatofitas.tif`: Cartografía final exportada.
  - `Grafico_Importancia_SHAP.png`: Visualización de la jerarquía de variables.
  - `Waterfall_Plot_General.png`: Análisis de contribución local.

👤 Contacto: Facundo Nicolas Ruarte Perez. Estudiante avanzado de Licenciatura en Biología (UNSJ).Diplomado en Ciencia de Datos (Mundos E - UNC).
📍 San Juan, Argentina.
📫 facuruarte1999@gmail.com.
