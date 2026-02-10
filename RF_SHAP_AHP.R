### Analisis de freatofitas - Entrenamiento con TWI y Variables Locales

#########################
rm(list=ls())
gc()

# ==============================================================================
# 1. CARGA DE LIBRERÍAS
# ==============================================================================
library(terra)      # Manejo de Rasters
library(sf)         # Manejo de Vectores
library(tidyverse)  # Manipulación de datos
library(ranger)     # Random Forest
library(fastshap)   # SHAP values

# ==============================================================================
# 2. CARGA DE DATOS (Sentinel y Topografía)
# ==============================================================================

# --- A. Variables de Vegetación (Sentinel) ---
dif_ndvi <- rast("S2_Diferencia_NDVI.tif") 
dif_ndwi <- rast("S2_Diferencia_NDWI.tif") 
inv_ndvi <- rast("S2_Invierno_NDVI_SanJuan.tif") 
inv_ndwi <- rast("S2_Invierno_NDWI_SanJuan.tif") 
ver_ndvi <- rast("S2_Verano_NDVI_SanJuan.tif") 
ver_ndwi <- rast("S2_Verano_NDWI_SanJuan.tif") 

# --- B. Variables Topográficas (Calculadas en QGIS/SAGA) ---
twi       <- rast("TWI (QGIS).tif")
dist_rios <- rast("distancia_horizontal_rios.sdat")
slope     <- rast("Slope.sdat")
hand      <- rast("proxy HAND.sdat")
dem       <- rast("dem_formatoraster_reproyectado.tif")

# ==============================================================================
# 2.C ALINEACIÓN Y REPROYECCIÓN (Priorizando el DEM - EPSG:32719)
# ==============================================================================

# 1. Definimos la referencia maestra (el DEM corregido)
ref_maestra <- dem

# 2. Reproyectamos todas las capas de Sentinel al CRS y rejilla del DEM
# Usamos bilinear para variables continuas (NDVI, NDWI)
dif_ndvi_res <- project(dif_ndvi, ref_maestra, method = "bilinear")
dif_ndwi_res <- project(dif_ndwi, ref_maestra, method = "bilinear")
inv_ndvi_res <- project(inv_ndvi, ref_maestra, method = "bilinear")
inv_ndwi_res <- project(inv_ndwi, ref_maestra, method = "bilinear")
ver_ndvi_res <- project(ver_ndvi, ref_maestra, method = "bilinear")
ver_ndwi_res <- project(ver_ndwi, ref_maestra, method = "bilinear")

# 3. Las capas topográficas ya están en el CRS correcto, pero usamos resample 
# por si tienen extensiones ligeramente diferentes tras el Fill Sinks
twi_res       <- resample(twi, ref_maestra, method = "bilinear")
dist_rios_res <- resample(dist_rios, ref_maestra, method = "bilinear")
slope_res     <- resample(slope, ref_maestra, method = "bilinear")
hand_res      <- resample(hand, ref_maestra, method = "bilinear")
# El DEM no necesita resample consigo mismo

# --- D. CREAR EL STACK FINAL ---
stack_predictores <- c(dif_ndvi_res, dif_ndwi_res, inv_ndvi_res, inv_ndwi_res, 
                       ver_ndvi_res, ver_ndwi_res, twi_res, dist_rios_res, 
                       slope_res, hand_res, ref_maestra)

names(stack_predictores) <- c("Dif_NDVI", "Dif_NDWI", "Inv_NDVI", "Inv_NDWI", 
                              "Ver_NDVI", "Ver_NDWI", "TWI", "Dist_Rios", 
                              "Slope", "HAND", "DEM")
print(paste("Stack listo con", nlyr(stack_predictores), "variables."))

# ==============================================================================
# 3. PREPARACIÓN DE DATOS DE ENTRENAMIENTO (Puntos de Freatofitas)
# ==============================================================================
# Es CRÍTICO que los puntos se muevan al EPSG:32719 antes de extraer
# 1. Cargar como sf para asegurar que leemos bien la tabla de atributos
puntos_sf <- st_read("pronis_datos_pa.geojson") %>% 
  st_transform(32719)

# 2. Convertir a objeto terra para el extract
puntos_terra <- vect(puntos_sf)


# 3. Extraer los valores de los rasters (ya obtenidos como puntos_terra)
datos_raster <- terra::extract(stack_predictores, puntos_terra, ID = FALSE)

# 4. Unir con la columna correcta: "PresFrea"
# Usamos puntos_sf que es donde reside la tabla de atributos original
datos_entrenamiento <- cbind(datos_raster, PresFrea = puntos_sf$PresFrea)

# 5. Limpieza de NAs (puntos fuera del mapa) y conversión a Factor
datos_entrenamiento <- na.omit(as.data.frame(datos_entrenamiento))
datos_entrenamiento$PresFrea <- as.factor(datos_entrenamiento$PresFrea)

print(paste("Entrenando con", nrow(datos_entrenamiento), "puntos validos."))
print(table(datos_entrenamiento$PresFrea)) # Verifica cuántas presencias (1) y ausencias (0) quedaron

# ==============================================================================
# 4. ENTRENAMIENTO RANDOM FOREST
# ==============================================================================

set.seed(123)
modelo_rf <- ranger(
  formula = PresFrea ~ .,           # Ahora la variable objetivo es PresFrea
  data = datos_entrenamiento, 
  num.trees = 500,
  importance = "permutation",
  probability = TRUE
)

# ==============================================================================
# 5. CÁLCULO DE VALORES SHAP
# ==============================================================================

# 1. Redefinimos la función wrapper de forma robusta
pfun <- function(object, newdata) {
  # Obtenemos todas las predicciones
  preds <- predict(object, data = newdata)$predictions
  
  # Buscamos cuál es la columna de la presencia (clase "1")
  # Si tus clases son 0 y 1, esto selecciona la probabilidad del 1
  return(preds[, "1"]) 
}

# 2. Si el error persiste, intenta esta versión que usa el índice de la columna
# (Asumiendo que la segunda columna es la presencia/clase 1)
pfun_index <- function(object, newdata) {
  predict(object, data = newdata)$predictions[, 2]
}
shap_values <- fastshap::explain(
  modelo_rf, 
  X = subset(datos_entrenamiento, select = -PresFrea), # Quitamos PresFrea de las variables
  pred_wrapper = pfun_index, 
  nsim = 20 
)

# Convertir los valores SHAP a una tabla de importancia
# Extraemos los nombres de las variables predictoras (sin la clase)
nombres_variables <- names(subset(datos_entrenamiento, select = -PresFrea))

# Creamos la tabla de importancia usando esos nombres
imp_df <- data.frame(
  Variable = nombres_variables,
  Importancia = colMeans(abs(as.matrix(shap_values))) # Forzamos a matriz para el cálculo
) %>% 
  arrange(desc(Importancia))

# Verificamos el resultado
print(imp_df)

# Graficar con ggplot2
ggplot(imp_df, aes(x = reorder(Variable, Importancia), y = Importancia)) +
  geom_col(fill = "steelblue", width = 0.7) +
  coord_flip() +
  labs(
    title = "Importancia de las Variables Predictoras",
    subtitle = "Basado en valores SHAP (Random Forest)",
    x = "Variable",
    y = "Importancia Media (|SHAP|)"
  ) +
  theme_minimal()

# ==============================================================================
# 4B. wATERFALL PLOT
library(shapviz)
library(ggplot2)

# 1. Definimos los datos X (solo predictores)
X_entrenamiento <- subset(datos_entrenamiento, select = -PresFrea)

# 2. Calculamos el valor base (baseline)
# Es el promedio de las predicciones de probabilidad de presencia (clase 1)
predicciones_prob <- predict(modelo_rf, data = X_entrenamiento)$predictions[, 2]
baseline_mu <- mean(predicciones_prob)

# 3. Creamos el objeto compatible
# Forzamos shap_values a ser matriz y le asignamos los nombres de variables
# 1. Aseguramos que S_mat sea una matriz y tenga nombres de columnas
S_mat <- as.matrix(shap_values)
colnames(S_mat) <- names(X_entrenamiento)

# 2. Constructor corregido: Pasamos la matriz como primer argumento (objeto)
sv_obj <- shapviz(
  S_mat,                # El objeto SHAP va primero y sin "S ="
  X = X_entrenamiento, 
  baseline = baseline_mu
)

# 3. Ahora el waterfall debería funcionar perfectamente
sv_waterfall(
  sv_obj, 
  row_id = 1,          # Fila de tu tabla de datos
  max_display = 11,    # Queremos ver todas tus variables (tienes 11)
  fill_colors = c("#2471A3", "#E74C3C") # Azul para negativo, Rojo para positivo (puedes invertirlos)
) + 
  labs(title = "Análisis Waterfall: Contribución a la Predicción",
       subtitle = "Punto de observación 1") +
  theme_minimal()
#Discrimnar por zonas
# ==============================================================================
# 5. RESULTADOS
# ==============================================================================
# 1. Generar la predicción sobre todo el stack de rasters
# En ranger con probability = TRUE, 'response' devuelve las probabilidades
prediccion_mapa <- predict(stack_predictores, modelo_rf, 
                           type = "response")

# 2. Seleccionar la banda de interés
# Ranger devuelve un raster con tantas capas como clases tengas (0 y 1)
# Si la clase 1 es freatofita, seleccionamos la capa 2
probabilidad_freatofitas <- prediccion_mapa[[2]]

# 3. Guardar el archivo para QGIS
writeRaster(probabilidad_freatofitas, "Mapa_Probabilidad_Freatofitas.tif", overwrite = TRUE)

# 4. Ver en R
plot(probabilidad_freatofitas, main = "Probabilidad de presencia de Freatofitas")

###
# Definir la ruta y nombre del archivo
ruta_salida <- "Mapa_Probabilidad_Freatofitas_UTM19S.tif"

# Guardar el raster
# Usamos datatype = 'FLT4S' para asegurar precisión decimal (0 a 1)
writeRaster(probabilidad_freatofitas, 
            filename = ruta_salida, 
            overwrite = TRUE,
            datatype = 'FLT4S', 
            gdal = c("COMPRESS=DEFLATE", "TFW=YES"))

print(paste("Mapa exportado exitosamente en:", getwd()))

#### Exportacion de graficos
# Primero generamos el gráfico y lo asignamos a un objeto
grafico_importancia <- ggplot(imp_df, aes(x = reorder(Variable, Importancia), y = Importancia)) +
  geom_col(fill = "steelblue", width = 0.7) +
  coord_flip() +
  labs(
    title = "Importancia de las Variables Predictoras",
    subtitle = "Basado en valores SHAP (Random Forest)",
    x = "Variable",
    y = "Importancia Media (|SHAP|)"
  ) +
  theme_minimal()

# Guardar el gráfico
ggsave("Grafico_Importancia_SHAP.png", 
       plot = grafico_importancia, 
       width = 8, height = 6, dpi = 300)

# 1. Encontramos la fila cuya predicción está más cerca de la media global
predicciones <- predict(modelo_rf, data = X_entrenamiento)$predictions[, 2]
id_promedio <- which.min(abs(predicciones - mean(predicciones)))

# 2. Creamos el gráfico para ese punto representativo
grafico_waterfall_gral <- sv_waterfall(
  sv_obj, 
  row_id = id_promedio, 
  max_display = 11,
  fill_colors = c("#2471A3", "#E74C3C") # Azul (negativo) y Rojo (positivo)
) + 
  labs(
    title = "Contribución de Variables (Waterfall Plot General)",
    subtitle = paste("Análisis basado en el punto representativo ID:", id_promedio),
    x = "Impacto SHAP (log-odds)"
  ) +
  theme_minimal()
print(grafico_waterfall_gral)
# 3. Exportar en alta resolución
ggsave("Waterfall_Plot_General.png", 
       plot = grafico_waterfall_gral, 
       width = 9, height = 7, dpi = 300)
####
# Tabla de pesos para análisis AHP (todavia sin debugg)
importancia_shap <- data.frame(
  Variable = names(shap_values),
  Mean_Abs_SHAP = colMeans(abs(shap_values))
) %>% 
  mutate(Peso_Relativo = Mean_Abs_SHAP / sum(Mean_Abs_SHAP)) %>% 
  arrange(desc(Peso_Relativo))

print(importancia_shap)
write.csv(importancia_shap, "Pesos_Variables_Finales.csv", row.names = FALSE)




