library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)
library(sp)
library(spacetime)
library(gstat)
library(readxl)

lista_librerías<-c("readr","sqldf","tidyr","tidyverse","ggplot2","patchwork","sp","fields","geoR","akima","shiny","gstat", "sf", "dplyr", "knitr", "raster", "scales", "reshape2", "broom") 


no_installs <- lista_librerías[!lista_librerías %in% installed.packages()]

if(length(no_installs) > 0) {
  cat("Los siguientes paquetes no están instalados :\n")
  cat(no_installs, sep = "\n")
  install.packages(no_installs)
} else {
  cat("Todos los paquetes están instalados. \n")
}

sapply(lista_librerías, require, character = TRUE)



colom_shp <- suppressMessages(st_read("mapa/loca.shp"))  # Ocultar mensajes de st_read
colom_shp <- st_transform(colom_shp[-9,], crs = "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")



CRS_UTM_CO <- CRS("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")


PM2_5 <- read_xlsx("PM2_5.xlsx")

# Para PM2.5
PM2_5_sf <- st_as_sf(PM2_5, coords = c("este", "norte"), crs = CRS_UTM_CO)

ggplot() +
  geom_sf(data = colom_shp, fill = "white", color = "black") +
  geom_sf(data = PM2_5_sf, aes(color = PM2.5), alpha = 0.8) +
  scale_color_gradient(name = "PM2.5 (µg/m³)", low = "blue", high = "red") +
  labs(title = "Mapa de estaciones PM2.5 - Bogotá") +
  theme_minimal()

# Gráfico 1: PM2.5 vs Este
p1 <- ggplot(PM2_5, aes(x = este, y = PM2.5)) +
  geom_point(color = "blue") +
  labs(title = "PM2.5 vs Este", x = "Este (X)", y = "PM2.5") +
  theme_minimal()

# Gráfico 2: PM2.5 vs Norte
p2 <- ggplot(PM2_5, aes(x = norte, y = PM2.5)) +
  geom_point(color = "blue") +
  labs(title = "PM2.5 vs Norte", x = "Norte (Y)", y = "PM2.5") +
  theme_minimal()

# Gráfico 3: PM2.5 vs Fecha
p3 <- ggplot(PM2_5, aes(x = fecha, y = PM2.5)) +
  geom_point(color = "darkred", alpha = 0.5) +
  labs(title = "PM2.5 a lo largo del tiempo", x = "Fecha", y = "PM2.5") +
  theme_minimal()

# Unir los tres gráficos en un solo panel
(p1 | p2) / p3


 fit1 <- lm(PM2.5 ~ sin(este) + sin(norte)+este+norte, data = PM2_5)
PM2_5$residuales <- residuals(fit1)

res1 <- PM2_5$residuales
par(mfrow = c(1, 3))

# Residuales estandarizados vs Coordenada Este
plot(PM2_5$este, res1,
     main = "Residuales vs Coordenada Este",
     xlab = "Este",
     ylab = "Residuales estandarizados",
     pch = 16, col = "darkblue")

# Residuales estandarizados vs Coordenada Norte
plot(PM2_5$norte, res1,
     main = "Residuales vs Coordenada Norte",
     xlab = "Norte",
     ylab = "Residuales estandarizados",
     pch = 16, col = "darkgreen")

# Boxplot de los residuales estandarizados
boxplot(res1,
        main = "Distribución de Residuales",
        ylab = "Residuales estandarizados",
        col = "red",
        ylim = c(-25, 25))

par(mfrow = c(1, 1))

PM2_5$residuals <- res1



# Preparar datos
PM2_5_clean <- PM2_5
PM2_5_clean$Fecha <- as.POSIXct(PM2_5$fecha)

# Convertir a SpatialPoints
coordinates(PM2_5_clean) <- ~este + norte

# Crear vectores de espacio y tiempo
spatial_points <- unique(PM2_5_clean@coords)  # puntos únicos
times <- sort(unique(PM2_5_clean$Fecha))      # tiempos únicos

# Crear SpatialPoints para la malla espacial
sp <- SpatialPoints(spatial_points)

# Crear data.frame con variable residuals ordenada para el STFDF:
# Primero, crear un data.frame vacío con filas = puntos espaciales, columnas = tiempos
# y rellenar con NA
residuals_matrix <- matrix(NA, nrow = length(sp), ncol = length(times),
                           dimnames = list(NULL, as.character(times)))

# Llenar residuals_matrix con los valores correspondientes
for (i in seq_len(nrow(PM2_5_clean))) {
  # Encontrar fila y columna
  point_index <- which(apply(sp@coords, 1, function(r) all(r == PM2_5_clean@coords[i, ])))
  time_index <- which(times == PM2_5_clean$Fecha[i])
  residuals_matrix[point_index, time_index] <- PM2_5_clean$residuals[i]
}

# Crear data.frame para STFDF (aplanar la matriz en un vector)
data_df <- data.frame(residuals = as.vector(residuals_matrix))

# Crear objeto STFDF
st_data <- STFDF(sp = sp, time = times, data = data_df)

# Ahora sí, calcular el variograma espaciotemporal
v_st <- variogramST(
  residuals ~ 1,
  data = st_data,
  tunit = "days",
  cutoff = 19000, #19000
  width = 2400, #2400
  tlags = 0:5,
  assumeRegular = TRUE
)
v_st <- subset(v_st, np > 0)


plot(v_st)

plot(v_st, wireframe = TRUE, colorkey = TRUE,
     main = "Semivariograma Espaciotemporal")



# Definir modelo inicial del variograma espacio-temporal (ajustable)
model_init <- vgmST(
  "separable",
  space = vgm(psill = 40, model = "Hol", range = 6000, nugget = 0.1),
  time  = vgm(psill = 1, model = "Pen", range = 8,       nugget = 0.1),
  sill  = 1
)

# Ajustar el modelo de variograma a los datos observados
model_fit <- fit.StVariogram(v_st, model_init, method = "L-BFGS-B")

plot(v_st,model_fit)

plot(v_st, model_fit, wireframe = TRUE, all = TRUE,
     colorkey = TRUE, main = "Ajuste del Modelo ST al Semivariograma Empírico")

# 1. Asegurarse de que el shapefile esté en CRS correcto
colom_shp <- st_transform(colom_shp, crs = CRS_UTM_CO)

# 2. Obtener bounding box del shapefile como objeto sf
bbox_colom <- st_bbox(colom_shp)

# 3. Crear malla de predicción en todo el país (en grilla regular)
x_range <- seq(bbox_colom["xmin"], bbox_colom["xmax"], by = 500)
y_range <- seq(bbox_colom["ymin"], bbox_colom["ymax"], by = 500)
grilla_espacio <- expand.grid(este = x_range, norte = y_range)

# 4. Convertir grilla a sf para filtrar solo puntos dentro de Colombia
grilla_sf <- st_as_sf(grilla_espacio, coords = c("este", "norte"), crs = CRS_UTM_CO)
grilla_dentro <- st_join(grilla_sf, colom_shp, join = st_within, left = FALSE)

# 5. Convertir de nuevo a Spatial (lo que requiere krigeST)
grilla_sp <- as(grilla_dentro, "Spatial")

gridded(grilla_sp) <- TRUE

# 6. Crear secuencia de tiempo
tiempos_pred <- seq(as.POSIXct("2020-10-16"), as.POSIXct("2021-03-31"), by = "10 days")

# 7. Malla espacio-temporal
grilla_st <- STF(sp = grilla_sp, time = tiempos_pred)

# EXTRA: Asegurar que el CRS de st_data@sp esté definido (si no lo tiene, asignarlo)
if (is.na(proj4string(st_data@sp))) {
  proj4string(st_data@sp) <- CRS_UTM_CO
}

# EXTRA: Asegurar que el CRS de st_data es el mismo que el de grilla_st
if (!identicalCRS(st_data@sp, grilla_sp)) {
  st_data@sp <- spTransform(st_data@sp, CRS_UTM_CO)
}
# 
gc()
# 8. Ejecutar kriging
kriging_result <- krigeST(
  residuals ~ 1,
  data = st_data,
  newdata = grilla_st,
  modelList = model_fit
)
gc()


graficar_kriging_con_tendencia <- function(fechas, kriging_result, fit1, colom_shp) {
  library(dplyr)
  library(ggplot2)
  
  # Paleta de 20 colores ordenados (oscuro a claro)
  colores_ordenados <- c(
    "#000000", "#2E0854", "#4B0082", "#0000CD", "#1E90FF",
    "#00CED1", "#00FF7F", "#7FFF00", "#ADFF2F", "#FFFF00",
    "#FFD700", "#FFA500", "#FF8C00", "#FF6347", "#FF4500",
    "#FF0000", "#DC143C", "#FF69B4", "#FFC0CB", "#FFFFFF"
  )
  
  fechas_xts <- index(kriging_result@time)
  
  # --- 1. Recolectar todos los valores PM25_estimado para calcular cuantiles ---
  todos_est <- data.frame()
  
  for (f in fechas) {
    fecha_objetivo <- as.POSIXct(f)
    if (!fecha_objetivo %in% fechas_xts) next
    kriging_en_fecha <- kriging_result[, fecha_objetivo]
    coords <- kriging_en_fecha@coords
    valores <- kriging_en_fecha@data$var1.pred
    df_tmp <- data.frame(x = coords[, 1], y = coords[, 2])
    df_tmp$este <- df_tmp$x
    df_tmp$norte <- df_tmp$y
    df_tmp$layer <- valores
    df_tmp$PM25_estimado <- df_tmp$layer + predict(fit1, newdata = df_tmp[, c("este", "norte")])
    todos_est <- bind_rows(todos_est, df_tmp)
  }
  
  # --- 2. Calcular los 20 cuantiles y etiquetas ---
  breaks <- quantile(todos_est$PM25_estimado, probs = seq(0, 1, length.out = 21), na.rm = TRUE)
  
  # Crear etiquetas del tipo "[a, b)"
  etiquetas <- sapply(1:20, function(i) {
    paste0("[", round(breaks[i], 1), ", ", round(breaks[i + 1], 1), ")")
  })
  
  # --- 3. Función para asignar tramo y color ---
  asignar_tramo <- function(valor) {
    idx <- findInterval(valor, breaks, rightmost.closed = TRUE, all.inside = TRUE)
    etiquetas[idx]
  }
  
  # --- 4. Graficar cada fecha ---
  for (f in fechas) {
    fecha_objetivo <- as.POSIXct(f)
    if (!fecha_objetivo %in% fechas_xts) {
      cat("Fecha no disponible:", f, "\n")
      next
    }
    
    cat("Graficando:", f, "\n")
    kriging_en_fecha <- kriging_result[, fecha_objetivo]
    coords <- kriging_en_fecha@coords
    valores <- kriging_en_fecha@data$var1.pred
    kriging_df <- data.frame(x = coords[, 1], y = coords[, 2], layer = valores)
    kriging_df$este <- kriging_df$x
    kriging_df$norte <- kriging_df$y
    kriging_df$PM25_estimado <- kriging_df$layer + predict(fit1, newdata = kriging_df[, c("este", "norte")])
    kriging_df <- na.omit(kriging_df)
    
    # Asignar tramo (cuantil) como factor
    kriging_df$tramo <- factor(asignar_tramo(kriging_df$PM25_estimado), levels = etiquetas)
    
    # Gráfico con leyenda manual

    p <- ggplot() +
      geom_point(data = kriging_df, aes(x = x, y = y, color = tramo), size = 1.8, alpha = 0.85) +
      scale_color_manual(
        values = setNames(colores_ordenados, etiquetas),
        name = "PM2.5 estimado (µg/m³)"
      ) +
      geom_sf(data = colom_shp, fill = NA, color = "black", linewidth = 0.3) +
      coord_sf() +
      labs(
        title = paste("PM2.5 estimado (por cuantiles) -", f),
        x = "Este", y = "Norte"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(size = 10)
      )
    print(p)
    Sys.sleep(2)
  }
}

fechas_deseadas <- c(
  "2020-10-16", "2020-10-26", "2020-11-05", "2020-11-15",
  "2020-11-25", "2020-12-05", "2020-12-15", "2020-12-25",
  "2021-01-04", "2021-01-14", "2021-01-24", "2021-02-03",
  "2021-02-13", "2021-02-23", "2021-03-05", "2021-03-15"
)

graficar_kriging_con_tendencia(fechas_deseadas, kriging_result, fit1, colom_shp)
