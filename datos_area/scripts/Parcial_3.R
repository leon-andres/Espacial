######################################################################
###############################
#######1. Reading dataset #######
###############################

library(sf)
library(spdep)
library(dbscan)
library(ggplot2)
library(dplyr)
library(sf)
library(spdep)
library(tibble)
library(dplyr)
library(knitr)
library(tmap)

# Lectura de los datos ####
columbus <- read_sf("columbus.shp")
columbus <- st_transform(columbus, crs = 32617)

# Matriz de vecindades ####

# Calcular los centroides
centroides <- st_centroid(columbus)

# Extraer coordenadas de los centroides
xy0 <- st_coordinates(centroides)

# Graficar el mapa
plot(st_geometry(columbus), border = "green", axes = T, col = "grey")
title(main = "Columbus Neighborhoods")
plot(st_geometry(centroides), add = TRUE, pch = 19, cex = 0.7, col = "blue")
grid() # Agregar grilla




############################################
# 2. Definir matrices de vecinos
############################################
# Matrices de pesos espaciales

# Vecinos físicos (adyacencia de polígonos)
adyacencia_rook_binaria     <- nb2listw(poly2nb(columbus, queen = FALSE), style = "B", zero.policy = TRUE)
adyacencia_rook_ponderada   <- nb2listw(poly2nb(columbus, queen = FALSE), style = "W", zero.policy = TRUE)

adyacencia_queen_binaria    <- nb2listw(poly2nb(columbus, queen = TRUE),  style = "B", zero.policy = TRUE)
adyacencia_queen_ponderada  <- nb2listw(poly2nb(columbus, queen = TRUE),  style = "W", zero.policy = TRUE)

# Vecinos por grafos
triangulacion_binaria       <- nb2listw(tri2nb(xy0), style = "B", zero.policy = TRUE)
triangulacion_ponderada     <- nb2listw(tri2nb(xy0), style = "W", zero.policy = TRUE)

soi_binaria                 <- nb2listw(graph2nb(soi.graph(tri2nb(xy0), xy0)), style = "B", zero.policy = TRUE)
soi_ponderada               <- nb2listw(graph2nb(soi.graph(tri2nb(xy0), xy0)), style = "W", zero.policy = TRUE)

relativa_binaria           <- nb2listw(graph2nb(relativeneigh(xy0), sym = TRUE), style = "B", zero.policy = TRUE)
relativa_ponderada         <- nb2listw(graph2nb(relativeneigh(xy0), sym = TRUE), style = "W", zero.policy = TRUE)

gabriel_binaria            <- nb2listw(graph2nb(gabrielneigh(xy0), sym = TRUE), style = "B", zero.policy = TRUE)
gabriel_ponderada          <- nb2listw(graph2nb(gabrielneigh(xy0), sym = TRUE), style = "W", zero.policy = TRUE)

# Vecinos por distancia (k vecinos más cercanos)
knn1_binaria               <- nb2listw(knn2nb(knearneigh(xy0, k = 1)), style = "B", zero.policy = TRUE)
knn1_ponderada             <- nb2listw(knn2nb(knearneigh(xy0, k = 1)), style = "W", zero.policy = TRUE)

knn2_binaria               <- nb2listw(knn2nb(knearneigh(xy0, k = 2)), style = "B", zero.policy = TRUE)
knn2_ponderada             <- nb2listw(knn2nb(knearneigh(xy0, k = 2)), style = "W", zero.policy = TRUE)

knn3_binaria               <- nb2listw(knn2nb(knearneigh(xy0, k = 3)), style = "B", zero.policy = TRUE)
knn3_ponderada             <- nb2listw(knn2nb(knearneigh(xy0, k = 3)), style = "W", zero.policy = TRUE)

knn4_binaria               <- nb2listw(knn2nb(knearneigh(xy0, k = 4)), style = "B", zero.policy = TRUE)
knn4_ponderada             <- nb2listw(knn2nb(knearneigh(xy0, k = 4)), style = "W", zero.policy = TRUE)

# Lista con todos los objetos
mat <- list(
  adyacencia_rook_binaria, adyacencia_rook_ponderada,
  adyacencia_queen_binaria, adyacencia_queen_ponderada,
  triangulacion_binaria, triangulacion_ponderada,
  soi_binaria, soi_ponderada,
  gabriel_binaria, gabriel_ponderada,
  relativa_binaria, relativa_ponderada,
  knn1_binaria, knn1_ponderada,
  knn2_binaria, knn2_ponderada,
  knn3_binaria, knn3_ponderada,
  knn4_binaria, knn4_ponderada
)


# Nombres para las matrices de vecinos
nombres_matrices <- c(
  "Adyacencia Rook Binaria", "Adyacencia Rook Ponderada",
  "Adyacencia Queen Binaria", "Adyacencia Queen Ponderada",
  "Triangulación Binaria", "Triangulación Ponderada",
  "SOI Binaria", "SOI Ponderada",
  "Gabriel Binaria", "Gabriel Ponderada",
  "Relativa Binaria", "Relativa Ponderada",
  "KNN 1 Binaria", "KNN 1 Ponderada",
  "KNN 2 Binaria", "KNN 2 Ponderada",
  "KNN 3 Binaria", "KNN 3 Ponderada",
  "KNN 4 Binaria", "KNN 4 Ponderada"
)

# Calcular índice de Moran y p-valores
resultados <- lapply(mat, function(w) moran.test(columbus$CRIME, w, alternative = "two.sided"))
valores_moran <- sapply(resultados, function(x) x$estimate[["Moran I statistic"]])
pvalores <- sapply(resultados, function(x) x$p.value)

# Crear la tabla
tabla_resultados <- tibble(
  Matriz = nombres_matrices,
  `Índice de Moran` = round(valores_moran, 4),
  `Valor p` = signif(pvalores, 4)
) |>
  arrange(desc(`Índice de Moran`))

# Mostrar tabla bonita en consola
kable(tabla_resultados, format = "simple", align = "lcc", caption = "Índices de Moran ordenados (de mayor a menor)")

columbus.crime <- knn2_ponderada

# Obtener lista de vecinos desde la matriz de pesos
vecinos <- knn2nb(knearneigh(xy0, k = 2))

# Crear líneas entre vecinos
segmentos <- list()
for (i in 1:length(vecinos)) {
  if (length(vecinos[[i]]) > 0) {
    for (j in vecinos[[i]]) {
      # Coordenadas de los puntos conectados
      coords_i <- xy0[i, ]
      coords_j <- xy0[j, ]
      # Crear líneas entre ellos
      segmentos <- append(segmentos, list(st_linestring(rbind(coords_i, coords_j))))
    }
  }
}

# Convertir a objeto sf para graficar
lineas_vecinas <- st_sfc(segmentos, crs = st_crs(columbus))

# Graficar mapa con vecindades
plot(st_geometry(columbus), col = "lightgrey", border = "darkgreen", main = "Conectividad KNN (k=2)")
plot(lineas_vecinas, col = "red", lwd = 1.5, add = TRUE)
plot(st_geometry(centroides), add = TRUE, pch = 19, col = "blue", cex = 0.6)


# # Convertir la matriz de vecinos (listw) a una matriz cuadrada clásica
# W <- listw2mat(best_matrix)
# 
# # Ahora ya tienes W_best como una matriz n × n clásica (por ejemplo 49 × 49)
# # Puedes revisarla con
# print(W)
# dim(W)
# 
# Z <- (as.vector(columbus[,9]))$CRIME
# WZ <- W %*% Z
# 
# # Asegurar que 'columbus_sf' sea un objeto sf válido
# columbus_sf <- st_as_sf(columbus)
# 
# # Crear una lista para almacenar los data.frames individuales
# mapas_df <- list()
# 
# # Mapa original: agregar la variable CRIME con etiqueta W0 para mantener el mismo formato
# columbus_sf$CRIME_W0 <- columbus$CRIME
# mapas_df[[1]] <- columbus_sf %>%
#   mutate(W = "Original", valor = CRIME_W0) %>%
#   dplyr::select(W, valor, geometry)
# 
# # Calcular y agregar medias móviles para cada matriz de pesos según la fórmula:
# # μ̂_i = (sum_j w_ij * y_j) / (sum_j w_ij)
# for (i in seq_along(mat)) {
#   W <- mat[[i]]
#   
#   # Numerador: suma ponderada de CRIME
#   num <- lag.listw(W, columbus$CRIME, zero.policy = TRUE)
#   
#   # Denominador: suma de pesos por fila
#   denom <- sapply(1:length(W$weights), function(i) sum(W$weights[[i]]))
#   
#   # Media móvil según fórmula dada
#   mu_hat <- num / denom
#   
#   # Crear data.frame con nombre de la matriz
#   temp_df <- columbus_sf %>%
#     mutate(W = paste0("W", i), valor = mu_hat) %>%
#     dplyr::select(W, valor, geometry)
#   
#   mapas_df[[i + 1]] <- temp_df
# }
# 
# # Unir todos los data.frames en uno solo
# mapas_completos <- do.call(rbind, mapas_df)
# 
# # Graficar con ggplot2
# ggplot(mapas_completos) +
#   geom_sf(aes(fill = valor), color = "white") +
#   facet_wrap(~W) +
#   scale_fill_viridis_c(option = "C", name = "Media Móvil") +
#   labs(title = "Mapas de medias móviles espaciales por matriz W") +
#   theme_minimal()

# Variables de interés para el test de Moran
vars_interes <- c("HOVAL", "INC", "OPEN", "PLUMB", "DISCBD")

# Lista para guardar las mejores matrices
mejores_matrices <- list()

# Test de Moran para cada variable y guardar su mejor matriz
for (var in vars_interes) {
  aux <- numeric(length(mat))  # Reiniciar aux en cada variable
  for (i in 1:length(mat)) {
    aux[i] <- moran.test(columbus[[var]], mat[[i]], alternative = "two.sided")$p.value
  }
  best_i <- which.min(aux)
  cat("Mejor matriz para", var, ":", best_i, "\n")
  print(moran.test(columbus[[var]], mat[[best_i]], alternative = "two.sided"))
  
  # Guardar la mejor matriz en la lista, usando el nombre de la variable
  mejores_matrices[[var]] <- mat[[best_i]]
}

# También puedes crear variables independientes si prefieres eso:
columbus.hoval  <- mejores_matrices[["HOVAL"]]
columbus.inc    <- mejores_matrices[["INC"]]
columbus.open   <- mejores_matrices[["OPEN"]]
columbus.plumb  <- mejores_matrices[["PLUMB"]]
columbus.discbd <- mejores_matrices[["DISCBD"]]

# Mapa de valores observados ####


# Activar modo de visualización estática
tmap_mode("plot")

# Mapa temático de crimen
tm_shape(columbus) +
  tm_fill("CRIME", 
          palette = "Reds", 
          style = "quantile", 
          title = "Delitos por 1000 hogares") +
  tm_borders() +
  tm_layout(main.title = "Crimen residencial y automotor en Columbus, Ohio (1980)", 
            legend.outside = TRUE)



# Dispersograma de Moran
moran.plot(columbus$CRIME, 
           columbus.crime, 
           labels=as.character(columbus$NEIG), 
           xlab="Delitos por 1000 hogares", 
           ylab="Delitos por 1000 hogares rezagado", 
           las=1, 
           pch=16, 
           cex=0.5)

legend("bottomright", 
       legend=c("I de Moran: 0.6073", "Valor P:      0.0000015"), 
       cex=1,
       bg='lightgreen')

title("Dispersograma de Moran para el Crimen residencial y automotor en Columbus, Ohio (1980)", cex.main=1)



# Calcular el estadístico Local G para la variable CRIME en el dataset columbus
localG_crime <- localG(columbus$CRIME, columbus.crime)

# Agregar el resultado como nueva columna al objeto sf
columbus$LocalG_CRIME <- localG_crime


# 2. Crear matriz vacía para guardar resultados simulados
sim.G <- matrix(0, nrow = 1000, ncol = length(localG_crime))

# 3. Repetir permutaciones
set.seed(123)  # para reproducibilidad
for (i in 1:1000) {
  sim.G[i, ] <- localG(sample(columbus$CRIME), columbus.crime)
}

# 4. Calcular p-valores empíricos
mc.pvalor.G <- (colSums(sweep(sim.G, 2, localG_crime, ">=")) + 1) / (nrow(sim.G) + 1)

# 5. Agregar resultado al dataset
columbus$MC_PVALOR_LOCALG <- mc.pvalor.G

# Paso 6: Visualizar Local G
tm_shape(columbus) +
  tm_fill("LocalG_CRIME", 
          palette = "-RdBu",
          style = "pretty", 
          title = "Estadístico Local G") +
  tm_borders() +
  tm_layout(main.title = "Hotspots y Coldspots del crimen en Columbus (Local G)", 
            legend.outside = TRUE)

tm_shape(columbus) +
  tm_fill("MC_PVALOR_LOCALG", 
          palette = "Reds",        # Valores bajos = más significativos
          style = "cont", 
          title = "p-valor Monte Carlo") +
  tm_borders() +
  tm_layout(main.title = "Significancia del Local G (Monte Carlo)", 
            legend.outside = TRUE)


######################################################################


## Modelo de regresión básico ####
# Ajustar modelo de regresión lineal clásico
modelo_ind <- lm(CRIME ~ HOVAL + INC + DISCBD + OPEN + PLUMB, data = columbus)




library(spdep)
library(spatialreg)
####Modelos SDEM, SDM, Manski, SARAR########
#reg.eq1=CAP_BAC ~ PIB + NBI + CAP_BOG + CAP_BC + CAP_OCC + CAP_CS + Población + IPM
reg.eq1=CRIME ~ HOVAL + INC + DISCBD + OPEN + PLUMB +CP
reg1=lm(reg.eq1,data=columbus)                                     #OLS            y=XB+e,    
reg2=lmSLX(reg.eq1,data=columbus, columbus.crime)                       #SLX            y=XB+WxT+e
reg3=lagsarlm(reg.eq1,data= columbus, columbus.crime)                   #Lag Y          y=XB+WxT+u,   u=LWu+e
reg4=errorsarlm(reg.eq1,data=columbus, columbus.crime)                  #Spatial Error  y=pWy+XB+e   
reg5=errorsarlm(reg.eq1, data=columbus, columbus.crime, etype="emixed") #SDEM Spatial Durbin Error Model y=XB+WxT+u,   u=LWu+e
reg6=lagsarlm(reg.eq1, data=columbus,columbus.crime, type="mixed")      #SDM Spatial Durbin Model (add lag X to SAR) y=pWy+XB+WXT+e 
reg7=sacsarlm(reg.eq1,data=columbus, columbus.crime, type="sacmixed")   #Manski Model: y=pWy+XB+WXT+u,   u=LWu+e (no recomendado)
reg8=sacsarlm(reg.eq1,data=columbus, columbus.crime, type="sac")         #SARAR o Kelejian-Prucha, Cliff-Ord, o SAC If all T=0,y=pWy+XB+u, u=LWu+e

columbus.crime <- knn4_ponderada
s=summary
s(reg1)
s(reg2)
s(reg3)
s(reg4)
s(reg5)
s(reg6)
s(reg7)
s(reg8)


lm.morantest(reg1, columbus.crime)

AIC(reg1, reg2, reg3, reg4, reg5, reg6, reg7, reg8)

modelo1 <- errorsarlm(reg.eq1,data=columbus, knn2_ponderada) 
modelo2 <- lagsarlm(reg.eq1,data=columbus, knn4_ponderada) 
modelo3 <- errorsarlm(reg.eq1,data=columbus, relativa_binaria) 
modelo4 <- errorsarlm(reg.eq1,data=columbus, knn1_ponderada) 
reg.eq1=CRIME ~ HOVAL + INC + DISCBD + OPEN + PLUMB +CP


s(modelo1)

reg.eq2= CRIME ~ HOVAL + INC + DISCBD + PLUMB +CP
modelo1=errorsarlm(reg.eq2,data=columbus, knn2_ponderada)                  #Spatial Error  y=pWy+XB+e   
s(modelo1)

reg.eq2= CRIME ~ HOVAL + INC + PLUMB +CP
modelo1=errorsarlm(reg.eq2,data=columbus, knn2_ponderada)                  #Spatial Error  y=pWy+XB+e   
s(modelo1)

reg.eq2= CRIME ~ HOVAL + INC +CP
modelo1=errorsarlm(reg.eq2,data=columbus, knn2_ponderada)                  #Spatial Error  y=pWy+XB+e   
s(modelo1)


# Modelo con filtrado espacial usando eigenvectores
filtro <- ME(CRIME ~ INC + HOVAL + CP, data = columbus, listw = knn2_ponderada)
summary(filtro)
str(filtro)

columbus$filtro1 <- filtro$vectors

reg.eq2= CRIME ~ HOVAL + INC +CP + filtro1
modelo1=errorsarlm(reg.eq2,data=columbus, knn2_ponderada)                  #Spatial Error  y=pWy+XB+e   
s(modelo1)

# Extrae los residuos del modelo SAC
resid_sac <- residuals(modelo1)

library(spatialreg)
bptest.Sarlm(modelo1)


# Prueba de Moran's I
moran.test(resid_sac, listw = columbus.crime)

library(lmtest)

# Extraer los residuos del modelo 3
residuos <- resid(modelo3)

# Crear el modelo para BP test
bp_model <- lm(I(residuos^2) ~ INC + HOVAL + PLUMB, data = columbus)

# Aplicar la prueba de Breusch-Pagan
bptest(bp_model)



# Gráficos clásicos
plot(residuos)

# QQ plot
qqnorm(residuos, main = "QQ plot de los residuos (modelo SAC)")
qqline(residuos, col = "red", lwd = 2)


plot(fitted(modelo2), columbus$CRIME, main = "Predichos vs Observados", xlab = "Predichos", ylab = "Observados")
abline(0,1,col="red")

shapiro.test(resid_sac)

plot(fitted(reg6), resid_sac,
     xlab = "Valores ajustados",
     ylab = "Residuos",
     main = "Residuos vs Ajustados (SAC)",
     pch = 20, col = "steelblue")
abline(h = 0, col = "red", lwd = 2)

library(lmtest)
bptest(reg6)
bptest.sarlm(modelo1)


library(spatialreg)
bptest.Sarlm(residuos)


# Modelo con filtrado espacial usando eigenvectores
filtro <- ME(CRIME ~ INC + HOVAL, data = columbus, listw = knn2_ponderada)
summary(filtro)
str(filtro)

# Extraer residuos
residuos <- residuals(modelo_ind)

# Agregar los residuos al objeto spatial
columbus$residuos <- residuos

# Mapa rápido de residuos
library(tmap)
tmap_mode("view")
tm_shape(columbus) +
  tm_fill("residuos", palette = "RdBu", style = "quantile") +
  tm_borders()


# Cargar paquetes si no los tienes
library(spdep)

# Crear matriz de vecinos por contigüidad
vecinos <- poly2nb(columbus)

# Crear lista de pesos espaciales estandarizada
lw <- nb2listw(vecinos, style = "W")

# Calcular el Moran's I de los residuos
moran.test(columbus$residuos, lw)


























# Cargar paquetes necesarios
library(sf)
library(spdep)
library(classInt)
library(RColorBrewer)
library(tmap)


# Obtener los valores ajustados del modelo 1
ajustados <- modelo1$fitted.values

# Agregar los valores ajustados al objeto sf
columbus$ajustados <- ajustados

# Calcular cuantiles para 4 grupos
breaks_q <- classIntervals(ajustados, n = 5, style = "quantile")$brks

# Crear paleta de colores
paleta <- brewer.pal(n = 5, name = "Reds")

# Mapa choropleth usando tmap
tm_shape(columbus) +
  tm_fill("ajustados",
          palette = paleta,
          breaks = breaks_q,
          title = "Valores ajustados") +
  tm_borders() +
  tm_layout(title = "Valores ajustados del Modelo 1 en Columbus",
            legend.outside = TRUE)

