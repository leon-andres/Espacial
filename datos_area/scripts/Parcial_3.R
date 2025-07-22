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


# Leer el shapefile
columbus <- read_sf("columbus.shp")

# Transformar al sistema CRS WGS84
columbus <- st_transform(columbus, crs = 4326)

# Calcular los centroides como objeto sf (¡sin sp!)
centroides <- st_centroid(columbus)

# Extraer coordenadas de los centroides
xy0 <- st_coordinates(centroides)

# Graficar el mapa
plot(st_geometry(columbus), border = "grey", axes = TRUE, col = "white")
title(main = "Columbus Neighborhoods")

# Agregar los puntos centrados (centroides)
plot(st_geometry(centroides), add = TRUE, pch = 19, cex = 0.7, col = "blue")

# Agregar grilla (opcional)
grid()



############################################
# 2. Definir matrices de vecinos (sin SMR)
############################################

# Vecinos físicos: rook y queen
rook_nb_b  <- nb2listw(poly2nb(columbus, queen = FALSE), style = "B", zero.policy = TRUE)
rook_nb_w  <- nb2listw(poly2nb(columbus, queen = FALSE), style = "W", zero.policy = TRUE)
queen_nb_b <- nb2listw(poly2nb(columbus, queen = TRUE),  style = "B", zero.policy = TRUE)
queen_nb_w <- nb2listw(poly2nb(columbus, queen = TRUE),  style = "W", zero.policy = TRUE)

# Vecinos basados en grafos
trinb         <- tri2nb(xy0)
tri_nb_b      <- nb2listw(trinb, style = "B", zero.policy = TRUE)
tri_nb_w      <- nb2listw(trinb, style = "W", zero.policy = TRUE)

soi_nb_b      <- nb2listw(graph2nb(soi.graph(trinb, xy0)), style = "B", zero.policy = TRUE)
soi_nb_w      <- nb2listw(graph2nb(soi.graph(trinb, xy0)), style = "W", zero.policy = TRUE)

relative_nb_b <- nb2listw(graph2nb(relativeneigh(xy0), sym = TRUE), style = "B", zero.policy = TRUE)
relative_nb_w <- nb2listw(graph2nb(relativeneigh(xy0), sym = TRUE), style = "W", zero.policy = TRUE)

gabriel_nb_b  <- nb2listw(graph2nb(gabrielneigh(xy0), sym = TRUE), style = "B", zero.policy = TRUE)
gabriel_nb_w  <- nb2listw(graph2nb(gabrielneigh(xy0), sym = TRUE), style = "W", zero.policy = TRUE)

# Vecinos por distancia (knn)
knn1_nb_b <- nb2listw(knn2nb(knearneigh(xy0, k = 1)), style = "B", zero.policy = TRUE)
knn1_nb_w <- nb2listw(knn2nb(knearneigh(xy0, k = 1)), style = "W", zero.policy = TRUE)
knn2_nb_b <- nb2listw(knn2nb(knearneigh(xy0, k = 2)), style = "B", zero.policy = TRUE)
knn2_nb_w <- nb2listw(knn2nb(knearneigh(xy0, k = 2)), style = "W", zero.policy = TRUE)
knn3_nb_b <- nb2listw(knn2nb(knearneigh(xy0, k = 3)), style = "B", zero.policy = TRUE)
knn3_nb_w <- nb2listw(knn2nb(knearneigh(xy0, k = 3)), style = "W", zero.policy = TRUE)
knn4_nb_b <- nb2listw(knn2nb(knearneigh(xy0, k = 4)), style = "B", zero.policy = TRUE)
knn4_nb_w <- nb2listw(knn2nb(knearneigh(xy0, k = 4)), style = "W", zero.policy = TRUE)

# List of all matrices
mat <- list(rook_nb_b, rook_nb_w,
            queen_nb_b, queen_nb_w,
            tri_nb_b, tri_nb_w,
            soi_nb_b, soi_nb_w,
            gabriel_nb_b, gabriel_nb_w,
            relative_nb_b, relative_nb_w,
            knn1_nb_b, knn1_nb_w,
            knn2_nb_b, knn2_nb_w,
            knn3_nb_b, knn3_nb_w,
            knn4_nb_b, knn4_nb_w)

# Finding best matrix for CRIME
aux <- numeric(length(mat))
for (i in 1:length(mat)) {
  aux[i] <- moran.test(columbus$CRIME, mat[[i]], alternative = "two.sided")$p.value
}
best_index <- which.min(aux)
best_matrix <- mat[[best_index]] #matriz binaria de vecinos generada por Delaunay triangulation

cat("Best matrix index for CRIME: ", best_index, "\n")
print(moran.test(columbus$CRIME, best_matrix, alternative = "two.sided"))


# Convertir la matriz de vecinos (listw) a una matriz cuadrada clásica
W <- listw2mat(best_matrix)

# Ahora ya tienes W_best como una matriz n × n clásica (por ejemplo 49 × 49)
# Puedes revisarla con
print(W)
dim(W)

Z <- (as.vector(columbus[,9]))$CRIME
WZ <- W %*% Z

# Asegurar que 'columbus_sf' sea un objeto sf válido
columbus_sf <- st_as_sf(columbus)

# Crear una lista para almacenar los data.frames individuales
mapas_df <- list()

# Mapa original: agregar la variable CRIME con etiqueta W0 para mantener el mismo formato
columbus_sf$CRIME_W0 <- columbus$CRIME
mapas_df[[1]] <- columbus_sf %>%
  mutate(W = "Original", valor = CRIME_W0) %>%
  dplyr::select(W, valor, geometry)

# Calcular y agregar medias móviles para cada matriz de pesos según la fórmula:
# μ̂_i = (sum_j w_ij * y_j) / (sum_j w_ij)
for (i in seq_along(mat)) {
  W <- mat[[i]]
  
  # Numerador: suma ponderada de CRIME
  num <- lag.listw(W, columbus$CRIME, zero.policy = TRUE)
  
  # Denominador: suma de pesos por fila
  denom <- sapply(1:length(W$weights), function(i) sum(W$weights[[i]]))
  
  # Media móvil según fórmula dada
  mu_hat <- num / denom
  
  # Crear data.frame con nombre de la matriz
  temp_df <- columbus_sf %>%
    mutate(W = paste0("W", i), valor = mu_hat) %>%
    dplyr::select(W, valor, geometry)
  
  mapas_df[[i + 1]] <- temp_df
}

# Unir todos los data.frames en uno solo
mapas_completos <- do.call(rbind, mapas_df)

# Graficar con ggplot2
ggplot(mapas_completos) +
  geom_sf(aes(fill = valor), color = "white") +
  facet_wrap(~W) +
  scale_fill_viridis_c(option = "C", name = "Media Móvil") +
  labs(title = "Mapas de medias móviles espaciales por matriz W") +
  theme_minimal()

# Variables de interés para el test de Moran
vars_interes <- c("HOVAL", "INC", "OPEN", "PLUMB", "DISCBD")

# Test de Moran para cada variable con diferentes matrices de pesos espaciales
for (var in vars_interes) {
  for (i in 1:length(mat)) {
    aux[i] <- moran.test(columbus[[var]], mat[[i]], alternative = "two.sided")$p.value
  }
  best_i <- which.min(aux)
  cat("Mejor matriz para", var, ":", best_i, "\n")
  print(moran.test(columbus[[var]], mat[[best_i]], alternative = "two.sided"))
}


######################################################################
