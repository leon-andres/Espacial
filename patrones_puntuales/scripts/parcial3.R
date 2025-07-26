library(sf)
library(readr)
library(spatstat)
library(dplyr)
library(ggplot2)
library(lubridate)

# Cargar de datos
datos <- read_csv("data/NYCAccidents2020.csv")

# Eliminar datos con ausencia de coordenadas y posible error de registro
datos_sf <- datos %>%
  filter(!is.na(LONGITUDE), !is.na(LATITUDE), LONGITUDE < -73.5, 
         LONGITUDE > -74.5, LATITUDE > 40.3, LATITUDE < 41.0) %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = FALSE)

# Se extrae una muestra de 1000 individuos
set.seed(123)
datos_sf <- sample_n(datos_sf, 1000)

# Shapefile de NYC
mapa <- read_sf("data/nybb.shp")

# Transformar coordenadas a UTM
datos_utm <- st_transform(datos_sf, crs = 26918)
mapa_utm  <- st_transform(mapa,  crs = 26918)

# Ventana de observación
bbox <- st_bbox(datos_utm)
win <- owin(xrange = c(bbox["xmin"], bbox["xmax"]),
            yrange = c(bbox["ymin"], bbox["ymax"]))

# Crear patrón de puntos
coords <- st_coordinates(datos_utm)
accidentes <- ppp(x = coords[,1], y = coords[,2], window = win)

# Gráfica de puntos
ggplot() +
  geom_sf(data = datos_utm, aes(), color = "black", size = 0.3, alpha = 0.5) +
  coord_sf(xlim = st_bbox(mapa_utm)[c("xmin", "xmax")],
           ylim = st_bbox(mapa_utm)[c("ymin", "ymax")],
           expand = FALSE) +
  labs(title = "Mapa de puntos - Accidentes en NYC") +
  theme_minimal()

# Mapa de accidente en NYC
ggplot() +
  geom_sf(data = mapa_utm, fill = NA, color = "blue", size = 1) +
  geom_sf(data = datos_utm, aes(), color = "black", size = 0.3, alpha = 0.5) +
  coord_sf(xlim = st_bbox(mapa_utm)[c("xmin", "xmax")],
           ylim = st_bbox(mapa_utm)[c("ymin", "ymax")],
           expand = FALSE) +
  labs(title = "Mapa de Accidentes en NYC") +
  theme_minimal()

# Prueba de aleatoriedad basada en cuadrantes
# Cuadrante 4x4
qcount4x4 <- quadratcount(accidentes, nx = 4, ny = 4)
plot(qcount4x4, main = "Cuadricula 4x4", cex = 0.5)
plot(accidentes, add = TRUE, pch = ".", cex = 0.3)
plot(mapa_utm$geometry, add = TRUE, border = "black", lwd = 1)
quadrat.test(accidentes, nx = 4, ny = 4)

# Cuadrante 6x6
qcount6x6 <- quadratcount(accidentes, nx = 6, ny = 6)
plot(qcount6x6, main = "Cuadricula 6x6", cex = 0.5)
plot(accidentes, add = TRUE, pch = ".", cex = 0.3)
plot(mapa_utm$geometry, add = TRUE, border = "black", lwd = 1)
quadrat.test(accidentes, nx = 6, ny = 6)

# Cuadrante 10x10
qcount10x10 <- quadratcount(accidentes, nx = 10, ny = 10)
plot(qcount10x10, main = "Cuadricula 10x10", cex = 0.5)
plot(accidentes, add = TRUE, pch = ".", cex = 0.3)
plot(mapa_utm$geometry, add = TRUE, border = "black", lwd = 1)
quadrat.test(accidentes, nx = 10, ny = 10)

# Pruebas de aleatoriedad basada en distancias
# Función G
env_G <- envelope(accidentes, Gest, nsim = 1000, correction = "border")
plot(env_G, main = "Función G")

# Test de Hopkins-Skellam
hopskel.test(accidentes, method = "asymptotic")

# Indice de Clark-Evans
clarkevans(accidentes, correction = "none")

# Estimación de intensidad
# Estimación de intensidad no paramétrica
est_no_par <- density.ppp(accidentes, sigma = bw.diggle(accidentes))
plot(est_no_par, main = "Estimación de intensidad no parámetrica")
plot(mapa_utm$geometry, add = TRUE, border = "black")

# Estimación de intensidad paramétrica
# Modelo media constante
modelo_homog <- ppm(accidentes, ~1)
summary(modelo_homog)

# Modelo media no constante
modelo_inhomog <- ppm(accidentes, ~x + y)
summary(modelo_inhomog)

# Pruebas para detectar interacción
# Función K de Ripley
# env_K <- envelope(accidentes, Kest, nsim = 100, correction = "border")
env_K <- envelope(accidentes, Kinhom, nsim = 100, correction = "border")
plot(env_K, main = "Función K de Ripley")

# Función L de Besag
# env_L <- envelope(accidentes, Lest, nsim = 100, correction = "border")
env_L <- envelope(accidentes, Linhom, nsim = 100, correction = "border")
plot(env_L, main = "Función L de Besag")

# Función de correlación por pares (g)
# env_G <- envelope(accidentes, pcf, nsim = 100, correction = "border")
env_PCF <- envelope(accidentes, pcfinhom, nsim = 100, correction = "Ripley")
plot(env_PCF, main = "Función de correlación por pares (g)")

# Patrón puntual espacial con marca
datos_utm <- datos_utm %>%
  mutate(heridos_cat = cut(`NUMBER OF PERSONS INJURED`,
                           breaks = c(0, 1, 2, Inf),
                           labels = c("0", "1-2", "3+"),
                           right = FALSE))

table(datos_utm$heridos_cat)

ggplot() +
  geom_sf(data = mapa_utm, fill = NA, color = "blue", size = 1) +
  geom_sf(data = datos_utm %>% filter(!is.na(heridos_cat)),
          aes(color = heridos_cat),
          size = 0.3, alpha = 0.7) +
  scale_color_brewer(palette = "Dark2", name = "No. de heridos") +
  coord_sf(xlim = st_bbox(mapa_utm)[c("xmin", "xmax")],
           ylim = st_bbox(mapa_utm)[c("ymin", "ymax")],
           expand = FALSE) +
  labs(title = "Mapa de Accidentes en NYC según No. de Heridos") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14))

# Crear patrón de puntos con marca (numero de personas heridas)
accidentes_marca <- ppp(x = coords[,1], y = coords[,2], window = win,
                        marks = datos_utm$heridos_cat)

accidentes_split <- split(accidentes_marca)

plot(accidentes_split, use.marks=FALSE,
     main="Accidentes de transito en NYC",
     pch=21, bg=2, cex=0.3)

summary(accidentes_marca)

# Cuadrante 4x4
par(mfrow = c(1, 3), mar = c(2, 2, 2, 1))
for (cat in names(accidentes_split)) {
  subpat <- accidentes_split[[cat]]
  q <- quadratcount(subpat, nx = 4, ny = 4) # Calcular cuadrante
  plot(q, main = paste("Marca:", cat), cex = 0.6)
  plot(subpat, add = TRUE, pch = ".", cex = 0.3)
  plot(mapa_utm$geometry, add = TRUE, border = "black", lwd = 1)
}
par(mfrow = c(1, 1))
quadrat.test.splitppp(accidentes_split, nx = 4, ny = 4)

# Cuadrante 6x6
par(mfrow = c(1, 3), mar = c(2, 2, 2, 1))
for (cat in names(accidentes_split)) {
  subpat <- accidentes_split[[cat]]
  q <- quadratcount(subpat, nx = 6, ny = 6) # Calcular cuadrante
  plot(q, main = paste("Marca:", cat), cex = 0.6)
  plot(subpat, add = TRUE, pch = ".", cex = 0.3)
  plot(mapa_utm$geometry, add = TRUE, border = "black", lwd = 1)
}
par(mfrow = c(1, 1))
quadrat.test.splitppp(accidentes_split, nx = 6, ny = 6)

# Cuadrante 10x10
par(mfrow = c(1, 3), mar = c(2, 2, 2, 1))
for (cat in names(accidentes_split)) {
  subpat <- accidentes_split[[cat]]
  q <- quadratcount(subpat, nx = 10, ny = 10) # Calcular cuadrante
  plot(q, main = paste("Marca:", cat), cex = 0.6)
  plot(subpat, add = TRUE, pch = ".", cex = 0.3)
  plot(mapa_utm$geometry, add = TRUE, border = "black", lwd = 1)
}
par(mfrow = c(1, 1))
quadrat.test.splitppp(accidentes_split, nx = 10, ny = 10)

# Pruebas de aleatoriedad basada en distancias
# Función G entre marcas
plot(alltypes(accidentes_marca, "G"))
# Función G de punto marca i a cualquier otra
plot(alltypes(accidentes_marca, "Gdot"))

# Estimación de intensidad
# Estimación de intensidad no paramétrica
est_no_par_mark <- density.splitppp(accidentes_split, 
                                    sigma = bw.diggle(accidentes))
plot(est_no_par_mark, main = "Estimación de intensidad no parámetrica")

# Estimación de intensidad paramétrica
# Modelo media constante
modelo_homog_mark <- ppm(accidentes_marca, ~1)
summary(modelo_homog_mark)

# Modelo media no constante
modelo_inhomog_mark <- ppm(accidentes_marca, ~marks)
summary(modelo_inhomog_mark)

# Pruebas para detectar interacción
# Función K de Ripley entre marcas
plot(alltypes(accidentes_marca, "Kcross.inhom"))
# Función K de Ripley de punto marca i a cualquier otra
plot(alltypes(accidentes_marca, "Kdot.inhom"))
# Función L de Besag entre marcas
plot(alltypes(accidentes_marca, "Lcross.inhom"))
# Función L de Besag de punto marca i a cualquier otra
plot(alltypes(accidentes_marca, "Ldot.inhom"))
# Función de correlación por pares (g) entre marcas
plot(alltypes(accidentes_marca, "pcfcross.inhom"))
# Función de correlación por pares (g) de punto marca i a cualquier otra
plot(alltypes(accidentes_marca, "pcfdot.inhom"))

datos_utm <- datos_utm %>%
  mutate(FECHA = ymd(`CRASH DATE`)) %>%
  filter(!is.na(FECHA)) %>%
  mutate(cuatrimestre = case_when(
    month(FECHA) %in% 1:4  ~ "Cuatri 1",
    month(FECHA) %in% 5:8  ~ "Cuatri 2",
    month(FECHA) %in% 9:12 ~ "Cuatri 3",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(cuatrimestre))

cuatri_split <- split(datos_utm, datos_utm$cuatrimestre)

analisis_estacion <- function(datos_sub, mapa_utm) {
  coords <- st_coordinates(datos_sub)
  win <- as.owin(st_geometry(mapa_utm))
  accidentes <- ppp(x = coords[,1], y = coords[,2], window = win)
  
  # Pruebas de homogeneidad basadas en cuadrantes
  for (n in c(4, 6, 10)) {
    qcount <- quadratcount(accidentes, nx = n, ny = n)
    plot(qcount, main = paste0("Cuadrícula ", n, "x", n))
    plot(accidentes, add = TRUE, pch = ".", cex = 0.3)
    plot(mapa_utm$geometry, add = TRUE, border = "black", lwd = 1)
    print(quadrat.test(accidentes, nx = n, ny = n))
  }
  
  # Prueba basada en distancias: Función G con corrección borde
  env_G <- envelope(accidentes, Gest, nsim = 99, correction = "border")
  plot(env_G, main = "Función G")
  
  # Estimación de intensidad no paramétrica
  sigma <- bw.diggle(accidentes)
  est_no_par <- density.ppp(accidentes, sigma = sigma)
  plot(est_no_par, main = paste0("Intensidad no paramétrica (sigma = ", round(sigma,1), ")"))
  plot(mapa_utm$geometry, add = TRUE, border = "black")
  
  # Estimación de intensidad paramétrica
  modelo_homog <- ppm(accidentes, ~1)
  print(summary(modelo_homog))
  
  modelo_inhomog <- ppm(accidentes, ~x + y)
  print(summary(modelo_inhomog))
  
  # Pruebas para detectar interacción: Función K, L y g (corrección borde)
  env_K <- envelope(accidentes, Kinhom, nsim = 99, correction = "border")
  plot(env_K, main = "Función K de Ripley")
  
  env_L <- envelope(accidentes, Linhom, nsim = 99, correction = "border")
  plot(env_L, main = "Función L de Besag")
  
  env_g <- envelope(accidentes, pcfinhom, nsim = 99, correction = "border")
  plot(env_g, main = "Función de correlación por pares (g)")
}

# Ejecutar análisis para cada cuatrimestre
for (esta in names(estacion_split)) {
  cat("\n\n### Análisis para:", esta, "###\n")
  analisis_estacion(esta[[esta]], mapa_utm)
}
