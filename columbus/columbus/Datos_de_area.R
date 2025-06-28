setwd("C:\\Users\\User\\Documents\\GitHub\\espacial-2025-1s\\Espacial\\columbus\\columbus")
library(readxl)
library(sf)
library(tmap)
library(ggplot2)


data <- read_xlsx("columbus.xlsx")
shp <- read_sf("columbus.shp")


# Activar modo de visualización estática
tmap_mode("plot")

# Mapa temático de crimen
tm_shape(shp) +
  tm_fill("CRIME", palette = "Reds", style = "quantile", title = "Crimen por vecindario") +
  tm_borders() +
  tm_layout(main.title = "Mapa de Crimen en Columbus, Ohio", legend.outside = TRUE)





ggplot(shp) +
  geom_sf(aes(fill = CRIME), color = "black") +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Crimen") +
  labs(title = "Crimen en vecindarios de Columbus") +
  theme_minimal()
