library(geoR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(sf)

### Carga y ajuste de datos ###
datos <- read.csv("calidad_aire.csv", header = T)

bogota <- datos %>% 
  filter(Código.del.Departamento == 11)

PM25 <- bogota %>%
  filter(Variable == "PM2.5") %>%
  rename(ID = ID.Estacion,
         Anio = Año) %>%
  group_by(ID, Anio, Variable) %>%
  mutate(Valor = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = c("ID", "Latitud", "Longitud", "Anio", "Valor"),
              names_from = "Variable",
              values_from = "Promedio") %>%
  select(-Valor)

PM10 <- bogota %>%
  filter(Variable == "PM10") %>%
  rename(ID = ID.Estacion,
         Anio = Año) %>%
  group_by(ID, Anio, Variable) %>%
  mutate(Valor = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = c("ID", "Latitud", "Longitud", "Anio", "Valor"),
              names_from = "Variable",
              values_from = "Promedio") %>%
  select(-Valor)

### Análisis Descriptivo (Año 2022) ###
### Variable PM2.5
PM25_2021 <- PM25 %>%
  filter(Anio == 2021)

PM10_2021 <- PM10 %>%
  filter(Anio == 2021)

#### Mapas variables PM2.5 y PM10
mapa <- st_read("mapa/MANZANA.shp")

ggplot(mapa) +
  geom_sf(fill = "lightgray", color = "black") +
  geom_point(data = PM25_2021, 
             aes(x = Longitud, 
                 y = Latitud, 
                 color = PM2.5, 
                 size = PM2.5), 
             alpha = 0.7) +
  scale_color_viridis_c(option = "plasma") +
  labs(title = "Particulas PM10 en Bogotá año 2021",
       x = "Longitud",
       y = "Latitud",
       color = "Valor",
       size = "Valor") +
  theme_minimal()

ggplot(mapa) +
  geom_sf(fill = "lightgray", color = "black") +
  geom_point(data = PM25_2021,
             aes(x = Longitud,
                 y = Latitud,
                 color = PM2.5,
                 size = PM2.5), 
             alpha = 0.7,
             show.legend = c(color = TRUE, size = FALSE))
scale_color_viridis_c(option = "plasma", name = "Valor") +
  labs(title = "Partículas PM2.5 en Bogotá año 2021",
       x = "Longitud",
       y = "Latitud") +
  coord_sf(xlim = c(-74.205, -74.04),
           ylim = c(4.4, 4.9)) +
  theme_minimal()

### Modelo de regresión para PM2.5
summary(lm(PM2.5 ~ Longitud + Latitud + I(Longitud*Latitud), data = PM25_2021))
summary(lm(PM2.5 ~ Longitud, data = PM25_2021))
summary(lm(PM2.5 ~ Longitud + I(Longitud*Latitud), data = PM25_2021))

# Es el modelo con mejor R^2 ajustado
summary(lm(PM2.5 ~ Longitud + I(Longitud^2), data = PM25_2021))
lm_PM2.5 <- lm(PM2.5 ~ Longitud + I(Longitud^2), data = PM25_2021)
plot(lm_PM2.5$residuals)
PM25_2021$residuales <- lm_PM2.5$residuals

### Modelo de regresión para PM10
summary(lm(PM10 ~ Longitud + Latitud + I(Longitud*Latitud), data = PM10_2021))
summary(lm(PM10 ~ Longitud + Latitud, data = PM10_2021))
summary(lm(PM10 ~ Longitud, data = PM10_2021))

# Es el modelo con mejor R^2 ajustado
summary(lm(PM10 ~ Longitud + I(Longitud^2), data = PM10_2021))
lm_PM10 <- lm(PM10 ~ Longitud + I(Longitud^2), data = PM10_2021)

PM10_2021$residuales <- lm_PM10$residuals

# Datos centrados 
PM25_2021resg <- as.geodata(PM25_2021, coords.col = 2:3, data.col = 6)
plot.geodata(PM25_2021resg, scatter3d = TRUE)

PM10_2021resg <- as.geodata(PM10_2021, coords.col = 2:3, data.col = 6)
plot.geodata(PM10_2021resg, scatter3d = TRUE)

# Variograma variable PM2.5
variog(PM25_2021resg)
plot(variog(PM25_2021resg))

# Variograma PM10
variPM10_2021 <- variog(PM10_2021resg)
plot(variog(PM10_2021resg))

expPM10_2021 <- variofit(variPM10_2021,
                         ini.cov.pars = c(51.68, 0.1),
                         cov.model = "exponential")

