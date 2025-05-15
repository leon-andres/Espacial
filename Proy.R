setwd("C:/Users/User/Documents/Espacial")

library(readr)
library(sqldf)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(sp)
library(fields)
library(geoR)
library(akima)

datos <- read_csv("CAC.csv", na = "N/A")

datos2 <- datos[,c(1,3:10,12,23:26,28)] 

# Promedio ponderado
datos2<-sqldf("select *, sum(Promedio*[No. de datos])/sum([No. de datos]) as [prom pond], sum([No. de datos]) as TotalDatos
      from datos2
      group by [ID Estacion], Variable, Año
      order by [ID Estacion], variable, Año")

datos2 <- datos2[, !(names(datos2) %in% c("Promedio", "No. de datos"))]

#Data de PM2.5
dat2.5 <- sqldf("select *
                from datos2
                where Variable = 'PM2.5' and Año = 2022 and `Código del Departamento`in (11)")


# Pasándol a coordenadas proyectadas
coords_geo <- SpatialPoints(cbind(dat2.5$Longitud, dat2.5$Latitud),
                            proj4string = CRS("+proj=longlat +datum=WGS84"))
plot(coords_geo, axes = TRUE, main = "Coordenadas Lat-Long")
data_aire_geo <- data.frame(coords_geo,PM2_5 = dat2.5$`prom pond`)
CRS_UTM_CO <- CRS("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")
coords_utm <- spTransform(coords_geo, CRS_UTM_CO)
plot(coords_utm, axes = TRUE, main = "Coordenadas UTM", col = "blue")

#PM2.5 con promedio ponderado y cantidad de datos con las que se hace el promedio (pesos)
PM2_5 <- data.frame(coords_utm,PM2.5 = dat2.5$`prom pond`, TotDat = dat2.5$TotalDatos)
PM2_5 <- PM2_5 %>%
  rename(
    este = coords.x1,
    norte = coords.x2
  )

# Gráficos de la variable vs coordenadas
p1 <- ggplot(PM2_5, aes(x = este, y = PM2.5)) +
  geom_point(color = "blue") +
  labs(title = "PM2.5 vs Este", x = "Este (X)", y = "PM2.5") +
  theme_minimal()
p2 <- ggplot(PM2_5, aes(x = norte, y = PM2.5)) +
  geom_point(color = "blue") +
  labs(title = "PM2.5 vs Norte", x = "Norte (Y)", y = "PM2_5") +
  theme_minimal()
(p1 + p2)

# Correlación contra las variables
cor(PM2_5[, c("este", "norte", "PM2.5")])

# MODELO PARA EXTRAER LA TENDENCIA, TENIENDO COMO PESO LOS TOTALES 
fit1 <- lm(PM2.5 ~ este, data = PM2_5, weights = TotDat)

summary(fit1)
anova(fit1)
res1 <- residuals(fit1)

# Gráficos de los residuales
par(mfrow = c(1, 2))
plot(PM2_5$este, res1)
plot(PM2_5$norte, res1)
par(mfrow = c(1, 1))
boxplot(res1)
PM2_5$residuals <- res1

# Objeto geodata
pmgs <- as.geodata(PM2_5, coords.col = 1:2, data.col = 6, covar.col = 5)
summary(pmgs)

# Gráficos de la geodata
plot(pmgs, qt.col = c("purple",
                          "pink",
                          "green",
                          "yellow"))
plot(pmgs, scatter3d = T)

# Variograma empírico para los residuales usando el estimador clásico
vari2 <- variog(pmgs, trend = "cte", estimator.type = "classical", bin.cloud = T)
vari2
plot(vari2)

# Boxplots para los rezagos al cuadrado
par(mfrow = c(2, 5))
for (i in 1:10) {
  boxplot(vari2$bin.cloud[[i]], xlab = paste("rezago", i))
}
par(mfrow = c(1, 1))

# Obtenemos "manualmente" el estimador resistente para datos atípicos
semivarianza <- numeric(13)  
distancia <- numeric(13)
for (i in 1:13) {
  semivarianza[i] <- (median(sqrt(sqrt(vari2$bin.cloud[[i]]))) / (2 * 0.457))^4
  distancia[i] <- vari2$u[i]
}

plot(distancia,semivarianza)






# datno2 <- sqldf("select *
#                 from datos2
#                 where Variable = 'NO2' and Año = 2022 and `Código del Departamento` = 11")
