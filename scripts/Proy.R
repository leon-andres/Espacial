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


# Pasándolo a coordenadas proyectadas
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
vari2 <- variog(pmgs, trend = "cte", estimator.type = "classical", bin.cloud = T, pairs.min = 7)
vari2
plot(vari2)

# Boxplots para los rezagos al cuadrado
par(mfrow = c(2, 5))
for (i in 1:10) {
  boxplot(vari2$bin.cloud[[i]], xlab = paste("rezago", i))
}
par(mfrow = c(1, 1))

# Obtenemos "manualmente" el estimador resistente para datos atípicos
for (i in 1:8) {
  vari2$v[i] <- (median(sqrt(sqrt(vari2$bin.cloud[[i+1]]))) / (2 * 0.457))^4
}
plot(vari2)
eyefit(vari2)

ini1 <- c(4, 4156.58)
fitvar1 <- variofit(vari2,
                    cov.model = "wave",
                    ini1,
                    fix.nugget = TRUE,
                    nugget = 1.72,
                    wei = "npairs")

ini2 <- c(7.91, 14251.14)
fitvar2 <- variofit(vari2,
                    cov.model = "wave",
                    ini2,
                    fix.nugget = FALSE,
                    wei = "npairs")
ini3 <- c(5.16, 5344.18)
fitvar3 <- variofit(vari2,
                    cov.model = "exponential",
                    ini3,
                    fix.nugget = FALSE,
                    wei = "npairs")

plot(vari2$u,vari2$v,
     xlab = "h",
     ylab = "semivarianza",
     cex.lab = 1.3,
     cex.axis = 1.2,
     main = "Estimación teórica del modelo de semivariograma",
     col.main = 4, cex.main =1.3)
lines(fitvar1, col = 1)
lines(fitvar2, col = 2)
lines(fitvar3, col = 3)
legend(130, 18000,
       c("MCO", "MCPnpairs", "MCPcressie"),
       lwd = 2,
       lty = 2:7,
       col = 2:7,
       box.col = 9,
       text.col = 2:7)

summary(fitvar1)
AIC(fitvar1)


# Función para predecir la semivarianza teórica en función de la distancia
predict_semivariogram <- function(model, h) {
  cov.model <- model$cov.model
  cov.pars <- model$cov.pars
  nugget <- model$nugget
  
  # γ(h) = sill - cov(h) + nugget
  sill <- sum(cov.pars, nugget)
  gamma.h <- sill - geoR::cov.spatial(h, cov.model = cov.model, cov.pars = cov.pars)
  return(gamma.h)
}

# Función para calcular el RMSE
calc_rmse <- function(emp, model) {
  pred <- predict_semivariogram(model, emp$u)  # emp$u son distancias
  sqrt(mean((emp$v - pred)^2))                # emp$v son semivarianzas empíricas
}

# Calcular RMSE para cada ajuste
rmse1 <- calc_rmse(vari2, fitvar1)
rmse2 <- calc_rmse(vari2, fitvar2)
rmse3 <- calc_rmse(vari2, fitvar3)

# Mostrar resultados
rmse_comparacion <- c(Wave_Fijo = rmse1, Wave_Libre = rmse2, Exponential = rmse3)
print(rmse_comparacion)

# 1. Graficar el semivariograma empírico
plot(vari2,
     main = "Semivariograma empírico con modelo ajustado",
     xlab = "Distancia (h)",
     ylab = expression(gamma(h)))

# 2. Superponer el modelo teórico ajustado (fitvar1)
lines.variomodel(fitvar1, col = "blue", lwd = 2)

# 3. Guardar el modelo final para interpolación posterior
final_model <- fitvar1



# datno2 <- sqldf("select *
#                 from datos2
#                 where Variable = 'NO2' and Año = 2022 and `Código del Departamento` = 11")

