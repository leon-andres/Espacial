################################################################
#### Punto 2: Modelo lineal de corregionalización sin Gstat ####
################################################################

# Variogramas directos y cruzados
calcular_variogramas <- function(datos, cutoff, width) {
  n_puntos <- nrow(datos)
  coords <- datos@coords
  var1 <- datos@data$PM10
  var2 <- datos@data$PM2.5
  
  # Matriz de distancias
  matriz_dist <- as.matrix(dist(coords))
  
  # Limites
  lags_limites <- seq(0, cutoff, by = width)
  n_lags <- length(lags_limites) - 1
  
  # Resultados por lag
  acumulador <- lapply(1:n_lags, function(i) {
    list(
      suma_diff_dir1 = 0,
      suma_diff_dir2 = 0,
      suma_prod_cruz = 0,
      num_pares = 0
    )
  })
  
  for (i in 1:(n_puntos - 1)) {
    for (j in (i + 1):n_puntos) {
      
      # Distancia
      distancia_ij <- matriz_dist[i, j]
      
      # Cutoff
      if (distancia_ij > cutoff) {
        next
      }
      
      # Pertenencia lag
      lag_idx <- findInterval(distancia_ij, lags_limites)
      
      if (lag_idx == 0 || lag_idx > n_lags) {
        next
      }
      
      # Diferencias al cuadrado y el producto
      diff1 <- var1[i] - var1[j]
      diff2 <- var2[i] - var2[j]
      
      sq_diff_dir1 <- diff1^2
      sq_diff_dir2 <- diff2^2
      prod_diff_cruz <- diff1 * diff2
      
      # Valores de los lags
      acumulador[[lag_idx]]$suma_diff_dir1 <- acumulador[[lag_idx]]$suma_diff_dir1 + sq_diff_dir1
      acumulador[[lag_idx]]$suma_diff_dir2 <- acumulador[[lag_idx]]$suma_diff_dir2 + sq_diff_dir2
      acumulador[[lag_idx]]$suma_prod_cruz <- acumulador[[lag_idx]]$suma_prod_cruz + prod_diff_cruz
      acumulador[[lag_idx]]$num_pares <- acumulador[[lag_idx]]$num_pares + 1
    }
  }
  
  # Calculo de valores
  resultado_final <- lapply(1:n_lags, function(k) {
    lag_data <- acumulador[[k]]
    
    if (lag_data$num_pares == 0) {
      return(NULL)
    }
    
    # Estimador de Matheron
    gamma_dir1 <- lag_data$suma_diff_dir1 / (2 * lag_data$num_pares)
    gamma_dir2 <- lag_data$suma_diff_dir2 / (2 * lag_data$num_pares)
    gamma_cruz <- lag_data$suma_prod_cruz / (2 * lag_data$num_pares)
    
    # Distancia representativa
    dist_media <- (lags_limites[k] + lags_limites[k+1]) / 2
    
    data.frame(
      dist = dist_media,
      gamma_dir1 = gamma_dir1,
      gamma_dir2 = gamma_dir2,
      gamma_cruz = gamma_cruz,
      np = lag_data$num_pares
    )
  })
  
  # Combinar los resultados
  resultado_df <- do.call(rbind, resultado_final)
  
  output_list <- list(
    directo1 = data.frame(dist = resultado_df$dist, 
                          gamma = resultado_df$gamma_dir1, 
                          np = resultado_df$np),
    directo2 = data.frame(dist = resultado_df$dist, 
                          gamma = resultado_df$gamma_dir2, 
                          np = resultado_df$np),
    cruzado  = data.frame(dist = resultado_df$dist, 
                          gamma = resultado_df$gamma_cruz, 
                          np = resultado_df$np)
  )
  
  return(output_list)
}

# Semivariogramas directos y cruzados
variogramases <- calcular_variogramas(
  datos = df_ok,
  cutoff = 15707.75,   
  width = 1208.288     
)

# # Visualizar los resultados
par(mfrow=c(2,2), mar = c(4, 4, 2, 1))
plot(variogramases$directo1$dist, variogramases$directo1$gamma, 
     xlab = "Distancia", ylab = "Semivarianza", main = "PM10", pch = 16,
     ylim = c(min(0, min(variogramases$directo1$gamma)), 
              max(variogramases$directo1$gamma) * 1.1))
plot(variogramases$directo2$dist, variogramases$directo2$gamma, 
     xlab = "Distancia", ylab = "Semivarianza", main = "PM2.5", pch = 16,
     ylim = c(min(0, min(variogramases$directo2$gamma)), 
              max(variogramases$directo2$gamma) * 1.1))
plot(variogramases$cruzado$dist, variogramases$cruzado$gamma, 
     xlab = "Distancia", ylab = "Semivarianza Cruzada", main = "Cruzado",
     pch = 16, ylim = c(min(0, min(variogramases$cruzado$gamma)), 
                        max(variogramases$cruzado$gamma) * 1.1))


modelo_lmc_gauss <- function(params, h) {
  b11_gau <- params[1]
  b12_gau <- params[2]
  b22_gau <- params[3]
  range_gau <- params[4]
  
  # Validez del modelo
  valido_gau <- (b11_gau>=0)&&(b22_gau>=0)&&((b11_gau*b22_gau-b12_gau^2) >= 0)
  
  if (!valido_gau) {
    return(list(gamma11 = rep(Inf, length(h)), 
                gamma22 = rep(Inf, length(h)), 
                gamma12 = rep(Inf, length(h))))
  }
  
  # Modelo Gaussiano
  cov_gau <- geoR::cov.spatial(h, cov.model = "gaussian", 
                               cov.pars = c(1, range_gau))
  gamma_gau <- 1 - cov_gau
  
  # semivariogramas teóricos
  gamma11_teorico <- b11_gau * gamma_gau
  gamma22_teorico <- b22_gau * gamma_gau
  gamma12_teorico <- b12_gau * gamma_gau
  
  return(list(gamma11 = gamma11_teorico, 
              gamma22 = gamma22_teorico, 
              gamma12 = gamma12_teorico))
}


valores_iniciales_gau <- c(
  b11_gau = 0.20,     
  b12_gau = 0.10,     
  b22_gau = 0.20,    
  range_gau = 8000
)

objetivo_wls <- function(params, variogramas_empiricos, modelo_func) {
  v_dir1_emp <- variogramas_empiricos$directo1
  v_dir2_emp <- variogramas_empiricos$directo2
  v_cruz_emp <- variogramas_empiricos$cruzado
  
  gamma_teorico <- modelo_func(params, h = v_dir1_emp$dist)
  
  if (is.infinite(gamma_teorico$gamma11[1])) {
    return(1e10)
  }
  
  error_sq_dir1 <- (v_dir1_emp$gamma - gamma_teorico$gamma11)^2
  error_sq_dir2 <- (v_dir2_emp$gamma - gamma_teorico$gamma22)^2
  error_sq_cruz <- (v_cruz_emp$gamma - gamma_teorico$gamma12)^2
  
  peso_dir1 <- v_dir1_emp$np
  peso_dir2 <- v_dir2_emp$np
  peso_cruz <- v_cruz_emp$np
  
  wsse <- sum(peso_dir1 * error_sq_dir1) + 
    sum(peso_dir2 * error_sq_dir2) + 
    sum(peso_cruz * error_sq_cruz)
  
  return(wsse)
}


ajustar_lmc <- function(variogramas_empiricos, 
                        valores_iniciales, modelo_func) {
  num_params <- length(valores_iniciales)
  
  limites_inf <- rep(0, num_params)
  limites_inf[c(2)] <- -Inf
  limites_sup <- rep(Inf, num_params)
  
  resultado_opt <- optim(
    par = valores_iniciales,
    fn = objetivo_wls,
    variogramas_empiricos = variogramas_empiricos,
    modelo_func = modelo_func,
    method = "L-BFGS-B",
    lower = limites_inf,
    upper = limites_sup,
    control = list(trace = 1, maxit = 1000)
  )
  
  return(resultado_opt)
}


ajuste_gauss <- ajustar_lmc(
  variogramas_empiricos = variogramases,
  valores_iniciales = valores_iniciales_gau,
  modelo_func = modelo_lmc_gauss
)

parametros_optimos_gau <- ajuste_gauss$par

v_pm10_emp <- variogramases$directo1
v_pm25_emp <- variogramases$directo2
v_cruzado_emp <- variogramases$cruzado

dist_teorica <- seq(0, 22000, length.out = 200)
curvas_ajustadas_puro <- modelo_lmc_gauss(parametros_optimos_gau, 
                                          h = dist_teorica)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

plot(v_pm10_emp$dist, v_pm10_emp$gamma, 
     xlab = "Distancia", ylab = "Semivarianza",
     main = "PM10", pch = 16,
     ylim = c(0, max(v_pm10_emp$gamma) * 1.1))
grid()
lines(dist_teorica, curvas_ajustadas_puro$gamma11, lwd = 2)

plot(v_pm25_emp$dist, v_pm25_emp$gamma, 
     xlab = "Distancia", ylab = "Semivarianza",
     main = "PM2.5", pch = 16,
     ylim = c(0, max(v_pm25_emp$gamma) * 1.1))
grid()
lines(dist_teorica, curvas_ajustadas_puro$gamma22, lwd = 2)

plot(v_cruzado_emp$dist, v_cruzado_emp$gamma,
     xlab = "Distancia", ylab = "Semivarianza Cruzada",
     main = "Cruzado", pch = 16,
     ylim = c(min(0, min(v_cruzado_emp$gamma)), max(v_cruzado_emp$gamma) * 1.1))
grid()
lines(dist_teorica, curvas_ajustadas_puro$gamma12, lwd = 2)

par(mfrow = c(1, 1))

################################################################
############ Punto 6: Método de pseudoverosimilitud ############
################################################################

library(readxl)
library(dplyr)
library(lubridate)

pm25_data <- read_excel("PM2_5.xlsx")

datos_st <- pm25_data %>%
  mutate(
    fecha = ymd(fecha),
    tiempo = as.numeric(fecha - min(fecha, na.rm = TRUE))
  ) %>%
  dplyr::select(este, norte, tiempo, PM2.5, estacion) %>%
  na.omit()

modelo_st_separable <- function(params, h, u) {
  nugget <- params[1]; sill_esp <- params[2]; range_esp <- params[3]
  sill_temp <- params[4]; range_temp <- params[5]
  
  # Componente espacial - modelo esférico
  gamma_esp <- ifelse(h > 0, sill_esp * ((1.5 * h / range_esp) - 0.5 * (h / range_esp)^3), 0)
  gamma_esp[h > range_esp] <- sill_esp
  
  # Componente temporal - modelo esférico
  gamma_temp <- ifelse(u > 0, sill_temp * ((1.5 * u / range_temp) - 0.5 * (u / range_temp)^3), 0)
  gamma_temp[u > range_temp] <- sill_temp
  
  return(nugget + gamma_esp + gamma_temp)
}

# Función de pseudoverosimilitud
#Negativo de la log-verosimilitud compuesta
funcion_objetivo_cl <- function(params, datos, modelo_st_func, cutoff_esp, cutoff_temp) {
  
  n_puntos <- nrow(datos)
  nll_total <- 0
  
  if (n_puntos < 2) {
    return(0)
  }
  
  for (i in seq_len(n_puntos - 1)) {
    for (j in (i + 1):n_puntos) {
      
      # Separación espacial (h) y temporal (u)
      h <- sqrt((datos$este[i] - datos$este[j])^2 + (datos$norte[i] - datos$norte[j])^2)
      u <- abs(datos$tiempo[i] - datos$tiempo[j])
      
      if (h > cutoff_esp || u > cutoff_temp) {
        next
      }
      
      V_sq <- (datos$valor[i] - datos$valor[j])^2
      gamma_teorico <- modelo_st_func(params, h, u)
      
      if (is.na(gamma_teorico) || gamma_teorico <= 1e-9) next
      
      nll_par <- log(gamma_teorico) + V_sq / (2 * gamma_teorico)
      nll_total <- nll_total + nll_par
    }
  }
  return(nll_total)
}


# Función de ajuste
ajustar_modelo_cl <- function(datos, valores_iniciales, cutoff_esp, cutoff_temp) {
  limites_inf <- rep(1e-6, length(valores_iniciales))
  resultado <- optim(
    par = valores_iniciales, fn = funcion_objetivo_cl, datos = datos,
    modelo_st_func = modelo_st_separable, cutoff_esp = cutoff_esp,
    cutoff_temp = cutoff_temp, method = "L-BFGS-B", lower = limites_inf,
    control = list(trace = 1, fnscale = 1e4)
  )
  return(resultado)
}

# Valores iniciales
valores_iniciales <- c(nugget=2, sill_esp=10, range_esp=15000, sill_temp=5, range_temp=60)

#ajuste_st <- ajustar_modelo_cl(
#    datos = datos_st,
#    valores_iniciales = valores_iniciales,
#    cutoff_esp = 25000, 
#    cutoff_temp = 180)

#parametros_estimados <- ajuste_st$par
#names(parametros_estimados) <- names(valores_iniciales)
#print(round(parametros_estimados, 4))

set.seed(123)
datos_prueba <- datos_st[sample(nrow(datos_st), 200), ]
ajuste_prueba <- ajustar_modelo_cl(
  datos = datos_prueba,
  valores_iniciales = valores_iniciales,
  cutoff_esp = 25000, 
  cutoff_temp = 180
)