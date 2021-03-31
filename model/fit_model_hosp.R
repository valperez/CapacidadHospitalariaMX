rm(list = ls())
set.seed(2342)
library(readr)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(lubridate)
library(dplyr)
library(ggplot2)
library(rstan)
library(zoo)
library(viridis)
library(bayestestR)
library(data.table)
library(R.utils) #install.packages('R.utils')
library(Metrics)
options(mc.cores = max(parallel::detectCores() - 2, 1))

#las fechas seran 10 de julio, 10 de octubre y 10 de enero
vector_fechas <- c(ymd("2020/07/10"), ymd("2020/10/10"), ymd("2021/01/10"))

for (i in 1:length(vector_fechas)){
  #-------------------- Lectura y seleccion de Datos -------------------
  githubURL <- ("https://github.com/RodrigoZepeda/CapacidadHospitalariaMX/blob/master/processed/HospitalizacionesMX_estatal.rds?raw=true")
  download.file(githubURL,"Hospitaizaciones_estatalMX.rds")
  message("Reading database")
  hospitalizaciones <- readRDS("Hospitaizaciones_estatalMX.rds")
  #hospitalizaciones <- readRDS("processed/Hospitaizaciones_estatalMX.rds")
  stan_fname        <- "BModel.stan"
  #stan_fname        <- "model/BModel.stan"
  
  
  #Vamos a pivotear para que me haga una tabla bonita
  message("Processing database")
  observados <- hospitalizaciones %>% 
    select(-`Ventilación (%)`,-`UCI y Ventilación (%)`) %>%
    mutate(`Hospitalizados (%)` = `Hospitalizados (%)`/100) %>%
    mutate(Fecha = ymd(Fecha)) %>%
    arrange(Fecha) %>%
    group_by(Fecha, Estado) %>%
    mutate(`Hospitalizados (%)` = 
             rollapply(`Hospitalizados (%)`, width = 7, 
                       FUN = function(x) mean(x, na.rm=TRUE), partial=TRUE, 
                       fill = NA, align="left")) %>%
    filter(Fecha > ymd("2020/04/10") & Fecha <= vector_fechas[i]) %>% 
    ungroup() %>%
    identity() 
  
  
  min_fecha <- min(observados$Fecha)
  max_fecha <- max(observados$Fecha)
  
  hospitalizaciones <- observados %>% #Para q sea + rápido
    pivot_wider(names_from = Fecha, values_from = `Hospitalizados (%)`, 
                id_cols = c("Estado")) 
  
  hospitalizaciones_match_estados <- hospitalizaciones %>% select(Estado) %>%
    mutate(EstadoNum = row_number())
  
  #Proporción de hospitalizados
  PHosp <- (hospitalizaciones %>% select(-Estado) %>% as.matrix()) 
  
  #Caracter?sticas del modelo 
  chains = 4; iter_warmup = 1000; nsim = 2000; pchains = 4; m = 15; 
  datos  <- list(m = m, 
                 dias_predict = 60,
                 ndias = ncol(PHosp) , nestados = nrow(PHosp), PHosp = PHosp,
                 sigma_mu_hiper = 0.1, sigma_kappa_hiper = 50,
                 mu_mu_temp = 0, sigma_sigma_temp = 10, mu_kappa = 200,
                 mu_mu_hiper = 0, sigma_sigma_hiper = 10, 
                 sigma_estado_hiper = 10) 
  
  # function form 2 with an argument named `chain_id`
  initf2 <- function(chain_id = 1, PHosp = PHosp) {
    list(kappa_raw    = runif(nrow(PHosp), 0, pi/2), 
         alpha_raw    = rnorm(nrow(PHosp), 0, 1),
         lambda_raw   = rep(0, m),
         mu_estado_raw    = rnorm(1, 0, 1),
         sigma_estado_raw = runif(1, 0, pi/2),
         mu_raw           = rnorm(1, 0, 1),
         sigma_raw        = runif(1, 0, pi/2),
         beta_1s_raw      = rnorm(nrow(PHosp),0,1),
         beta_1c_raw      = rnorm(nrow(PHosp),0,1),
         beta_2s_raw      = rnorm(nrow(PHosp),0,1),
         beta_2c_raw      = rnorm(nrow(PHosp),0,1),
         mu_time_raw      = rnorm(1,0,1),
         sigma_time_raw   = runif(1, 0, pi/2)
    )}
  
  # generate a list of lists to specify initial values
  init_ll <- lapply(1:chains, function(id) initf2(chain_id = id, PHosp = PHosp))
  
  #Vamos a intentar con rstan
  message("Fitting model. Go grab a coffee this will take A LOT")
  hosp_model <- cmdstan_model(stan_fname, cpp_options = list(
    cxx_flags = "-O3 -march=native", stan_threads = TRUE
  ))
  
  model_sample <- hosp_model$sample(data = datos, chains = chains, seed = 47, 
                                    iter_warmup = iter_warmup,
                                    adapt_delta = 0.95, iter_sampling = nsim - iter_warmup,
                                    init = init_ll,
                                    output_dir = "cmdstan",
                                    threads_per_chain = 4)
  
  message("Saving results")
  model_sample$save_object(file = paste0("ModelBetafit_fecha", vector_fechas[i],"m", m, ".rds"))
  model_sample$cmdstan_diagnose()
  
  #---------------------------
  #Correr el modelo
  
  #Guardamos las simulaciones por si las dudas
  #modelo_ajustado <- readRDS("model_fit1.rds")
  modelo_ajustado <- summarise_draws(model_sample$draws(), 
                                     ~ quantile(., probs = c(0.005, 0.025, 0.05, 
                                                             0.125, 0.25, 0.325,0.4, 0.5,
                                                             0.6, 0.675,0.75, 0.875, 0.95, 
                                                             0.975, 0.995)))
  Hosp            <- modelo_ajustado %>% filter(str_detect(variable, "Hosp"))
  
  Hosp <- Hosp %>%
    mutate(EstadoNum = str_extract(variable, "\\[.*,")) %>%
    mutate(DiaNum    = str_extract(variable, ",.*\\]")) %>%
    mutate(EstadoNum = str_remove_all(EstadoNum,"\\[|,")) %>%
    mutate(DiaNum    = str_remove_all(DiaNum,"\\]|,")) %>%
    mutate(EstadoNum = as.numeric(EstadoNum)) %>%
    mutate(DiaNum    = as.numeric(DiaNum)) %>%
    select(-variable) %>%
    left_join(hospitalizaciones_match_estados, by = "EstadoNum") %>%
    mutate(Fecha = !!ymd(min_fecha) + DiaNum) %>% 
    full_join(observados, by = c("Fecha","Estado")) %>%
    arrange(Estado, Fecha)
  
  Hosp %>% drop_na(`50%`) %>% 
    write_csv(file = gzfile(paste0("Hosp_", vector_fechas[i],"m", m, ".csv.gz")))
  
  ggplot(Hosp, aes(x = Fecha)) +
    geom_ribbon(aes(ymin = `0.5%`, ymax = `99.5%`, fill = "99%"), alpha = 0.2) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = "95%"), alpha = 0.2) +
    geom_ribbon(aes(ymin = `5%`, ymax = `95%`, fill = "90%"), alpha = 0.2) +
    geom_ribbon(aes(ymin = `12.5%`, ymax = `87.5%`, fill = "75%"), alpha = 0.2) +
    #geom_ribbon(aes(ymin = `25%`, ymax = `75%`, fill = "50%"), alpha = 0.2) +
    geom_line(aes(y = `50%`, color = "Predichos"), size = 0.1) +
    geom_vline(aes(xintercept = today() + 30), linetype = "dashed") +
    geom_point(aes(y = `Hospitalizados (%)`, color = "Observados"), shape = 1,
               size = 0.1,data = Hosp %>% filter(Fecha <= !!max_fecha)) +
    facet_wrap(~Estado, ncol = 8) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_date(breaks = "2 months") +
    theme_classic() +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1)) +
    scale_color_manual("Modelo", 
                       values = c("Observados" = viridis(5)[1], 
                                  "Predichos" = "gray25")) +
    scale_fill_manual("Probabilidad\ndel escenario", 
                      values = c("75%" = viridis(5)[1],
                                 "90%" = viridis(5)[2],
                                 "95%" = viridis(5)[3],
                                 "99%" =  viridis(5)[4])) +
    labs(
      x = "Fecha",
      y = "Capacidad Hospitalaria (%)",
      title = "Escenarios a largo plazo capacidad hospitalaria a partir de la RED-IRAG",
      caption = "*Predicciones después de 30 días (línea vertical) son sólo escenarios",
      subtitle = "Modelo Beta-Bayesiano | Github: @CapacidadHospitalariaMX | Datos de https://www.gits.igg.unam.mx/red-irag-dashboard"
    ) +
    ggsave(paste0("Hosp_predict_v2_",vector_fechas[i],"m", m,".pdf"), 
           width = 20, height = 10)
  
  for (estado in unique(Hosp$Estado)){
    Hosp %>% filter(Estado == !!estado) %>%
      ggplot(aes(x = Fecha)) +
      geom_ribbon(aes(ymin = `0.5%`, ymax = `99.5%`, fill = "99%"), alpha = 0.2) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = "95%"), alpha = 0.2) +
      geom_ribbon(aes(ymin = `5%`, ymax = `95%`, fill = "90%"), alpha = 0.2) +
      geom_ribbon(aes(ymin = `12.5%`, ymax = `87.5%`, fill = "75%"), alpha = 0.2) +
      #geom_ribbon(aes(ymin = `25%`, ymax = `75%`, fill = "50%"), alpha = 0.2) +
      geom_line(aes(y = `50%`, color = "Predichos"), size = 0.1) +
      geom_vline(aes(xintercept = today() + 30), linetype = "dashed") +
      geom_point(aes(y = `Hospitalizados (%)`, color = "Observados"),
                 size = 4,data = Hosp %>% filter(Fecha <= !!max_fecha & Estado == !!estado)) +
      geom_point(aes(y = `Hospitalizados (%)`), color = "white",
                 size = 1,data = Hosp %>% filter(Fecha <= !!max_fecha & Estado == !!estado)) +
      facet_wrap(~Estado, ncol = 8) +
      scale_y_continuous(labels = scales::percent) +
      scale_x_date(breaks = "2 months") +
      theme_classic() +
      theme(axis.text.x  = element_text(angle = 90, hjust = 1)) +
      scale_color_manual("Modelo", 
                         values = c("Observados" = viridis(5)[1], 
                                    "Predichos" = "gray25")) +
      scale_fill_manual("Probabilidad\ndel escenario", 
                        values = c("75%" = viridis(5)[1],
                                   "90%" = viridis(5)[2],
                                   "95%" = viridis(5)[3],
                                   "99%" =  viridis(5)[4])) +
      labs(
        x = "Fecha",
        y = "Capacidad Hospitalaria (%)",
        title = "Escenarios a largo plazo capacidad hospitalaria a partir de la RED-IRAG",
        caption = "*Predicciones después de 30 días (línea vertical) son sólo escenarios",
        subtitle = "Modelo Beta-Bayesiano | Github: @CapacidadHospitalariaMX | Datos de https://www.gits.igg.unam.mx/red-irag-dashboard"
      ) +
      theme(strip.text = element_text(size=50), axis.text = element_text(size = 16),
            axis.title=element_text(size=20,face="bold"),
            plot.title = element_text(size = 20), plot.subtitle = element_text(size = 14),
            plot.caption = element_text(size = 14)) +
      ggsave(paste0("Hosp_predict_v2_",estado,"_",vector_fechas[i],"m", m,".pdf"), width = 20, height = 10)
  }
}
if (file.exists("predictions/PREDICCIONES_HOSP.pdf")){file.remove("predictions/PREDICCIONES_HOSP.pdf")}

# Para obtener las probas verdaderas

# Obtenemos otra vez las hospitalizaciones
githubURL <- ("https://github.com/RodrigoZepeda/CapacidadHospitalariaMX/blob/master/processed/HospitalizacionesMX_estatal.rds?raw=true")
download.file(githubURL,"Hospitaizaciones_estatalMX.rds")
message("Reading database")
hospitalizaciones <- readRDS("Hospitaizaciones_estatalMX.rds")

# Arreglamos para que el data frame esté bonito
observados <- hospitalizaciones %>% 
  select(-`Ventilación (%)`,-`UCI y Ventilación (%)`) %>%
  mutate(`Hospitalizados (%)` = `Hospitalizados (%)`/100) %>%
  mutate(Fecha = ymd(Fecha)) %>%
  arrange(Fecha) %>%
  group_by(Fecha, Estado) %>%
  mutate(`Hospitalizados (%)` = 
           rollapply(`Hospitalizados (%)`, width = 7, 
                     FUN = function(x) mean(x, na.rm=TRUE), partial=TRUE, 
                     fill = NA, align="left")) %>%
  filter(Fecha > ymd("2020/04/10")) %>% 
  ungroup() %>%
  identity() 

#Leemos los resultados del modelo
Hosp_julio <- fread("Hosp_2020-07-10m15.csv.gz")
Hosp_oct <- fread("Hosp_2020-10-10m15.csv.gz")
Hosp_enero <- fread("Hosp_2021-01-10m15.csv.gz")

#Hacemos un vector de meses
meses <- c("julio", "oct", "enero")

true_prob <- function(mes, estado){
  
  assign(paste0("Hosp_", mes, estado), get(paste0("Hosp_", mes)) %>%
           filter(Estado == !!estado) %>%
           select(`50%`,
                  `12.5%`, `87.5%`, 
                  `5%`, `95%`, 
                  `2.5%`, `97.5%`,
                  `0.5%`, `99.5%`,
                  EstadoNum, DiaNum, Estado, Fecha, 
                  `Hospitalizados (%)`) %>%
           mutate(Fecha = as.Date(Fecha)) %>%
           left_join(observados, by = c("Estado", "Fecha")) %>%
           rename(Hospitalizados_obs = `Hospitalizados (%).y`) %>%
           mutate(dentro_75 = if_else(Hospitalizados_obs >= `12.5%` & Hospitalizados_obs <= `87.5%`,  TRUE, FALSE)) %>%
           mutate(dentro_90 = if_else(Hospitalizados_obs >= `5%` & Hospitalizados_obs <=  `95%`,  TRUE, FALSE)) %>%
           mutate(dentro_95 = if_else(Hospitalizados_obs >= `2.5%` & Hospitalizados_obs <= `97.5%`,  TRUE, FALSE)) %>%
           mutate(dentro_99 = if_else(Hospitalizados_obs >= `0.5%` & Hospitalizados_obs <= `99.5%`,  TRUE, FALSE)) )
  
  assign(paste0("Hosp_resultados", mes, estado), get(paste0("Hosp_", mes, estado)) %>%
           summarise(sum(dentro_75, na.rm = TRUE) / nrow(get(paste0("Hosp_", mes, estado))), 
                     sum(dentro_90, na.rm = TRUE) / nrow(get(paste0("Hosp_", mes, estado))), 
                     sum(dentro_95, na.rm = TRUE) / nrow(get(paste0("Hosp_", mes, estado))),
                     sum(dentro_99, na.rm = TRUE) / nrow(get(paste0("Hosp_", mes, estado)))))
  
  lista <- list("Hosp" = get(paste0("Hosp_", mes, estado)), 
                "resultados" = get(paste0("Hosp_resultados", mes, estado)))
  return(lista)
  
}

estados <- c(unique(observados$Estado))

# Usamos la función para calcular para todos los estados 
for (mes in meses){
  for (estado in estados){
    aux <- true_prob(mes, estado)
    assign(paste0("Hosp_", mes, estado), aux$Hosp)
    assign(paste0("Hosp_resultados_", mes, estado), aux$resultados)
  }
}

#Y ahora haremos un for para unir todo
julio <- `Hosp_resultados_julioCiudad de México`
oct <- `Hosp_resultados_octCiudad de México`
enero <- `Hosp_resultados_eneroCiudad de México`

for (mes in meses){
  for (i in 2:length(estados)){
    assign(mes, bind_rows(get(mes), get(paste0("Hosp_resultados_", mes, estados[i]))) )
  }
}

rownames(julio) <- estados
rownames(oct) <- estados
rownames(enero) <- estados

# Función para hacerle bien los nombres de las columnas
nombres_columnas <- function(mes){
  vector_nombres <- c(paste0("Intervalo_75", mes), 
                      paste0("Intervalo_90", mes),
                      paste0("Intervalo_95", mes),
                      paste0("Intervalo_99", mes))
}

colnames(julio) <- nombres_columnas("julio")
colnames(oct) <- nombres_columnas("octubre")
colnames(enero) <- nombres_columnas("enero")

#finalmente bindeamos todo
resultados <- bind_cols(julio, oct, enero)



calculo_mse <- function(mes){
  mse_aux <- data.frame()
  for (estado in estados){
    aux <- get(paste0("Hosp_", mes, estado)) 
    aux <- aux %>%
      select(`50%`, `Hospitalizados_obs`)
    mse_aux <- rbind(mse_aux, mean((aux$`50%` - aux$Hospitalizados_obs)^2, na.rm = T))
  }
  return(mse_aux)
}

mse_julio <- calculo_mse("julio")
mse_oct <- calculo_mse("oct")
mse_enero <- calculo_mse("enero")

mse <- cbind(mse_julio, mse_oct, mse_enero)
colnames(mse) <- meses
rownames(mse) <- estados


calculo_mae <- function(mes){
  mae_aux <- data.frame()
  for (estado in estados){
    aux <- get(paste0("Hosp_", mes, estado)) 
    aux <- aux %>%
      select(`50%`, `Hospitalizados_obs`)
    mae_aux <- rbind(mae_aux, mean(abs(aux$`50%` - aux$Hospitalizados_obs), na.rm = T))
  }
  return(mae_aux)
}

mae_julio <- calculo_mae("julio")
mae_oct <- calculo_mae("oct")
mae_enero <- calculo_mae("enero")

mae <- cbind(mae_julio, mae_oct, mae_enero)
colnames(mae) <- meses
rownames(mae) <- estados

write.csv(resultados, file = "resultados_int.csv", row.names = T)
write.csv(mse, file = "mse.csv", row.names = T)
write.csv(mae, file = "mae.csv", row.names = T)



#Vamos a generar las gráficas que no se generaron y mejorarlas
# OJO: el único que no se genero así nada fue empezar en enero 
# con m=30, 60

for (k in 1:length(vector_fechas)){
  Hosp <-  fread(paste0("Hosp_",vector_fechas[k], "m15.csv.gz"))
  
  ggplot(Hosp, aes(x = Fecha)) +
    geom_ribbon(aes(ymin = `0.5%`, ymax = `99.5%`, fill = "99%"), alpha = 0.2) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = "95%"), alpha = 0.2) +
    geom_ribbon(aes(ymin = `5%`, ymax = `95%`, fill = "90%"), alpha = 0.2) +
    geom_ribbon(aes(ymin = `12.5%`, ymax = `87.5%`, fill = "75%"), alpha = 0.2) +
    geom_ribbon(aes(ymin = `25%`, ymax = `75%`, fill = "50%"), alpha = 0.2) +
    geom_line(aes(y = `50%`, color = "Predichos"), size = 0.1) +
    geom_point(aes(y = `Hospitalizados (%)`, color = "Observados"),
               size = 0.1, data = observados) +
    facet_wrap(~Estado, ncol = 8) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_date(breaks = "2 months") +
    theme_classic() +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1)) +
    scale_color_manual("Modelo", 
                       values = c("Observados" = "red", 
                                  "Predichos" = "gray25")) +
    scale_fill_manual("Probabilidad\ndel escenario", 
                      values = c("50%" = viridis(5)[1],
                                 "75%" = viridis(5)[2],
                                 "90%" = viridis(5)[3],
                                 "95%" = viridis(5)[4],
                                 "99%" =  viridis(5)[5])) +
    labs(
      x = "Fecha",
      y = "Capacidad Hospitalaria (%)",
      title = "Predicciones a largo plazo capacidad hospitalaria a partir de la RED-IRAG", 
      caption = "Datos de https://www.gits.igg.unam.mx/red-irag-dashboard/reviewHome",
      subtitle = vector_fechas[k],
      tag = m
      #subtitle = "Modelo Beta-Bayesiano | Github: @CapacidadHospitalariaMX"
    ) +
    ggsave(paste0("Hosp_predict_v2_",vector_fechas[k],"m", m, "_conObservados.pdf"), width = 20, height = 10)
  dev.off
}
