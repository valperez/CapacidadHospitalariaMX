rm(list = ls())
library(readr)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(lubridate)
library(dplyr)
library(ggplot2)
library(rstan)
library(zoo)
library(bayestestR)
options(mc.cores = max(parallel::detectCores() - 2, 1))

# Vectores de valores iniciales
fecha_min <- ymd("2020/04/10")
#las fechas seran 10 de julio, 10 de octubre y 10 de enero
vector_fechas <- c(ymd("2020/07/10"), ymd("2020/10/10"), ymd("2021/01/10"))

#Lectura y seleccion de Datos 
githubURL <- ("https://github.com/RodrigoZepeda/CapacidadHospitalariaMX/blob/master/processed/HospitalizacionesMX_estatal.rds?raw=true")
download.file(githubURL,"Hospitalizaciones_estatalMX.rds")

for (i in 1:length(vector_fechas)){

  hospitalizaciones <- readRDS("~/GitHub/CapacidadHospitalariaMX/processed/HospitalizacionesMX_estatal.rds", Encoding("UTF-8"))
  setwd("~/GitHub/CapacidadHospitalariaMX/model")
  stan_fname        <- "BModel.stan"
  
  colnames(hospitalizaciones) <- c("Estado", "Hospitalizados (%)", "Ventilación (%)", "UCI y Ventilación (%)", "Fecha")
  
  #Vamos a pivotear para que me haga una tabla bonita
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
    filter(Fecha > fecha_min & Fecha < vector_fechas[i]) %>% 
    ungroup() %>%
    identity() 
  
  min_fecha <- min(observados$Fecha)
  max_fecha <- max(observados$Fecha)
  
  hospitalizaciones <- observados %>% #Para q sea + rápido
    pivot_wider(names_from = Fecha, values_from = `Hospitalizados (%)`, id_cols = c("Estado")) 
  
  hospitalizaciones_match_estados <- hospitalizaciones %>% select(Estado) %>%
    mutate(EstadoNum = row_number())
  
  
  #Proporción de hospitalizados
  PHosp <- (hospitalizaciones %>% select(-Estado) %>% as.matrix()) 
  
  dias_prediccion <- as.numeric(vector_fechas[i] - ymd("2020/02/05"))
  
  #Caracter?sticas del modelo 
  chains = 4; iter_warmup = 500; nsim = 1000; pchains = 4; m = 7; 
  datos  <- list(m = m, 
                 dias_predict = dias_prediccion,
                 ndias = ncol(PHosp) , nestados = nrow(PHosp), PHosp = PHosp,
                 sigma_mu_hiper = 0.1, sigma_kappa_hiper = 25,
                 mu_mu_hiper = 0, sigma_sigma_hiper = 0.1, sigma_estado_hiper = 0.1) 
  
  #Vamos a intentar con rstan
  sc_model <- stan(file = stan_fname, model_name = "Modelo_beta_covid", iter = nsim, 
                   warmup = iter_warmup, data = datos, chains = chains, seed = 47,
                   control = list(adapt_delta = 0.95))
  
  
  #---------------------------
  #Salvar el modelo
  save(sc_model, file = paste0("ModelBetafit", vector_fechas[i], ".RData"))
  
  #Guardamos las simulaciones por si las dudas
  modelo_ajustado <- summarise_draws(sc_model, median)
  Hosp            <- modelo_ajustado %>% filter(str_detect(variable, "Hosp"))
  
  #analisis <- modelo_ajustado %>% filter(str_detect(variable, "kappa"))
  
  #Calculamos el intervalo de confianza mínima densidad
  for (confint in c(0.5, 0.75, 0.9, 0.95)){
    message(paste0("Calculando el intervalo de ",  100*confint,"%"))
    cime  <- hdi(sc_model, ci = confint)
    cime  <- cime %>% rename(!!sym(paste0("Lower",100*confint)) := CI_low) %>%
      rename(!!sym(paste0("Upper", 100*confint)) := CI_high) %>%
      rename(variable = Parameter)
    Hosp  <- Hosp %>% left_join(cime, by = "variable")
  }
  
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
  save(Hosp, file = paste0("Hosp", vector_fechas[i], ".RData"))
  
  #pdf(paste0("Hosp_predict_",today(),".pdf"), width = 20, height = 10)
  grafica <- ggplot(Hosp, aes(x = Fecha)) +
    geom_ribbon(aes(ymin = `Lower95`, ymax = `Upper95`, fill = "95%")) +
    geom_ribbon(aes(ymin = `Lower90`, ymax = `Upper90`, fill = "90%")) +
    geom_ribbon(aes(ymin = `Lower75`, ymax = `Upper75`, fill = "75%")) +
    geom_ribbon(aes(ymin = `Lower50`, ymax = `Upper50`, fill = "50%")) +
    geom_line(aes(y = `median`, color = "Predichos"), size = 0.1) +
    geom_point(aes(y = `Hospitalizados (%)`, color = "Observados"),
               size = 0.1,data = Hosp %>% filter(Fecha <= !!max_fecha)) +
    facet_wrap(~Estado, ncol = 8) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_date(breaks = "2 months") +
    theme_classic() +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1)) +
    scale_color_manual("Modelo", 
                       values = c("Observados" = "red", 
                                  "Predichos" = "gray75")) +
    scale_fill_manual("Intervalo", 
                      values = c("50%" = "#008b8b",
                                 "75%" = "#00748b",
                                 "90%" = "#005c8b","#00458b"))
}

#Cargamos los model fits
fitjulio <- load("~/GitHub/CapacidadHospitalariaMX/model/ModelBetafit2020-07-10.RData")
sc_model_julio <- sc_model
draws_julio <- as_draws_df(sc_model_julio)


fitoctubre <- load("~/GitHub/CapacidadHospitalariaMX/model/ModelBetafit2020-10-10.RData")
sc_model_octubre <- sc_model
draws_octubre <- as_draws_df(sc_model_octubre)

fitenero <- load("~/GitHub/CapacidadHospitalariaMX/model/ModelBetafit2021-01-10.RData")
sc_model_enero <- sc_model
draws_enero <- as_draws_df(sc_model_enero)

## Vamos a ver la eficiencia de las predicciones

# De los observados, la verdadera proba de los intervalos 75, 90, 95. 
# Es decir, cuantos observados caen en mis intervalos


## ---- Vamos a trabajar con las hospitalizaciones observadas

fecha_min <- ymd("2020/04/10")

hospitalizaciones <- readRDS("~/GitHub/CapacidadHospitalariaMX/processed/HospitalizacionesMX_estatal.rds", Encoding("UTF-8"))
setwd("~/GitHub/CapacidadHospitalariaMX/model")

colnames(hospitalizaciones) <- c("Estado", "Hospitalizados (%)", "Ventilación (%)", "UCI y Ventilación (%)", "Fecha")

#Vamos a pivotear para que me haga una tabla bonita
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
  filter(Fecha > fecha_min) %>%
  ungroup() %>%
  identity()

# Estas son las hospitalizaciones totales de todos los estados para todos los días del
# 11 de Abril 2020 al 5 de Febrero 2021
hospitalizaciones <- observados %>% #Para q sea + rápido
  pivot_wider(names_from = Fecha, values_from = `Hospitalizados (%)`, id_cols = c("Estado"))

hospitalizaciones_match_estados <- hospitalizaciones %>% select(Estado) %>%
  mutate(EstadoNum = row_number())

observados <- observados %>%
  left_join(hospitalizaciones_match_estados, by = "Estado")

hospitalizaciones <- hospitalizaciones %>%
  left_join(hospitalizaciones_match_estados, by = "Estado")

# summary_julio <- summary_julio %>%
#   left_join(hospitalizaciones_match_estados, by = "EstadoNum") %>%
#   left_join(observados, by = c("Estado", "EstadoNum", "Fecha")) %>%
#   rename(Observados = `Hospitalizados (%)`)

## -- Intento 2 ---
# Primero cargamos todo
load("~/GitHub/CapacidadHospitalariaMX/model/Hosp2020-07-10.RData")
Hosp_julio <- Hosp

load("~/GitHub/CapacidadHospitalariaMX/model/Hosp2020-10-10.RData")
Hosp_oct <- Hosp

load("~/GitHub/CapacidadHospitalariaMX/model/Hosp2021-01-10.RData")
Hosp_enero <- Hosp

meses <- c("julio", "oct", "enero")

true_prob <- function(mes){
  assign(paste0("Hosp_", mes), get(paste0("Hosp_", mes))) <- assign(paste0("Hosp_", mes), get(paste0("Hosp_", mes))) %>%
    select("Lower50", "Upper50", 
           "Lower75", "Upper75", 
           "Lower90", "Upper90", 
           "Lower95", "Upper95",
           EstadoNum, DiaNum, Estado, Fecha, 
           `Hospitalizados (%)`) %>%
    left_join(observados, by = c("Estado", "Fecha", "EstadoNum")) %>%
    rename(Hospitalizados_obs = `Hospitalizados (%).y`) %>%
    mutate(dentro_50 = if_else(Hospitalizados_obs >= Lower50 & Hospitalizados_obs <= Upper50, TRUE, FALSE)) %>%
    mutate(dentro_75 = if_else(Hospitalizados_obs >= Lower75 & Hospitalizados_obs <= Upper75,  TRUE, FALSE)) %>%
    mutate(dentro_90 = if_else(Hospitalizados_obs >= Lower90 & Hospitalizados_obs <= Upper90,  TRUE, FALSE)) %>%
    mutate(dentro_95 = if_else(Hospitalizados_obs >= Lower95 & Hospitalizados_obs <= Upper95,  TRUE, FALSE)) 
  
  assign(paste0("Hosp_resultados", mes), get(paste0("Hosp_", mes))) %>%
    summarise(sum(dentro_50, na.rm = TRUE) / nrow(paste0("Hosp_", mes)), 
              sum(dentro_75, na.rm = TRUE) / nrow(paste0("Hosp_", mes)), 
              sum(dentro_90, na.rm = TRUE) / nrow(paste0("Hosp_", mes)), 
              sum(dentro_95, na.rm = TRUE) / nrow(paste0("Hosp_", mes)))
  
}

# Para Julio
dias_prediccion <- as.numeric(vector_fechas[1] - ymd("2020/02/05"))
summary_julio <- summary(sc_model_julio)
#summary son todas las cadenas mrged mientras que c_summary es el summary de cada cadena
summary_julio <- summary_julio$summary

summary_julio <- as.data.frame(summary_julio)

summary_julio <- summary_julio %>%
  mutate(variable = row.names(summary_julio)) %>%
  mutate(EstadoNum = str_extract(variable, "\\[.*,")) %>%
  mutate(DiaNum    = str_extract(variable, ",.*\\]")) %>%
  mutate(EstadoNum = str_remove_all(EstadoNum,"\\[|,")) %>%
  mutate(DiaNum    = str_remove_all(DiaNum,"\\]|,")) %>%
  mutate(EstadoNum = as.numeric(EstadoNum)) %>%
  mutate(DiaNum    = as.numeric(DiaNum)) %>%
  select(-variable) %>%
  filter(!is.na(EstadoNum)) %>%
  filter(!is.na(DiaNum)) %>%
  mutate(Fecha = fecha_min + DiaNum)

summary_julio <- summary_julio %>%
  full_join(Hosp_julio, by = c("EstadoNum", "DiaNum", "Fecha"), na.rm = TRUE)

cor(summary_julio$hosp_obs, summary_julio$mean, method = c( "spearman"))
cor.test(summary_julio$hosp_obs, summary_julio$mean, method=c("spearman"))


# Ahora vamos a obtener la correlacion de Spearman


#dev.off()