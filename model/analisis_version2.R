rm(list = ls())
library(readr)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(lubridate)
library(dplyr)
library(ggplot2)
library(rstan)

#-------------------- Lectura y seleccion de Datos -------------------
#githubURL <- ("https://github.com/RodrigoZepeda/CapacidadHospitalariaMX/blob/master/processed/HospitalizacionesMX_estatal.rds?raw=true")
#download.file(githubURL,"Hospitaizaciones_estatalMX.rds")
hospitalizaciones <- readRDS("Hospitaizaciones_estatalMX.rds")
stan_fname        <- "~/Dropbox/HospitalizacionesCOVIDMX/model/BModel.stan"

#Vamos a pivotear para que me haga una tabla bonita
observados <- hospitalizaciones %>% 
  select(-`Ventilación (%)`,-`UCI y Ventilación (%)`) %>%
  mutate(`Hospitalizados (%)` = `Hospitalizados (%)`/100) %>%
  mutate(Fecha = ymd(Fecha)) %>%
  arrange(Fecha) %>%
  filter(Fecha > ymd("2020/04/10")) %>% #Para que ya no tenga NA's
  #filter(Estado %in% sample(unique(hospitalizaciones$Estado), 5)) %>%
  identity() 

min_fecha <- min(observados$Fecha)
max_fecha <- max(observados$Fecha)

hospitalizaciones <- observados %>% #Para q sea + rápido
  pivot_wider(names_from = Fecha, values_from = `Hospitalizados (%)`, id_cols = c("Estado")) 

hospitalizaciones_match_estados <- hospitalizaciones %>% select(Estado) %>%
  mutate(EstadoNum = row_number())

#Proporción de hospitalizados
PHosp <- (hospitalizaciones %>% select(-Estado) %>% as.matrix()) 

#Caracter?sticas del modelo 
chains = 4; iter_warmup = 500; nsim = 500; pchains = 4; m = 7; 
datos  <- list(m = m, 
               dias_predict = 30,
               ndias = ncol(PHosp) , nestados = nrow(PHosp), PHosp = PHosp,
               sigma_mu_hiper = 0.1, sigma_kappa_hiper = 25,
               mu_mu_hiper = 0, sigma_sigma_hiper = 0.1, sigma_estado_hiper = 0.1) 

#Vamos a intentar con rstan
#sc_model <- rstan:: stan(file = stan_fname, data = data, chains = 1, warmup = 90, iter = 100)
sc_model <- cmdstanr::cmdstan_model(stan_fname, pedantic = F, force_recompile = F)


#---------------------------
#Correr el modelo
fit      <- sc_model$sample(data = datos, 
                            refresh=0, iter_warmup = iter_warmup,
                            iter_sampling = nsim,  parallel_chains = pchains,
                            chains = chains, adapt_delta = 0.8)  
save(fit, file = "ModelBetafit.RData")

#Guardamos las simulaciones por si las dudas
modelo_ajustado <- fit$summary()

Hosp <- modelo_ajustado %>% filter(str_detect(variable, "Hosp"))
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

ggplot(Hosp, aes(x = Fecha)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.75, fill = "gray75") +
  geom_line(aes(y = mean, color = "Predichos")) +
  geom_point(aes(y = `Hospitalizados (%)`, color = "Observados"), alpha = 0.5, 
             data = Hosp %>% filter(Fecha <= !!max_fecha)) +
  facet_wrap(~Estado) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_date(breaks = "2 months") +
  theme_classic() +
  scale_color_manual("Modelo", values = c("Observados" = "tomato3", "Predichos" = "gray50"))
  
  
  

