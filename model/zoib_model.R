rm(list = ls())
library(readr)
library(tidyverse)
library(lubridate)
library(zoib)
library(dplyr)
library(ggplot2)
library(zoo)

autoregresive_order <- 7

#-------------------- Lectura y seleccion de Datos -------------------
#Lectura de la base de datos 
githubURL <- ("https://github.com/RodrigoZepeda/CapacidadHospitalariaMX/blob/master/processed/HospitalizacionesMX.rds?raw=true")
download.file(githubURL,"HospitaizacionesMX.rds")
hospitalizaciones <- readRDS("HospitaizacionesMX.rds")

#Formateo de la base de datos
hospitalizaciones <- hospitalizaciones %>%
  mutate(Fecha  = as.Date(Fecha))      %>%                          
  mutate(Dias   = as.numeric(Fecha - ymd("2020/04/01"))) %>%
  select(Dias, Estado, CLUES, `Hospitalizados (%)`)     %>%
  arrange(Dias) %>%
  distinct(Dias, CLUES, .keep_all = T) %>%
  mutate(hosp = replace_na(`Hospitalizados (%)`, 0)) %>%
  mutate(hosp = hosp/100) %>% 
  mutate(hosp = case_when(hosp < 0 ~ 0, hosp > 1 ~ 1, TRUE ~ hosp)) 

#Agregamos las variables del pasado
for (i in 1:autoregresive_order){
  hospitalizaciones <- hospitalizaciones %>%
    arrange(Dias) %>%
    group_by(CLUES) %>%
    mutate(!!sym(paste0("hosp_", i)) := lag(hosp, i))
}
hospitalizaciones <- hospitalizaciones %>% 
  filter(Dias > autoregresive_order) %>%
  drop_na()




#Great explanation of zoib
#https://stats.stackexchange.com/a/380158/81181
formula <- paste0("hosp ~ ", paste0("hosp_", 1:autoregresive_order, collapse = " + "))
formula <- paste0(formula, " | as.factor(Estado) + as.factor(CLUES) | hosp_1 | hosp_1 ")
#formula <- paste0(formula, " | Dias | hosp_1 | hosp_1 ") #TODO Vale cambiar para usar esta

modelo <- zoib(as.formula(formula),
               random = 0, 
               n.chain = 2,
               n.iter = 500,
               EUID = as.numeric(as.factor(hospitalizaciones$CLUES)),
               data = as.data.frame(hospitalizaciones), 
               zero.inflation = T, one.inflation = T,
               joint = F)

predict_model <- pred.zoib(modelo, hospitalizaciones)

predats <- predict_model$summary %>% 
  as.data.frame() %>% bind_cols(hospitalizaciones) 



