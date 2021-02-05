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
hospitalizaciones <- readRDS("/Users/rod/Dropbox/HospitalizacionesCOVIDMX/processed/HospitalizacionesMX.rds")

#Formateo de la base de datos
hospitalizaciones <- hospitalizaciones %>%
  mutate(Fecha  = as.Date(Fecha))      %>%                          
  mutate(Estado = as.factor(Estado))   %>%
  mutate(CLUES  = as.factor(CLUES))    %>%
  mutate(Dias   = as.numeric(Fecha - ymd("2020/04/01"))) %>%
  select(Dias, Estado, CLUES, `Hospitalizados (%)`)     %>%
  arrange(Dias) %>%
  distinct(Dias, CLUES, .keep_all = T) %>%
  mutate(hosp = replace_na(`Hospitalizados (%)`, 0)) %>%
  mutate(hosp = hosp/100) %>% 
  mutate(hosp = case_when(hosp < 0 ~ 0, hosp > 1 ~ 1, TRUE ~ hosp))  %>%
  filter(as.numeric(Estado) < 5) %>%
  filter(Dias < max(Dias) - 50)

#Agregamos las variables del pasado
for (i in 1:autoregresive_order){
  hospitalizaciones <- hospitalizaciones %>%
    mutate(!!sym(paste0("hosp_", i)) := lag(hosp, i))
}
hospitalizaciones <- hospitalizaciones %>% 
  mutate(Dias >= min(Dias) + autoregresive_order)

#Great explanation of zoib
#https://stats.stackexchange.com/a/380158/81181
formula <- paste0("hosp ~ ", paste0("hosp_", 1:autoregresive_order, collapse = " + "))
formula <- paste0(formula, " | Estado + CLUES | 1 | 1")

modelo <- zoib(as.formula(formula),
               random = 0, 
               EUID = as.numeric(hospitalizaciones$CLUES), 
               n.chain = 4,
               n.iter = 100,
               data = hospitalizaciones, zero.inflation = TRUE, one.inflation = TRUE,
               joint = F)




