#install.packages("pacman")
#install.packages("forecast")
library(tidyverse)
library(vroom)
library(sf)
library(lubridate)
library(plyr)
library(readxl)
library(janitor)
library(dplyr)
library(forecast)


################################################################################

#Database containing all the results of the SARS-CoV-2 N1 gene quantification 
#and the reported cases for each analyzed zone (A, B, C), interpolated and smooted. 


datos_gcor_tot <- read_csv("C:/R/WBE_Covid_CMX/Data/ZA_B_C_TB1_2_3_sum_fenew_21_22_Smooth.csv")



datos_gcor_tot <-
  datos_gcor_tot %>%
  janitor::clean_names() %>%
  dplyr::select(no, proximidad,fenew,sars_ma_7d, n_ma_7d) %>%
  mutate(fecha = fenew)%>%
  mutate(fecha = as.Date(fecha,format="%d/%m/%Y"))%>%
  mutate(fe = fenew)%>%
  mutate(fe = as.Date(fecha,format="%d/%m/%Y"))



####grouping data  ####

#### N1 SARS-CoV-2 RT-qPCR gene counts 
datos_sars_tot <-
  datos_gcor_tot %>%
  janitor::clean_names() %>%
  dplyr::select(no, proximidad,sars_ma_7d,fe)%>%
  filter(!is.na(sars_ma_7d),!is.na(proximidad))

#view(datos_sars_tot)

#### Case reported by zone
datos_n_tot <-
  datos_gcor_tot %>%
  janitor::clean_names() %>%
  dplyr::select(no, proximidad,n_ma_7d,fecha)%>%
  filter(!is.na(n_ma_7d))

#view(datos_n_tot)


###### Making Lag timing,  ##### 

## TB1 21 april to 22 feb 

exp_tot_tot<- seq(from = as.Date("2021-04-08"), 
              to = as.Date ("2022-03-01"), 
              by=1)

mi_expansion_tot_tot <- expand.grid(fecha = exp_tot_tot,lag = -8:-6)

#view(mi_expansion_ene_tot)

## TB2 21 april to 21 sep 

exp_sep_tot<- seq(from = as.Date("2021-04-08"), 
              to = as.Date ("2021-09-15"), 
              by=1)

mi_expansion_sep_tot <- expand.grid(fecha = exp_sep_tot,lag = -9:-7)


## TB3 21 nov to 22 feb 


exp_ene_tot<- seq(from = as.Date("2021-11-22"), 
              to = as.Date ("2022-03-01"), 
              by=1)

mi_expansion_ene_tot <- expand.grid(fecha = exp_ene_tot,lag = -7:-5)




#Sample dates are created here with the quantifications of
#SARS-CoV-2 RT-qPCR N1 gene counts and the proposed lag dates

## TB1 21 april to 22 feb 

cuant_iter_smo_tot_tot <- 
  datos_n_tot %>%
  left_join(y = mi_expansion_tot_tot) %>% 
  mutate(fecha.true = fecha,
         fecha.lag = fecha+lag,
         proximidad =proximidad,
         bloque = "tot") %>%
  filter(!is.na(fecha.lag))


#View(cuant_iter_smo_tot_tot)


## TB2 21 april to 21 sep 


cuant_iter_smo_sep_tot <- 
  datos_n_tot %>%
  left_join(y = mi_expansion_sep_tot) %>% 
  mutate(fecha.true = fecha,
         fecha.lag = fecha+lag,
         proximidad =proximidad,
         bloque = "sep") %>%
  filter(!is.na(fecha.lag))


#View(cuant_iter_smo_sep)

## TB3 21 nov to 22 feb 
cuant_iter_smo_ene_tot <- 
  datos_n_tot %>%
  left_join(y = mi_expansion_ene_tot) %>% 
  mutate(fecha.true = fecha,
         fecha.lag = fecha+lag,
         proximidad =proximidad,
         bloque = "ene") %>%
  filter(!is.na(fecha.lag))


#View(cuant_iter_smo_ene)



####Making data lags


## TB1 21 april to 22 feb 

final_prox123_tot_tot<- left_join(cuant_iter_smo_tot_tot, datos_sars_tot, 
                              by =c('fecha.lag'='fe'))%>%
  filter(!is.na(sars_ma_7d))


#View(final_prox123_tot_tot)

## TB2 21 april to 21 sep 

final_prox123_sep_tot<- left_join(cuant_iter_smo_sep_tot, datos_sars_tot, 
                              by =c('fecha.lag'='fe'))%>%
  filter(!is.na(sars_ma_7d))


#View(final_prox123_sep)

## TB3 21 nov to 22 feb

final_prox123_ene_tot<- left_join(cuant_iter_smo_ene_tot, datos_sars_tot, 
                              by =c('fecha.lag'='fe'))%>%
  filter(!is.na(sars_ma_7d))


#View(final_prox123_ene)


##### Making correlation graphs #####

## TB1 21 april to 22 feb 
grafinal_tot_tot <-
  final_prox123_tot_tot %>%
  arrange(fecha, lag) %>% 
  ggplot() +
  theme_bw() +
  aes(sars_ma_7d, n_ma_7d, color = as.factor(proximidad.x)) +
  geom_point() +
  scale_color_manual(name = "",
                     labels = c("Zone A",
                                "Zone B",
                                "Zone C"),
                     values = c("#750275", "#3198FF", "#FE6100")) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(lag)) +
  labs(x = "Gene N1 cop/L SARS-CoV-2, SMA 7 days",
       y = "Reported cases, SMA 7 days") +
  
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold") # Tamaño del título de cada faceta
  ) +
  ggtitle("", subtitle = "TB1 (2021-apr to 2022-feb)")

ggsave("grafinal_tot_tot.jpeg", plot = grafinal_tot_tot, dpi = 300, width = 12, height = 5, units = "in")


X11()
grafinal_tot_tot

## TB2 21 april to 21 sep
grafinal_sep_tot <-
  final_prox123_sep_tot %>%
  arrange(fecha, lag) %>% 
  ggplot() +
  theme_bw() +
  aes(sars_ma_7d, n_ma_7d, color = as.factor(proximidad.x)) +
  geom_point() +
  scale_color_manual(name = "",
                     labels = c("Zone A",
                                "Zone B",
                                "Zone C"),
                     values = c("#750275", "#3198FF", "#FE6100")) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(lag)) +
  labs(x = "Gene N1 cop/L SARS-CoV-2, SMA 7 days",
       y = "Reported cases, SMA 7 days") +
  
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold") # Tamaño del título de cada faceta
  ) +
  ggtitle("", subtitle = "TB2 (2021-apr to 2021-sep)")

ggsave("grafinal_sep_tot.jpeg", plot = grafinal_sep_tot, dpi = 300, width = 12, height = 5, units = "in")


X11()
grafinal_sep_tot

## TB3 21 nov to 22 feb
grafinal_ene_tot <-
  final_prox123_ene_tot %>%
  arrange(fecha, lag) %>% 
  ggplot() +
  theme_bw() +
  aes(sars_ma_7d, n_ma_7d, color = as.factor(proximidad.x)) +
  geom_point() +
  scale_color_manual(name = "",
                     labels = c("Zone A",
                                "Zone B",
                                "Zone C"),
                     values = c("#750275", "#3198FF", "#FE6100")) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(lag)) +
  labs(x = "Gene N1 cop/L SARS-CoV-2, SMA 7 days",
       y = "Reported cases, SMA 7 days") +
  
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold") # Tamaño del título de cada faceta
  ) +
  ggtitle("", subtitle = "TB3 (2021-nov to 2022-feb)")

ggsave("grafinal_ene_tot.jpeg", plot = grafinal_ene_tot, dpi = 300, width = 12, height = 5, units = "in")


X11()
grafinal_ene_tot

