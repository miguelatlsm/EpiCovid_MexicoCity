# # install.packages("pacman")
# #install.packages("forecast")
pacman::p_load(rio,          # File import
               here,         # File locator
               tidyverse,    # data management + ggplot2 graphics
               tsibble,      # handle time series datasets
               slider,       # for calculating moving averages
               imputeTS,     # for filling in missing values
               feasts,       # for time series decomposition and autocorrelation
               forecast,     # fit sin and cosin terms to data (note: must load after feasts)
               trending,     # fit and assess models
               tmaptools,    # for getting geocoordinates (lon/lat) based on place names
               ecmwfr,       # for interacting with copernicus sateliate CDS API
               stars,        # for reading in .nc (climate data) files
               units,        # for defining units of measurement (climate data)
               yardstick,    # for looking at model accuracy
               surveillance  # for aberration detection
)
library(tidyverse)
library(vroom)
library(sf)
library(lubridate)
library(plyr)
library(readxl)
library(janitor)
library(dplyr)
library(forecast)
library(forecast)
library(stars)
library(ggpubr)
library(lsr)
library(rstatix)
library(datarium)
library(multcomp)
library(car)
library(ggplot2)


################################################################################

#Database containing all the results of the SARS-CoV-2 N1 gene quantification 
#and the reported cases for each analyzed zone (A, B, C)

datos_prueba <- read_excel("C:/R/WBE_Covid_CMX/Data/ZA_B_C_TB1_2_3_sum_fenew_21_22.xlsx")


# Clining data

datos_prueba <- 
  datos_prueba %>% 
  janitor::clean_names()%>% 
  mutate(fenew = fenew %>% ymd(), 
         fecinisi = fecinisi %>% ymd(),
         fesars = fesars %>% ymd(),
         sars_cop_l=sars_cop_l)

#### Tendenci mesures ####


datos_tendencia <- 
  datos_prueba %>% 
  janitor::clean_names()%>% 
  mutate(fenew = fenew %>% ymd(), 
         fecinisi = fecinisi %>% ymd(),
         fesars = fesars %>% ymd(),
         sars_cop_l=sars_cop_l,
         muestreo=as.factor(muestreo))%>%
  filter(!is.na(sars_cop_l))

#View(datos_tendencia)

tapply(datos_tendencia$sars_cop_l,datos_tendencia$muestreo, mean)
tapply(datos_tendencia$sars_cop_l,datos_tendencia$muestreo, max)
tapply(datos_tendencia$sars_cop_l,datos_tendencia$muestreo, min)
tapply(datos_tendencia$sars_cop_l,datos_tendencia$muestreo, var)
tapply(datos_tendencia$sars_cop_l,datos_tendencia$muestreo, sd)



### Initial graph gen N1 gen/L SARS CoV 2 vs date 


Sys.setlocale("LC_TIME", "English")
graf_sars_cop <- 
  datos_prueba %>% 
  janitor::clean_names() %>% 
  dplyr::select(fenew, sars_cop_l, de_sars_cop_l, muestreo) %>%
  filter(!is.na(sars_cop_l)) %>%
  mutate(fenew = as.Date(fenew, format="%d/%b/%y"), 
         muestra="Copilco") %>%
  ggplot(aes(fenew, sars_cop_l, color = as.factor(muestreo)))  +
  theme_bw() +
  
  geom_line(size = 0.5) +  
  geom_errorbar(aes(ymin = sars_cop_l - de_sars_cop_l,
                    ymax = sars_cop_l + de_sars_cop_l), 
                color = "black", width = 0.2) +
  geom_point(size = 2) +
  

  scale_y_continuous(
    breaks = seq(0, 35000, 5000),
    limits = c(0, 35000),  
    name = "N1 gene copies/L SARS-CoV-2"
  ) +
  
  scale_x_date(
    date_breaks = "1 month", date_labels = "%b-%y",
    name = "Sampling Date"
  ) +

  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),  
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 12, face = "bold"),
    legend.background = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = 16, face = "bold", colour = "brown"),
    plot.subtitle = element_text(size = 14, face = "bold")
  ) +
  
  guides(colour = guide_legend(title = "Sampling season"))

#Saving graph

#ggsave("graf_sars_cop.jpeg", plot = graf_sars_cop, dpi = 300, width = 8, height = 6, units = "in")


X11()
graf_sars_cop



###Boxplot to compare data with Wilcoxon test


box_plot_datos_prueba <-
  t_datos_prueba %>%
  ggplot(aes(x = muestreo, y = sars_cop_l, color = muestreo)) +
  theme_bw() +
  geom_boxplot() +
  stat_summary(fun.y = mean, 
               geom = "point", 
               shape = 17, 
               size = 3, 
               color = "red") +
  geom_jitter(shape = 16, 
              position = position_jitter(0.2)) +
  scale_y_continuous(breaks = seq(0, 35000, 5000)) +
  scale_colour_discrete(
    name = "Sampling Date",
    labels = c("21/apr/2021 to 8/sep/2021", "30/nov/2021 to 24/jan/2022")) +
  labs(
    x = "Sampling season",
    y = "N1 gene copies/L SARS-CoV-2"
  ) +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = rel(1.2), lineheight = .9, face = "bold", colour = "brown"),
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 11, color = "black"), 
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 11, face = "bold")   
  ) +
  stat_compare_means(method = "wilcox.test", 
                     label.x = 1.5, 
                     label.y = 20000, 
                     size = 5, 
                     fontface = "bold")  
X11()
box_plot_datos_prueba

#Saving Graph
#ggsave("box_plot_datos_prueba.jpeg", plot = box_plot_datos_prueba, dpi = 300, width = 8, height = 6, units = "in")




#### For graphs of reported cases by each zone (A, B or C) ###



colors_i <- c("Zone A (cases level 1)" = "#750275", 
              "Zone B (cases level 1+2)" = "#3198FF",
              "Zone C (cases level 1+2+3)" = "#FE6100")

Sys.setlocale("LC_TIME", "English")

graf_n_123_12_1_cop <- 
  datos_prueba_123_12_1 %>% 
  janitor::clean_names() %>%
  dplyr::select(fenew, n_123, n_12, n_1) %>%
  filter(!is.na(n_123)) %>%
  mutate(fenew = as.Date(fenew, format="%d/%b/%y"), 
         muestra="Copilco") %>%
  ggplot(aes(x = fenew)) +
  theme_bw() +
  

  geom_line(aes(y = n_123, color = "Zone C (cases level 1+2+3)"), linewidth = 0.7) +
  geom_line(aes(y = n_12, color = "Zone B (cases level 1+2)"), linewidth = 0.7) +
  geom_line(aes(y = n_1, color = "Zone A (cases level 1)"), linewidth = 0.7) +
  

  geom_point(aes(y = n_123, color = "Zone C (cases level 1+2+3)"), size = 1) +
  geom_point(aes(y = n_12, color = "Zone B (cases level 1+2)"), size = 1) +
  geom_point(aes(y = n_1, color = "Zone A (cases level 1)"), size = 1) +

  scale_y_continuous(
    breaks = seq(0, 50, 10),
    name = "Cases reported by zip code"
  ) +
  

  scale_x_date(
    date_breaks = "1 month", date_labels = "%b-%y",
    name = "Sampling date"
  ) +

  scale_color_manual(values = colors_i) +

  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.text = element_text(size = 14, face = "bold"),  # Texto de la leyenda en tamaño 16 y negrita

    axis.title = element_text(size = 12, face = "bold"),   # Tamaño 11 y negrita en los títulos de los ejes

    axis.text = element_text(size = 11),                   # Tamaño 11 en los textos de los ejes

    plot.title = element_text(
      size = rel(1.2), lineheight = .9,
      face = "bold", colour = "brown"
    )
  ) +

  theme(axis.text.x = element_text(angle = 90)) +
  guides(colour = guide_legend(title = "")) +
  
  ggtitle("", subtitle = "")


X11()
graf_n_123_12_1_cop

#Saving graph
#ggsave("graf_n_123_12_1_cop.jpeg", plot = graf_n_123_12_1_cop, dpi = 300, width = 16, height = 6, units = "in")

