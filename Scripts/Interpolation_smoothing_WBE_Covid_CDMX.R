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


################################################################################

#You need to select the database that corresponds to the area and period you want to analyze.


datos_prueba <- read_excel("C:/R/WBE_Covid_CMX/Data/ZA_sum_fenew_21_22_TB1.xlsx") #Zone A 21 april to 22 feb 
# datos_prueba <- read_excel("C:/R/WBE_Covid_CMX/Data/ZA_sum_fenew_21_22_TB2.xlsx") #Zone A 21 april to 21 sep 
# datos_prueba <- read_excel("C:/R/WBE_Covid_CMX/Data/ZA_sum_fenew_21_22_TB3.xlsx") #Zone A 21 nov to 22 feb 
# 
# datos_prueba <- read_excel("C:/R/WBE_Covid_CMX/Data/ZB_sum_fenew_21_22_TB1.xlsx") #Zone B 21 april to 22 feb
# datos_prueba <- read_excel("C:/R/WBE_Covid_CMX/Data/ZB_sum_fenew_21_22_TB2.xlsx") #Zone B 21 april to 21 sep 
# datos_prueba <- read_excel("C:/R/WBE_Covid_CMX/Data/ZB_sum_fenew_21_22_TB3.xlsx") #Zone B 21 nov to 22 feb 
# 
# datos_prueba <- read_excel("C:/R/WBE_Covid_CMX/Data/ZC_sum_fenew_21_22_TB1.xlsx") #Zone C 21 april to 22 feb
# datos_prueba <- read_excel("C:/R/WBE_Covid_CMX/Data/ZC_sum_fenew_21_22_TB2.xlsx") #Zone C 21 april to 21 sep
# datos_prueba <- read_excel("C:/R/WBE_Covid_CMX/Data/ZC_sum_fenew_21_22_TB3.xlsx") #Zone C 21 nov to 22 feb 



# Clean data

datos_prueba <- 
  datos_prueba %>% 
  janitor::clean_names()%>% 
  mutate(fenew = fenew %>% ymd(), 
         fecinisi = fecinisi %>% ymd(),
         fesars = fesars %>% ymd(),
         sars_cop_l=sars_cop_l)

#View(datos_prueba_ts)

#### Making serial time object

datos_prueba_ts <- tsibble(datos_prueba, index = fenew)


#Data interpolation


datos_prueba_ts <- datos_prueba_ts %>% 
  mutate(n_int = imputeTS::na_interpolation(n),
         sars_cop_l_int = imputeTS::na_interpolation(sars_cop_l),
         fenew2 = as.Date(fenew,format="%b-%y")
  )


#### Interpolated graphs of SARS-CoV 2 samples from Water

graf_di_agua <- ggplot_na_imputations(
  datos_prueba_ts$sars_cop_l, datos_prueba_ts$sars_cop_l_int,
  x_axis_labels = datos_prueba_ts$fenew2,
  color_points = "steelblue",
  color_imputations = "indianred",
  size_imputations = 2, size_points = 2,
  size_lines = 0.5,
  label_known = "Known data",
  label_imputations = "Imputed data"
) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(
    breaks = seq(0, 30000, 5000),
    name = "N1 gene copies/L SARS-CoV-2"
  ) +
  scale_x_date(
    date_breaks = "1 month", date_labels = "%b-%y",
    name = "Sampling season"
  ) +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(
      size = rel(1.2), lineheight = .9,
      face = "bold", colour = "brown"
    ),
    axis.title = element_text(size = 12, face = "bold"), 

    axis.text = element_text(size = 11), 
    legend.text = element_text(size = 14, face = "bold")  
  ) +
  
  ggtitle("", subtitle = "Data interpolation, gene copies/L, Season 1 and 2")+
  theme(plot.subtitle = element_text(size = 14, face = "bold"))


x11()
graf_di_agua

#Saving graph
#ggsave("graf_di_agua.jpeg", plot = graf_di_agua, dpi = 300, width = 16, height = 6, units = "in")



#### Interpolated graphs of SARS-CoV 2 samples from reported data

graf_di_pruebas <- ggplot_na_imputations(
  datos_prueba_ts$n, datos_prueba_ts$n_int,
  x_axis_labels = datos_prueba_ts$fenew2,
  color_points = "steelblue",
  color_imputations = "indianred",
  size_imputations = 2, size_points = 2,
  size_lines = 0.5,
  label_known = "Known data",
  label_imputations = "Imputed data"
) +
  theme_bw() +

  theme(axis.text.x = element_text(angle = 90)) +

  scale_y_continuous(
    breaks = seq(0, 50, 10),
    name = "Cases reported by zip code"
  ) +

  scale_x_date(
    date_breaks = "1 month", date_labels = "%b-%y",
    name = "Sampling season"
  ) +

  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", colour = "black"),

    plot.title = element_text(
      size = rel(1.2), lineheight = .9,
      face = "bold", colour = "brown"
    ),

    axis.title = element_text(size = 12, face = "bold"), 

    axis.text = element_text(size = 11),  

    legend.text = element_text(size = 14, face = "bold")  
  ) +
  
  ggtitle("", subtitle = "Data imputation, Zone C (Level 1+2+3)")+
  theme(plot.subtitle = element_text(size = 14, face = "bold"))


x11()
graf_di_pruebas

#Saving graph
#ggsave("graf_di_pruebas.jpeg", plot = graf_di_pruebas, dpi = 300, width = 16, height = 6, units = "in")


#### Data smoothing using 7-day moving averages


datos_prueba_ts <- datos_prueba_ts %>% 
  ## create the ma_4w variable 
  ## slide over each row of the case variable
  mutate(n_ma_7d = slider::slide_dbl(n_int, 
                                     ## for each row calculate the mean
                                     ~ mean(.x, na.rm = TRUE),
                                     ## use the four previous weeks
                                     .before = 3, .after = 3, .complete = TRUE),
         sars_ma_7d = slider::slide_dbl(sars_cop_l_int, 
                                        ## for each row calculate the mean
                                        ~ mean(.x, na.rm = TRUE),
                                        ## use the four previous weeks
                                      .before = 3, .after = 3, .complete = TRUE)
  )



#View(datos_prueba_ts)

#Escribir la nueva base de datos con la base de datos inputada y ya suavisada 
# el nombre final quedo como  ¿?_It_Smo_cop_tot con ¿? como 1, 12, 123

#write.csv(datos_prueba_ts, "C:/R/Master/Niveles/Inter_Smooth/¿?_It_Smo_cop_tot.csv")


##### Visualizar datos ya interpolados y suavisados #####

#### para todo


# graf_suv_todo <- ggplot(datos_prueba_ts, aes(x = fenew)) +
#    theme_light()+
#    geom_line(aes(y = sars_ma_7d), color= "black")+
#    geom_line(aes(y = n_ma_7d*600), color= "green") +
#    scale_y_continuous(name = "cop/L gen N1 SARS-CoV-2, MA 7 días",
#      sec.axis = sec_axis(~./400, name = "Casos reportados, MA 7 días"))+
#    scale_x_date(
#      date_breaks = "1 month", date_labels = "%b-%y",
#      name = "Fecha de muestreo")+
#    theme(axis.text.x = element_text(angle = 90))+
#    ggtitle("",
#             subtitle = "cop/L gen N1 SARS-CoV-2 Bloque 1,2 / casos reportados Bloque C (nivel 1+2+3)" )
#  X11()
#  graf_suv_todo

# 
# #### para enero
# 
 colors_ene <- c("Season 2 (gene N1 cop/L)" = "red", 
                 "Zone C (cases level 1+2+3)" = "black")
 
 graf_suv_ene <- ggplot(datos_prueba_ts, aes(x = fenew)) +
   theme_bw()+
   geom_line(aes(y = sars_ma_7d, color = "Season 2 (gene N1 cop/L)"),linewidth = 0.7)+
   geom_line(aes(y = n_ma_7d*600, color = "Zone C (cases level 1+2+3)"),linewidth = 0.7) +
   
   scale_y_continuous(
     name = "Gene N1 cop/L SARS-CoV-2, SMA 7 days",
     sec.axis = sec_axis(~./400, name = "Cases reported, SMA 7 days")
   )+
   scale_x_date(
     date_breaks = "1 month", 
     date_labels = "%b-%y",
     name = "Sampling Date"
   )+
   # Configuración de colores manuales
   scale_color_manual(values = colors_ene) +
   theme(axis.text = element_text(angle = 90))+
   theme(
     legend.position = "bottom",
     legend.background = element_rect(fill = "white", colour = "black"),
     
     # Tamaño y negrita de los títulos de los ejes
     axis.title = element_text(size = 12, face = "bold"),  
     
     # Tamaño del texto en los ejes
     axis.text = element_text(size = 11),  
     
     # Tamaño y negrita del texto en la leyenda
     legend.text = element_text(size = 14, face = "bold"),  
     
     # Tamaño y negrita del subtítulo
     plot.subtitle = element_text(size = 14, face = "bold")
   ) +
   
   # Remover título de la leyenda
   guides(color = guide_legend(title = NULL)) +
   
   # Título de la gráfica
   ggtitle("", subtitle = "Season 2 (Gene N1 cop/L) / Cases reported Zone C")
 
 
 x11()
 graf_suv_ene
 
 ggsave("graf_suv_ene.jpeg", plot = graf_suv_ene, dpi = 300, width = 8, height = 6, units = "in")

#### para septiembre
  
colors_sep <- c("Season 1 (gene N1 cop/L)" = "blue", 
              "Zone C (cases level 1+2+3)" = "black")

graf_suv_sep <- ggplot(datos_prueba_ts, aes(x = fenew)) +
  theme_bw()+
  geom_line(aes(y = sars_ma_7d, color = "Season 1 (gene N1 cop/L)"),linewidth = 0.7)+
  geom_line(aes(y = n_ma_7d*2000, color = "Zone C (cases level 1+2+3)"),linewidth = 0.7) +
  
  scale_y_continuous(
    name = "Gene N1 cop/L SARS-CoV-2, SMA 7 days",
    sec.axis = sec_axis(~./400, name = "Cases reported, SMA 7 days")
    )+
  scale_x_date(
    date_breaks = "1 month", 
    date_labels = "%b-%y",
    name = "Sampling Date"
    )+
  # Configuración de colores manuales
  scale_color_manual(values = colors_sep) +
  theme(axis.text = element_text(angle = 90))+
theme(
  legend.position = "bottom",
  legend.background = element_rect(fill = "white", colour = "black"),
  
  # Tamaño y negrita de los títulos de los ejes
  axis.title = element_text(size = 12, face = "bold"),  
  
  # Tamaño del texto en los ejes
  axis.text = element_text(size = 11),  
  
  # Tamaño y negrita del texto en la leyenda
  legend.text = element_text(size = 14, face = "bold"),  
  
  # Tamaño y negrita del subtítulo
  plot.subtitle = element_text(size = 14, face = "bold")
) +
  
  # Remover título de la leyenda
  guides(color = guide_legend(title = NULL)) +
  
  # Título de la gráfica
  ggtitle("", subtitle = "Season 1 (Gene N1 cop/L) / Cases reported Zone C")



x11()
graf_suv_sep

##Saving graph
#ggsave("graf_suv_sep.jpeg", plot = graf_suv_sep, dpi = 300, width = 8, height = 6, units = "in")



