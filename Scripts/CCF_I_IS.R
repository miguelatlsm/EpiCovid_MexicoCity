# # install.packages("pacman")
# #install.packages("forecast")
# pacman::p_load(rio,          # File import
#                here,         # File locator
#                tidyverse,    # data management + ggplot2 graphics
#                tsibble,      # handle time series datasets
#                slider,       # for calculating moving averages
#                imputeTS,     # for filling in missing values
#                feasts,       # for time series decomposition and autocorrelation
#                forecast,     # fit sin and cosin terms to data (note: must load after feasts)
#                trending,     # fit and assess models
#                tmaptools,    # for getting geocoordinates (lon/lat) based on place names
#                ecmwfr,       # for interacting with copernicus sateliate CDS API
#                stars,        # for reading in .nc (climate data) files
#                units,        # for defining units of measurement (climate data)
#                yardstick,    # for looking at model accuracy
#                surveillance  # for aberration detection
# )
# library(tidyverse)
# library(vroom)
# library(sf)
# library(lubridate)
# library(plyr)
# library(readxl)
# library(janitor)
# library(dplyr)
# library(forecast)
# library(forecast)
# library(stars)
# library(ggpubr)
# library(lsr)
# library(rstatix)
# library(datarium)
# library(multcomp)
# library(car)
# library(pacman)


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



#### Making serial time object

datos_prueba_ts <- tsibble(datos_prueba, index = fenew)

#Data interpolation

datos_prueba_ts <- datos_prueba_ts %>% 
  mutate(n_int = imputeTS::na_interpolation(n),
         sars_cop_l_int = imputeTS::na_interpolation(sars_cop_l),
         fenew2 = as.Date(fenew,format="%b-%y")
  )

# view(datos_prueba_ts)

#### Data smoothing using 7-day moving averages

datos_prueba_ts_s <- datos_prueba_ts %>% 
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

#View(datos_prueba_ts_s)


###Making CCF graph with interpolated values only
grafccf<-
  datos_prueba_ts %>% 
  CCF( sars_cop_l_int, n_int,
       lag_max = 12,
       type = "correlation")%>%
  autoplot()+
  theme_bw()+
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1),
    name = "cross-correlation coefficient (r)")+
  scale_x_continuous(
    name = "Lag in days")+
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold") 
  ) +
  ggtitle("ZA and TB3 with interpolated values only" )

x11()
grafccf

#Saving graph
ggsave("graf_ccf_i_za_tb3.jpeg", plot = grafccf, dpi = 300, width = 8, height = 6, units = "in")


### Making CCF graph with interpolation and smoothing process
grafccf_s<-
  datos_prueba_ts_s %>% 
  CCF( sars_ma_7d, n_ma_7d,
       lag_max = 12,
       type = "correlation")%>%
  autoplot()+
  theme_bw()+
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1),
    name = "cross-correlation coefficient (r)")+
  scale_x_continuous(
    name = "Lag in days")+
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold") # Tamaño del título de cada faceta
  ) +
  ggtitle("ZC and TB1 with interpolation and smoothing process" )

x11()
grafccf_s

#Saving graph
#ggsave("graf_ccf_zc_tb1.jpeg", plot = grafccf_s, dpi = 300, width = 8, height = 6, units = "in")


### Calculation of CCF values with interpolated values only
ccf <- datos_prueba_ts_s %>% 
  CCF( sars_cop_l_int, n_int,
       lag_max = 15,
       type = "correlation")%>%
  arrange(-ccf) %>% 
  slice_head(n = 15)

###Calculation of CCF values with interpolation and smoothing process

ccf_s <- datos_prueba_ts_s %>% 
  CCF( sars_ma_7d, n_ma_7d,
       lag_max = 15,
       type = "correlation")%>%
  arrange(-ccf) %>% 
  slice_head(n = 10)

### Calculation of CCF values with interpolated values only
ccf
###Calculation of CCF values with interpolation and smoothing process
ccf_s



