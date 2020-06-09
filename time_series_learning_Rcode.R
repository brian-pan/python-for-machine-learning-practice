library(astsa)
library(fGarch)
library(magrittr)
library(tidyverse)
library(vars)
library(forecast)
library(Metrics)
library(tseries)
#1
my_data = read_csv( "ontario_marriages.csv" ) 

toronto_data = my_data %>% 
  filter(YEAR > 1980, CITY == "TORONTO" ) %>% 
  arrange( YEAR, MONTH ) 