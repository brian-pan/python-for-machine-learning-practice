library(astsa)
library(fGarch)
library(magrittr)
library(tidyverse)
library(vars)
library(forecast)
library(Metrics)
library(tseries)

# get unadjusted (raw) series
ua = get_cansim_vector( "v2057814", start_time = "1980-01-01", end_time = "1999-12-01") %>%
  pull(VALUE) %>% 
  ts( start = c(1980,1), frequency = 12)

plot(ua)
acf(ua)
pacf(ua)