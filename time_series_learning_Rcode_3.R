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

# apply the 12-point MA to the original series:
library(forecast)
trend_ua = ma(ua, order = 12, centre = T)
plot(ua)
lines(trend_ua, col="red")

# get detrended data:
detrend_ua = ua/trend_ua

mat_ua = t(matrix(data = detrend_ua, nrow = 12))
seasonal_ua = colMeans(mat_ua, na.rm = TRUE)
mean(seasonal_ua)
plot(as.ts(rep(seasonal_ua,20)), ylab="")
abline(h=mean(seasonal_ua), col="red")