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

random_ua = ua / (trend_ua*seasonal_ua)

# combine all terms and rename them
observed = ua
trend = trend_ua
seasonal = rep(seasonal_ua,20)
random = random_ua
# plot the renamed terms
plot(cbind(observed, trend, seasonal, random),
     main = "multiplicative time series decomposition",
     xlab = "Years")


## Q3
# get StatCan's adjusted data:
sc_adj_data = get_cansim_vector( "v2057605", start_time = "1980-01-01", end_time = "1999-12-01") %>%
  pull(VALUE) %>% 
  ts( start = c(1980,1), frequency = 12)
# make unadjusted data (ua) to be seaonally adjusted:
mat_ua_2 = t(matrix(data = detrend_ua, nrow = 12))
seasonal_ua_2 = colMeans(mat_ua_2, na.rm = TRUE)
adj_data = ua / seasonal_ua_2