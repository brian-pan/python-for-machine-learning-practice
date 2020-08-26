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
# two adjusted data plots:
plot(sc_adj_data, col = "blue")
lines(adj_data, col="red")
# find Mae:
library(Metrics)
mae(sc_adj_data, adj_data)


##Q4
library(seasonal)
library(ggplot2)

# X11 method:
ua %>% seas(x11="") -> fit_1
autoplot(ua, series="Original Data") + 
  autolayer(seasadj(fit_1), series="X11 Seasonally Adjusted") +
  autolayer(sc_adj_data, series = "StaCan Seasonally Adjusted") +
  xlab("Year") + ylab("data without seasonality") +
  ggtitle("Comparation of X11 method with StaCan's method") +
  scale_colour_manual(values = c("grey", "blue", "red"))
x11_method = seasadj(fit_1)
r_X11 = remainder(fit_1)

# seats method:
ua %>% seas() -> fit_2
autoplot(ua, series="Original Data") + 
  autolayer(seasadj(fit_2), series="Seats Seasonally Adjusted") +
  autolayer(sc_adj_data, series = "StaCan Seasonally Adjusted") +
  xlab("Year") + ylab("data without seasonality") +
  ggtitle("Comparation of Seats method with StaCan's method") +
  scale_colour_manual(values = c("grey", "red", "blue"))
seats_method = seasadj(fit_2)
r_Seats = remainder(fit_2)

# STL method:
log(ua) %>% stl(t.window=13, s.window="periodic", robust=TRUE)-> fit_3
autoplot(ua, series="Original Data") + 
  autolayer(exp(seasadj(fit_3)), series="STL Seasonally Adjusted") +
  autolayer(sc_adj_data, series = "StaCan Seasonally Adjusted") +
  xlab("Year") + ylab("data without seasonality") +
  ggtitle("Comparation of STL method with StaCan's method") +
  scale_colour_manual(values = c("grey", "blue", "red"))
STL_method = exp(seasadj(fit_3))
r_STL = remainder(fit_3)

library(Metrics)
mae(x11_method,sc_adj_data)
mae(seats_method,sc_adj_data)
mae(STL_method,sc_adj_data)

autoplot(ua, series="Original Data") + 
  autolayer(seasadj(fit_1), series="X11 Seasonally Adjusted") +
  autolayer(seasadj(fit_2), series="Seats Seasonally Adjusted") +
  autolayer(exp(seasadj(fit_3)), series="STL Seasonally Adjusted") +
  autolayer(sc_adj_data, series = "StaCan Seasonally Adjusted") +
  xlab("Year") + ylab("data without seasonality") +
  ggtitle("Comparation of all 3 methods with StaCan's method") +
  scale_colour_manual(values = c("light grey", "blue", "grey", "light blue", "red"))

# using unadjusted data:
plot(random_ua)
acf(na.omit(random_ua))
pacf(na.omit(random_ua))
auto.arima(random_ua)

# using unadjusted data with X11 method:
plot(r_X11)
acf(na.omit(r_X11))
pacf(na.omit(r_X11))
auto.arima(r_X11)

# using unadjusted data with STL method:
plot(r_Seats)
acf(r_Seats)
pacf(r_Seats)
auto.arima(r_Seats)

# using StaCan's seasonally adjusted data:
sc_adj_data = get_cansim_vector( "v2057605", start_time = "1980-01-01", end_time = "1999-12-01") %>%
  pull(VALUE) %>% 
  ts( start = c(1980,1), frequency = 12)
library(forecast)
trend_sc = ma(sc_adj_data, order = 12, centre = T)
random_sc = sc_adj_data/trend_sc