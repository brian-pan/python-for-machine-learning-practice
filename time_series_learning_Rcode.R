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

marriages = ts( toronto_data$COUNT, start=c(1980,1), frequency=12)
plot(marriages, type='o', pch=20)
acf(marriages)
# Use the ```window``` function to extract a subset of the series from Jan 1990 till Dec 2000, and plot it.
window( marriages, start = c(1990,1), end = c(2000, 12) )

N=500
W = ts(rnorm(N))
acf(W, main = "")

R = cumsum(W)
acf(R)

X = W + c( 0, W[-N])
acf(X)

#2
library(astsa)
plot( jj, type="o", pch =20)
Tm = as.vector( time(jj) )
Qr = as.factor( cycle(jj) )

out = lm( log(jj) ~ Tm + Qr -1 )
summary(out)

exp( .25*out$coefficients["Tm"] + out$coefficients["Qr4"] - out$coeff["Qr3"] ) - 1
plot(jj, main ="Original"); lines( Tm, exp( fitted(out) ), col = 2 )
plot( log(jj), main ="Transformed"); lines( Tm, fitted(out), col = 2 )

plot(out$residuals, type = 'o', pch = 20)
acf(out$residuals)

library(tidyverse)
read_csv( url( "https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.csv" ), skip = 1, na = "***"  ) %>% 
  pivot_longer( "Jan":"Dec", names_to = "Month", values_to = "Temp") %>% 
  drop_na(Temp) %>% 
  # uncomment below to check data 
  # View() 
  pull( "Temp" ) %>% 
  ts( start = c(1880,1), frequency = 12 ) -> temp
plot( temp ) 

plot(temp)
stats::filter(temp, c(.5, rep(1,99), .5)/100) %>% 
  lines( x = time(temp), y = ., type = "l", lwd = 2, col = 2 )

# Fit cubic trend
Tm = as.vector(time(temp)); Tm2 = Tm^2; Tm3 = Tm^3
out = lm( temp ~ Tm + Tm2 + Tm3  )

# Create predictions and find first time for which it is > 2
new_Tm = tail(Tm,1) + cumsum( rep(1/12, 500 ) )
pred = predict(out, newdata = data.frame(Tm = new_Tm, Tm2 = new_Tm^2, Tm3 = new_Tm^3) )
new_Tm[ which(pred>2) ] %>% head(1)

# Plot
plot(temp, ylim = c(-1, 2.3), xlim = range( c(Tm,new_Tm)))
abline(h=2, lty = 2, col = 2)
lines( Tm, out$fitted, type = "l", col = 2)
lines( new_Tm, pred, type = "l", col = 2)


read_csv( "ontario_marriages.csv" ) %>%  
  filter(YEAR >= 1980 ) %>% 
  group_by(YEAR, MONTH) %>% 
  summarise(COUNT = sum(COUNT)) %>% 
  arrange(YEAR, MONTH) %>%  
  pull(COUNT) %>% 
  ts( start=c(1980,1), frequency=12) -> marriages

plot(marriages, type='o', pch=20)

dcmp = decompose(marriages)
plot(dcmp)

V_X = var(marriages)
V_T = var(dcmp$trend, na.rm = T)
V_S = var(dcmp$seasonal)
V_R = var(dcmp$random, na.rm = T)

V_T/V_X 
V_S/V_X 
V_R/V_X 

#F_T 
max(0, 1 - V_R / var( dcmp$trend + dcmp$random, na.rm = T ) )

#F_S
max(0, 1-V_R / var( dcmp$seasonal + dcmp$random, na.rm = T ))


#3
n=100
x = stats::filter( rnorm(n), c(0,-.9), method = "recursive")
par(mfrow=c(1,2))
plot(x); acf(x, main = "")

# ma
x_ma = stats::filter( x, rep(1/4,4), side = 1 )
plot(x_ma); acf(x_ma, main = "", na.action = na.pass )

# ar simulate
n = 1000
par(mfrow=c(1,2))
x = filter( rnorm(n), c(.75, -.5), method = "recursive")
plot(x); acf(x, main = "")
# ar simulate (method 2)
x = arima.sim( model = list( ar = c(.75, -.5) ), n )
par(mfrow=c(1,2))
plot(x); acf(x, main = "")

psi = ARMAtoMA( ar = c(.75,-.5), lag.max = 50)
plot( psi, type = "h" )
# Repeat the previous steps
x = filter( rnorm(n), c(.75, .5), method = "recursive")
par(mfrow=c(1,2))
plot(x); acf(x, main = "")

try(
  arima.sim( model = list( ar = c(.75, .5) ), n )
)
psi = ARMAtoMA( ar = c(.75,.5), lag.max = 50)
plot( psi, type = "h" )

x = filter( rnorm(n), c(1, .75, .5), side = 1, method = "convolution")
par(mfrow = c(1,2))
plot(x); acf(x, main = "", na.action = na.pass )
# way 2
x = arima.sim( model = list( ma = c(.75, .5) ), n )
par(mfrow=c(1,2))
plot(x); acf(x, main = "")

Pi = astsa::ARMAtoAR( ma = c(.75,.5), lag.max = 50)
plot( Pi, type = "h" )
# Repeat the previous steps for MA(2)
x = filter( rnorm(n), c(1, .75, -.5), side = 1, method = "convolution")
par(mfrow=c(1,2))
plot(x); acf(x, main = "", na.action = na.pass )
# method 2
x = arima.sim( model = list( ma = c(.75, -.5) ), n )
par(mfrow=c(1,2))
plot(x); acf(x, main = "")

# ARMA to AR function
Pi = astsa::ARMAtoAR( ma = c(.75, -.5), lag.max = 50)
plot( Pi, type = "h" )


#4
h_max = 30; h = 0:h_max
ACF = ARMAacf( ar = c(-1.6, -.64), lag.max = h_max )
plot( h, ACF, type = "h", main = "(a)" ); abline(h=0, lty=2)

ACF = ARMAacf( ar = c(.4,.45), lag.max = h_max )
plot( h, ACF, type = "h", main = "(b)" ); abline(h=0, lty=2)

ACF = ARMAacf( ar = c(1.2, -.85), lag.max = h_max )
plot( h, ACF, type = "h", main = "(c)" ); abline(h=0, lty=2)


x_AR = arima.sim( list(ar=.6), 100)
x_MA = arima.sim( list(ma=.9), 100)
x_ARMA = arima.sim( list(ar=.6, ma=.9), 100)

par(mfrow=c(1,2))
acf(x_AR, ylim = c(-1,1)); pacf(x_AR, ylim = c(-1,1))
acf(x_MA, ylim = c(-1,1)); pacf(x_MA, ylim = c(-1,1))
acf(x_ARMA, ylim = c(-1,1)); pacf(x_ARMA, ylim = c(-1,1))

X = arima.sim( list(ar=.75), 1000)
Y = arima.sim( list(ar=-.75), 1000)
Z = X + Y
par(mfrow = c(1,2))
acf(Z); pacf(Z)