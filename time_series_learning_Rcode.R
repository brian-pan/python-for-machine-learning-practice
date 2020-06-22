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

# ar
X = arima.sim( list(ar=.75), 1000)
Y = arima.sim( list(ar=-.75), 1000)
Z = X + Y
par(mfrow = c(1,2))
acf(Z); pacf(Z)
# ma
X = arima.sim( list(ma=.75), 1000)
Y = arima.sim( list(ma=-.75), 1000)
Z = X + Y
par(mfrow = c(1,2))
acf(Z); pacf(Z)


#5
(out = arima(x, order = c(1,0,1), include.mean = F) )
(pred= predict(out, n.ahead = 1))

phi = out$coef[1]; theta = out$coef[2]
theor_acf = ARMAacf(ar = phi, ma = theta, lag.max = 100)
DL = acf2AR( theor_acf )
phi_100 = DL[100,]

(my_pred = sum( phi_100 * rev(x) ))

my_pred - as.numeric(pred$pred)

(gamma0 = (1 + 2*phi*theta + theta^2) /(1-phi^2)*out$sigma2)
(my_se = sqrt( gamma0 * prod( 1 -  diag(DL)^2 ) ))
(my_se = sqrt( gamma0 * ( 1 - sum( phi_100 * theor_acf[-1] ) ) ))

my_se - as.numeric(pred$se)

library(forecast)
(fore = forecast(out, h = 30))
plot(fore)

sqrt( var(x) )*1.96


#6
library(astsa)
library(magrittr)
plot(cmort)

## Yule Walker
fit_yw = ar.yw(cmort, order.max = 2)
## Linear Regression
p = 2 # AR order
n = length(cmort) # series length
Y = cmort[(p+1):n] # IV 
X = matrix(0,(n-p),p) # DVs placeholder
for(i in 1:p){
  X[,i] = cmort[(p+1-i):(n-i)] # lagged values of Y (by lag i)
}
fit_lm = lm(Y ~ X)

## (a) 

# AR coef's
fit_yw$ar
fit_lm$coefficients[-1]

# Series mean \mu
fit_yw$x.mean
# corresponding intercept (see SS p.103)
fit_yw$x.mean * ( 1 - sum(fit_yw$ar) )
# lm intercept
fit_lm$coef[1]

# sigma_w^2
fit_yw$var.pred  
sigma(fit_lm)^2

## (b)
# from fit_yw$asy.var.coef 
fit_yw$asy.var.coef %>% diag() %>% sqrt()

# from summary(fit_lm)
summary(fit_lm)$coef[-1,2]



set.seed(123)
for(i in 1:5){
  X = arima.sim( model = list( ar = .9, ma = -.9 ), n = 200 )
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  plot(X); acf(X); pacf(X)
  print( arima(X, order = c(1,0,1) ) )
}

set.seed(123)
for(i in 1:10){
  X = arima.sim( model = list( ar = .9, ma = .5 ), n = 200 )
  print( arima(X, order = c(1,0,1) ) )
}



library(astsa)
gnpgr = diff(log(gnp)) # GNP growth rate
sarima( gnpgr, 1, 0, 0)

# Fit AR(1) 
fit_ar = arima(gnpgr ,order=c(1,0,0)) 
# or fit_ar = ar.mle(gnpgr ,order=1) 

stdres = fit_ar$residuals / sqrt(fit_ar$sigma2)

layout( matrix( c(1,1,2,3,4,4), 3, 2, byrow = TRUE) )

plot(stdres) # residual plot
qqnorm(stdres); qqline(stdres) # normal QQ-plot
acf(stdres) # residual ACF
# Ljung-Box test (H=20) 
lags=3:20; p.values=rep(0,length(lags))
for(i in 1:length(lags)){
  p.values[i]=Box.test(stdres, lags[i], type = "Ljung-Box", fitdf = 2)$p.value
}
plot(lags, p.values, ylim=c(0,1), main="Ljung-Box p-values"); abline(h=.05, lty=2)


#7
library(astsa)
plot(globtemp)
astsa::acf2(globtemp) # neat ACF/PACF wrapper (Note: ACF starts @ 1)

# RW-type behavior, perform a unit-root test:
library(tseries)
adf.test(globtemp)
acf2(diff(globtemp))

(out_AR3 = sarima(globtemp, 3,1,0))
(out_MA4 = sarima(globtemp, 0,1,4))

sarima.for(globtemp, 3,1,0, n.ahead = 10)  
out_AR3$fit$coef["constant"] 

out_auto = forecast::auto.arima(globtemp)
AIC(out_auto)
AIC(out_AR3$fit)  

sarima(globtemp, 1,1,3)
sarima.for(globtemp, 1,1,3, n.ahead = 10)

plot( ARMAacf(ar = c(rep(0,11), .9), ma = .5, lag.max = 90 ),
      type = 'h', lwd = 2, ylab = "ACF", xlab = "lag" )

plot(chicken)
acf2(chicken)
plot(diff(chicken))
acf2(diff(chicken))
library(forecast)
auto.arima(chicken)
sarima.for(chicken, 2,1,1,0,0,1,12, n.ahead = 12)

ljj = log(jj)
plot( cbind(jj,ljj), type = 'o', pch = 20)
auto.arima( ljj)

sarima(ljj,2,0,0,1,1,0,4) # fit w/ sarima for diagnostics 

# log-scale forecasts 
log_forecasts = sarima.for( ljj, 2,0,0,1,1,0,4, n.ahead = 12 )

# Original scale forecasts
log_upper = log_forecasts$pred + 1.96 * log_forecasts$se
log_lower = log_forecasts$pred - 1.96 * log_forecasts$se
U = exp( log_upper )
L = exp( log_lower )

plot(jj, xlim = c(1960,1984), ylim = c(0,35), type = 'o', main = "Original Scale")
lines( exp( log_forecasts$pred ), col = 2, type = 'o' )
xx = c(time(U), rev(time(U)))
yy = c(L, rev(U))
polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))


#8
library(astsa)
plot(sales)
acf2(sales)
tseries::adf.test(sales) # ADF test doesn't reject RW
plot(diff(sales))
acf2(diff(sales))

library(forecast)
auto.arima(sales)
library(astsa)
sarima( sales, p = 1, d = 1, q = 1 , details = F )
ccf( diff(sales), diff(lead) )