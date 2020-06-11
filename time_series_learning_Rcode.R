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