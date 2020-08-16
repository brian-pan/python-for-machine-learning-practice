# F19 FINAL
rm(list=ls())
load("STAD57_W19_Final.RData")

# Q4
set.seed(12345)
n = 100; h = 1:n - 1
X = cumsum( rnorm(n))
acf( X, demean = FALSE, type = "covariance", lag.max = n-1 )
lines( h, (n-h) * (n-h+1) / 2 / n, col = 2 )

N=1000
acf.mat = matrix(0,N,n)
for(i in 1:N){
  Xt = cumsum( rnorm(n))
  acf.mat[i,] = acf( Xt, demean = FALSE, type = "covariance", lag.max = n-1, plot = FALSE )$acf
}
plot( h, colMeans(acf.mat), type = "h", xlab = "mean ACF" )
lines( h, (n-h) * (n-h+1) / n / 2, col = 2 )
lines( h, colMeans(acf.mat) + 2*apply(acf.mat, 2, sd) / sqrt(N), col = 3)
lines( h, colMeans(acf.mat) - 2*apply(acf.mat, 2, sd) / sqrt(N), col = 3)
legend( "topright", legend = c("theoretical", "95% CI"), col = 2:3, lwd=2 )

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(tourists)
acf(tourists, 60); pacf(tourists, 60)

d12tourists = diff(tourists, lag = 12)
acf(d12tourists, 60); pacf(d12tourists, 60)

out = arima( tourists, order = c(3,0,1), seasonal = list( order = c(0,1,2), period = 12 ) ) # excludes mean

# alternatively, can use forecast::Arima
library(forecast)
out2 = Arima( tourists, order = c(3,0,1), seasonal = list( order = c(0,1,2), period = 12 ), include.drift = T ) # includes mean

# or astsa::sarima()
library(astsa)
out3 = sarima( tourists, 3,0,1, 0,1,2, 12, details = F ) # includes mean

mean( abs( out$residuals ) )
summary(out2)
mean( abs(out3$fit$residuals) )

plot( out$residuals )
acf(out$residuals)
qqnorm(out$residuals); qqline(out$residuals)

# astsa::sarima plots diagnostics by default
sarima( tourists, 3,0,1, 0,1,2, 12 )