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