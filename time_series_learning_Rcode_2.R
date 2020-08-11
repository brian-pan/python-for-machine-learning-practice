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