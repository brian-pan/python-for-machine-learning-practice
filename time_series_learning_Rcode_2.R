# F19 FINAL
rm(list=ls())
load("STAD57_W19_Final.RData")

# Q4
set.seed(12345)
n = 100; h = 1:n - 1
X = cumsum( rnorm(n))
acf( X, demean = FALSE, type = "covariance", lag.max = n-1 )
lines( h, (n-h) * (n-h+1) / 2 / n, col = 2 )