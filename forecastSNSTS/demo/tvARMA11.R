
## In this demo we demonstrate how the data analysis in Section 5 of
## Kley et al (2017) can be performed, but with artificial data.

## To this end, we generate artificial data for the example analysis.
## We use a tvARMA(1,1) process with the following coefficient functions:
a     <- list( function(u) { return(u - 0.5) } )
b     <- list( function(u) { return( 1 + sin(pi * u) ) } )
sigma <- function(u) { return( 1 + 3 * u ) }

## The function tvARMA from this package can then be used to simulate
X <- tvARMA(T = 250, a, b, sigma)

## Show what the path looks like.
plot(X, type="l")

## We now choose the maximum model order p_max and the size of the
## validation and end of the training set m. 
p_max <- 7
m <- 20

## We also choose the maximum forecasting horizon
H <- 6

## Given the length L of the time series we simulated
L <- length(X)

## We have
## the end point m0 of the main part of the training set,
## the end point m1 of the final part of the training set,
## the end point m2 of the validation set and
## the end point m3 of the test set: 
m0 <- L - 3*m # 190
m1 <- L - 2*m # 210
m2 <- L - m   # 230
m3 <- L       # 250

## define the set \mathcal{N}; i.e., which N to compute the predictions for
Ns <- c(0,seq(min(p_max+1, ceiling((L/2)^(4/5))), ceiling(L^(4/5))))

## compute all the linear prediction coefficients needed
coef0 <- predCoef(X, p_max, H, (m0-H+1):(m3-1), Ns)

## compute the MSPE on the final part of the training set
mspe <- MSPE(X, coef0, m0+1, m1, p_max, H, Ns)

## compute the minima for all N >= N_min
N_min <- 35
M <- mspe$mspe
N <- mspe$N

res_e <- matrix(0, nrow = H, ncol = 5)
for (h in 1:H) {
  res_e[h, 1] <- idx1_s <- which(M[h, , N == 0] == min(M[h, , N == 0]),
                                 arr.ind = TRUE)[1]
  res_e[h, 2] <- min(M[h, , N == 0])
  idx1_ls <- which(M[h, , N != 0 & N >= N_min] == min(M[h, , N != 0 & N >= N_min]),
                   arr.ind = TRUE)[1,]
  res_e[h, 3] <- idx1_ls[1]
  res_e[h, 4] <- N[N != 0 & N >= N_min][idx1_ls[2]]
  res_e[h, 5] <- min(M[h, , N != 0 & N >= N_min])
}

## A table similar to the top row from Table 5 in Kley et al (2017):
res_e

## compute the MSPE of the null predictor
vr <- sum(X[(m0 + 1):m1]^2) / (m1 - m0)

## A plot similar to the top plot from Figure 4 in Kley et al (2017)
plot(mspe, vr = vr, N_min = N_min, h = 1, add.for.legend=15)

## A plot similar to the bottom plot from Figure 4 in Kley et al (2017)
plot(mspe, vr = vr, N_min = N_min, h = 6, add.for.legend=15)

## compute MSPE on the validation set 
mspe <- MSPE(X, coef0, m1 + 1, m2, p_max, H, Ns)
M <- mspe$mspe
N <- mspe$N

## Compare the stationary approach with the locally stationary approach,
## both for the optimally chosen p_s and (p_ls, N_ls) that achieved the
## minimal MSPE on the final part of the training set.

res_v <- matrix(0, nrow = H, ncol = 3)
for (h in 1:H) {
  res_v[h, 1] <- M[h, res_e[h, 1], N == 0]
  res_v[h, 2] <- M[h, res_e[h, 3], N == res_e[h, 4]]
  res_v[h, 3] <- res_v[h, 1] / res_v[h, 2]
}

## compute MSPE on the validation set 
mspe <- MSPE(X, coef0, m2 + 1, m3, p_max, H, Ns)
M <- mspe$mspe
N <- mspe$N

## Compare the stationary approach with the locally stationary approach,
## both for the optimally chosen p_s and (p_ls, N_ls) that achieved the
## minimal MSPE on the final part of the training set.

res_t <- matrix(0, nrow = H, ncol = 3)
for (h in 1:H) {
  res_t[h, 1] <- M[h, res_e[h, 1], N == 0]
  res_t[h, 2] <- M[h, res_e[h, 3], N == res_e[h, 4]]
  res_t[h, 3] <- res_t[h, 1] / res_t[h, 2]
}

## A table similar to the bottom rows from Table 5 in Kley et al (2017)
cbind(res_v, res_t)
