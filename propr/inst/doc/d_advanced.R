## ------------------------------------------------------------------------
library(propr)
N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

## ------------------------------------------------------------------------
Y <- X / rowSums(X) * abs(rnorm(N))

## ---- dpi = 66, fig.width = 8, fig.height = 8, fig.show = "hold"---------
pairs(X) # absolute data

## ---- dpi = 66, fig.width = 8, fig.height = 8, fig.show = "hold"---------
pairs(Y) # relative data

## ---- warning = FALSE----------------------------------------------------
suppressWarnings(cor(X)) # absolute correlation
cor(Y) # relative correlation

## ------------------------------------------------------------------------
propr:::proprVLR(Y[, 1:4]) # relative VLR
propr:::proprVLR(X) # absolute VLR

## ---- dpi = 66, fig.width = 8, fig.height = 8, fig.show = "hold"---------
pairs(propr:::proprCLR(Y[, 1:4])) # relative clr-transformation

## ---- dpi = 66, fig.width = 8, fig.height = 8, fig.show = "hold"---------
pairs(propr:::proprCLR(X)) # absolute clr-transformation

## ---- message = FALSE----------------------------------------------------
propr(Y[, 1:4])@matrix # relative proportionality with clr

## ---- message = FALSE----------------------------------------------------
propr(X)@matrix # absolute proportionality with clr

## ---- dpi = 66, fig.width = 8, fig.height = 8, fig.show = "hold"---------
pairs(propr:::proprALR(Y, ivar = 5)) # relative alr

## ---- dpi = 66, fig.width = 8, fig.height = 8, fig.show = "hold"---------
pairs(X[, 1:4]) # absolute data

## ---- message = FALSE----------------------------------------------------
propr(Y, ivar = 5)@matrix # relative proportionality with alr

