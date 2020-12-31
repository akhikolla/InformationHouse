## ----setup, include = FALSE---------------------------------------------------
library(cort)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----clayton-frailty-model----------------------------------------------------
psi <- function(t,alpha) (1 + sign(alpha)*t) ^ (-1/alpha) # generator
rClayton <- function(n,dim,alpha){
  val <- matrix(runif(n * dim), nrow = n)
  gam <- rgamma(n, shape = 1/alpha, rate = 1)
  gam <- matrix(gam, nrow = n, ncol = dim)
  psi(- log(val) / gam,alpha)
}

## ----dataset------------------------------------------------------------------
set.seed(12)
n = 200 # taken small to reduce runtime of the vignette.
d = 4
n_trees = 5 # taken small to reduce runtime of the vignette.
number_max_dim_forest = 2 # taken small to reduce runtime of the vignette.

data <- matrix(nrow=n,ncol=d)
data[,c(1,4,3)] = rClayton(n=n,dim=d-1,alpha=7)
data[,2] = runif(n)
data[,3] <- 1 - data[,3]

pairs(data,cex=0.6)
  

## ----run_cort-----------------------------------------------------------------
set.seed(12)
(model = Cort(data,verbose_lvl = 1))

## ----fig.cap="Pairs-plot of original data (in black, bottom-left corner) versus a simulation from the model (in red, top-right corner)"----
pairs(model)

## ----fig.cap="Gray boxes representing 2-d projections of the fitted density. In red, the imputed data points."----
plot(model)

