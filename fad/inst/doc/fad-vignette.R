## ----setup, include=FALSE, cache=FALSE,echo=FALSE---------------------------------------
library(knitr)
# set global chunk options
opts_chunk$set(fig.path='figure/fad-', fig.align='center', fig.show='hold')
options(formatR.arrow=TRUE,width=90)
Sys.setenv(RSTUDIO_PDFLATEX = "/Library/TeX/texbin/latexmk")

## ----load-------------------------------------------------------------------------------
library(fad)

## ----fad-arg----------------------------------------------------------------------------
args(fad)

## ----set-seed, include=FALSE------------------------------------------------------------
set.seed(1)

## ----est-sim-data,cache=TRUE------------------------------------------------------------
n <- 50 # number of observations
p <- 100 # size of coordinates
q <- 3 # true number of factors
mu = rnorm(p) 
L = matrix(rnorm(p*q),p,q)
D = runif(p,0.2,0.8)
X = matrix(rnorm(n*q),n,q)
X = tcrossprod(X,L) + matrix(rnorm(n*p),n,p)%*%diag(sqrt(D)) + rep(mu, each = n)

BICs = rep(0,11) # store BIC values
for(i in 1:11){
  out = fad(x = X,factors = i-1)
  BICs[i] = out$BIC
}
BICs
ind = which.min(BICs)-1# obtain the optimal factor favored by BIC
ind

## ----est-sim-cov1,cache=TRUE------------------------------------------------------------
BICs = rep(0,11) # store BIC values
for(i in 1:11){
  out = fad(covmat = cov(X), n.obs = n, factors = i-1)
  BICs[i] = out$BIC
}
BICs
ind = which.min(BICs)-1
ind

## ----est-sim-cov2,cache=TRUE------------------------------------------------------------
out1 = fad(x = X,factors = 3)
# or
out2 = fad(covmat = cov(X), n.obs = n, factors = 3)

