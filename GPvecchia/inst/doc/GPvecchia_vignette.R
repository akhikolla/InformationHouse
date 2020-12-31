## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim=c(7,5)
)

## ------------------------------------------------------------------------
library(GPvecchia)
library(Matrix)
library(fields)

## ------------------------------------------------------------------------
set.seed(1988)
spatial.dim=2
n=50
if(spatial.dim==1){
  locs=matrix(runif(n),ncol=1)
} else {
  locs <- cbind(runif(n),runif(n))
}

## ------------------------------------------------------------------------
beta=2
sig2=1; range=.1; smooth=1.5
covparms =c(sig2,range,smooth)
covfun <- function(locs) sig2*MaternFun(fields::rdist(locs),covparms)
nuggets=rep(.1,n)

## ----fig4, out.width = '400px'-------------------------------------------
Om0 <- covfun(locs)+diag(nuggets)
z=as.numeric(t(chol(Om0))%*%rnorm(n))
data=z+beta

# plot simulated data
if(spatial.dim==1) {
  plot(locs,data)
} else {
  fields::quilt.plot(locs,data, nx=n, ny=n)
}

## ------------------------------------------------------------------------
n.p=100
if(spatial.dim==1){  #  1-D case
  locs.pred=matrix(seq(0,1,length=n.p),ncol=1)
} else {   # 2-D case
  grid.oneside=seq(0,1,length=round(sqrt(n.p)))
  locs.pred=as.matrix(expand.grid(grid.oneside,grid.oneside)) # grid of pred.locs
}
n.p=nrow(locs.pred)

## ------------------------------------------------------------------------
preds=vecchia_pred(vecchia.est,locs.pred)

## ------------------------------------------------------------------------
##  exact prediction
mu.exact=as.numeric(covfun(rbind(locs,locs.pred))[,1:n]%*%solve(Om0,data-beta))+beta
cov.exact=covfun(rbind(locs,locs.pred))-
  covfun(rbind(locs,locs.pred))[,1:n]%*%solve(Om0,t(covfun(rbind(locs,locs.pred))[,1:n]))
var.exact=diag(cov.exact)
cov.exact.pred=cov.exact[n+(1:n.p),n+(1:n.p)]


### plot Vecchia and exact predictions
if(spatial.dim==1) {
  plot(locs,z)
  lines(locs.pred,preds$mean.pred,col='blue')
  lines(locs.pred,preds$mean.pred-1.96*sqrt(preds$var.pred),col='blue',lty=2)
  lines(locs.pred,preds$mean.pred+1.96*sqrt(preds$var.pred),col='blue',lty=2)
  lines(locs.pred,mu.exact[n+(1:n.p)],col='red')
  lines(locs.pred,mu.exact[n+(1:n.p)]-1.96*sqrt(var.exact[n+(1:n.p)]),col='red',lty=2)
  lines(locs.pred,mu.exact[n+(1:n.p)]+1.96*sqrt(var.exact[n+(1:n.p)]),col='red',lty=2)
} else {
  sdrange=range(sqrt(c(preds$var.pred,var.exact[n+(1:n.p)])))
  defpar = par(mfrow=c(2,3))
  fields::quilt.plot(locs,z, nx=sqrt(n.p), ny=sqrt(n.p))
  fields::quilt.plot(locs.pred,preds$mean.pred, nx=sqrt(n.p), ny=sqrt(n.p))
  fields::quilt.plot(locs.pred,sqrt(preds$var.pred),zlim=sdrange, nx=sqrt(n.p), ny=sqrt(n.p))
  fields::quilt.plot(locs,z, nx=sqrt(n.p), ny=sqrt(n.p))
  fields::quilt.plot(locs.pred,mu.exact[n+(1:n.p)], nx=sqrt(n.p), ny=sqrt(n.p))
  fields::quilt.plot(locs.pred,sqrt(var.exact[n+(1:n.p)]),zlim=sdrange, nx=sqrt(n.p), ny=sqrt(n.p))
  par(defpar)
}

## ------------------------------------------------------------------------
m=20
vecchia.approx=vecchia_specify(locs,m)
vecchia_likelihood(z,vecchia.approx,covparms,nuggets)

## ------------------------------------------------------------------------
library(mvtnorm)
dmvnorm(z,mean=rep(0,n),sigma=Om0,log=TRUE)

## ------------------------------------------------------------------------
m=30
vecchia.approx=vecchia_specify(locs,m,locs.pred=locs.pred)
preds=vecchia_prediction(z,vecchia.approx,covparms,nuggets)
# returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord

## ------------------------------------------------------------------------
Sigma=V2covmat(preds)$Sigma.pred
cov.range=quantile(rbind(Sigma,cov.exact.pred),c(.01,.99))
defpar = par(mfrow=c(1,2))
fields::image.plot(cov.exact.pred,zlim=cov.range)
fields::image.plot(Sigma,zlim=cov.range)
par(mfrow=c(defpar))

## ------------------------------------------------------------------------
H=Matrix::sparseMatrix(i=1:(n+n.p),j=1:(n+n.p),x=1)[(n+1):(n+n.p),]

# compute variances of Hy
lincomb.vars=vecchia_lincomb(H,preds$U.obj,preds$V.ord)
plot(preds$var.pred,lincomb.vars)

## ------------------------------------------------------------------------
mean(preds$mu.pred)

# compute entire covariance matrix of Hy (here, 1x1)
H=Matrix::sparseMatrix(i=rep(1,n.p),j=n+(1:n.p),x=1/n.p)
lincomb.cov=vecchia_lincomb(H,preds$U.obj,preds$V.ord,cov.mat=TRUE)

## ------------------------------------------------------------------------
m=20
mra.options.fulls=list(M=1)
blockFS = vecchia_specify(locs, m, 'maxmin', conditioning='mra', mra.options=mra.options.fulls, verbose=TRUE)

## ------------------------------------------------------------------------
mra.options.mpproc=list(r=c(m,1))
MPP = vecchia_specify(locs, m, 'maxmin', conditioning='mra', mra.options=mra.options.mpproc, verbose=TRUE)

## ------------------------------------------------------------------------
mra.options.mra = list(r=c(10, 5, 5), M=2, J=2)
MRA_rJM = vecchia_specify(locs, m, 'maxmin', conditioning='mra', mra.options=mra.options.mra, verbose=TRUE)

## ------------------------------------------------------------------------
NNGP = vecchia_specify(locs, m, cond.yz='y')

## ------------------------------------------------------------------------
vecchia_likelihood(z,blockFS,covparms,nuggets)
vecchia_likelihood(z,MPP,covparms,nuggets)
vecchia_likelihood(z,MRA_rJM,covparms,nuggets)
vecchia_likelihood(z,NNGP,covparms,nuggets)
vecchia_likelihood(z, vecchia_specify(locs, m), covparms, nuggets)
dmvnorm(z,mean=rep(0,n),sigma=Om0,log=TRUE)

## ------------------------------------------------------------------------
# simulate latent process
y=as.numeric(t(chol(Om0))%*%rnorm(n))

## ------------------------------------------------------------------------
data.model = "logistic"

# simulate data
if(data.model=='poisson'){
  z = rpois(n, exp(y))
} else if(data.model=='logistic'){
  z = rbinom(n,1,prob = exp(y)/(1+exp(y)))
} else if(data.model=='gamma'){
  z = rgamma(n, shape = default_lh_params$alpha, rate = default_lh_params$alpha*exp(-y))
}else{
  print('Error: Distribution not implemented yet.')
}

# plot simulated data, 1 or 2D
defpar = par(mfrow=c(1,2))
if(spatial.dim==1) {
  plot(locs,y, main = "latent")
  plot(locs,z, main = "observed")
} else {
  fields::quilt.plot(locs,y, main = "Latent")
  fields::quilt.plot(locs,z, main = "Observed")
}
par(defpar)


## ------------------------------------------------------------------------
m=10
if(spatial.dim==1){
  vecchia.approx=vecchia_specify(locs,m) #IW ordering
} else {
  vecchia.approx=vecchia_specify(locs,m,cond.yz='zy') #RF ordering
}

## ------------------------------------------------------------------------
posterior = calculate_posterior_VL(z,vecchia.approx,likelihood_model=data.model,
                                   covparms = covparms)
if (spatial.dim==1){
  par(mfrow=c(1,1))
  ord = order(locs) # order so that lines appear correctly
  y_limits = c(min(y, posterior$mean[ord]), max(y, posterior$mean[ord]))
  plot(locs[ord], y[ord], type = "l", ylim = y_limits )
  lines(locs[ord], posterior$mean[ord], type = "l", col=3, lwd=3)
  legend("bottomright", legend = c("Latent", "VL"), col= c(1,3), lwd=c(1,3))
} else if (spatial.dim==2){
  dfpar = par(mfrow=c(1,2))
  # ordering unnecessary; we are using a scatter plot rather than lines
  quilt.plot(locs, y, main= "Truth")
  quilt.plot(locs, posterior$mean,  main= "VL m=10")
  par(defpar)
  
}


## ------------------------------------------------------------------------
######  specify prediction locations   #######
n.p=30^2
if(spatial.dim==1){  #  1-D case
  locs.pred=matrix(seq(0,1,length=n.p),ncol=1)
} else {   # 2-D case
  grid.oneside=seq(0,1,length=round(sqrt(n.p)))
  locs.pred=as.matrix(expand.grid(grid.oneside,grid.oneside)) # grid of pred.locs
}
n.p=nrow(locs.pred)

######  specify Vecchia approximation   #######
vecchia.approx.pred = vecchia_specify(locs, m=10, locs.pred=locs.pred)
###  carry out prediction
preds = vecchia_laplace_prediction(posterior, vecchia.approx.pred, covparms)

# plotting predicitions
if (spatial.dim==1){
  defpar = par(mfrow=c(1,1))
  ord = order(locs) # order so that lines appear correctly
  plot(locs[ord], y[ord], type = "l", xlim=c(0,1.2), ylim = c(-1,3))
  lines(locs, posterior$mean, type = "p", col=4, lwd=3, lty=1)
  lines(locs.pred, preds$mu.pred, type = "l", col=3, lwd=3, lty=1)
  lines(locs.pred,preds$mu.pred+sqrt(preds$var.pred), type = "l", lty = 3, col=3)
  lines(locs.pred,preds$mu.pred-sqrt(preds$var.pred), type = "l", lty = 3, col=3)
  legend("topleft", legend = c("Latent", "VL: Pred", "VL: 1 stdev"), 
         col= c(1,3,3), lwd=c(1,2,1), lty = c(1,1,3))
  par(defpar)
} else if (spatial.dim==2){
  defpar =  par(mfrow=c(1,2))
  # ordering unnecessary; we are using a scatter plot rather than lines
  quilt.plot(locs, y, main= "True Latent", 
             xlim = c(0,1), ylim = c(0,1), nx=64, ny=64)
  quilt.plot(locs.pred, preds$mu.pred,  main= "VL Prediction",nx = 30, ny=30)
  par(defpar)
}

## ------------------------------------------------------------------------
vecchia_laplace_likelihood(z,vecchia.approx,likelihood_model=data.model,covparms = covparms)

## ---- eval = FALSE-------------------------------------------------------
#  # currently set up for covariance estimation
#  vecchia.approx=vecchia_specify(locs, m=10, cond.yz = "zy") # for posterior
#  vecchia.approx.IW = vecchia_specify(locs, m=10) # for integrated likelihood
#  if (spatial.dim==1) vecchia.approx=vecchia.approx.IW
#  
#  vl_likelihood = function(x0){
#    theta = exp(x0)
#    covparms=c(theta[1], theta[2], theta[3]) # sigma range smoothness
#    prior_mean = 0 # can be a parameter as well
#    # Perform inference on latent mean with Vecchia Laplace approximation
#    vll = vecchia_laplace_likelihood(z,vecchia.approx, likelihood_model=data.model,
#                                     covparms, return_all = FALSE,
#                                     likparms = default_lh_params, prior_mean = prior_mean,
#                                     vecchia.approx.IW = vecchia.approx.IW)
#    return(-vll)
#  
#  }
#  x0 = log(c(.07,1.88, 1.9))
#  vl_likelihood(x0)
#  # Issues with R aborting, maxit set to 1
#  res = optim(x0, vl_likelihood, method = "Nelder-Mead", control = list("trace" = 1, "maxit" = 1))
#  exp(res$par[1:3])
#  vl_likelihood(x0)

