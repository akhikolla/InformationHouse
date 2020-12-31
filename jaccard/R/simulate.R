sim.coeff.p <- function(nobs, nvar, baseP=.5, d=1, pi0=.9, seed=NULL)  {
  if(!is.null(seed)) set.seed(seed)
  nobs0 = round(nobs*pi0)
  nobs1 = nobs-round(nobs*pi0)

  baseB = log(baseP/(1-baseP))
  B0 = matrix(baseB, nrow=nobs0, ncol=nvar)
  B1 = matrix(runif(nobs1*d, min=1, max=2), nrow=nobs1, ncol=d)
  L = matrix(rnorm(nvar*d), nrow=d, ncol=nvar)
  B1L = B1 %*% L
  BL = rbind(B1L,B0)
  prob = exp(BL)/(1+exp(BL))

  Y = matrix(rbinom(nobs*nvar, 1, as.numeric(prob)), nobs, nvar)
  H = c(rep(1, nobs1), rep(0, nobs0))
  return(list(Y=Y, H=H, prob=prob))
}
