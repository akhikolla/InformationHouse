negll <- function(p, y, x, offset) {
  lp <- x %*% matrix(p[1:(length(p)-1)],ncol=1)+offset
  theta <- p[length(p)]
  negll <- -sum(dnbinom(y,mu=exp(lp),size=theta,log=TRUE))
  if (negll>1.0e100) negll <- 1.0e100
  if (!is.finite(negll)) browser()
  return(negll)
} 

fitnegbin <- function(y,x,offset) {
  # obtain starting values
  pois.glm <- glm(y~x[,-1], family=poisson, offset=offset)
  initp <- c(coef(pois.glm),1)
  # fit <- optim(initp, negll, method =  "L-BFGS-B",
  #              lower = c(rep(-Inf,length(initp)-1),0),  hessian = FALSE, x=x, y=y, offset=offset)
  fit <- nlminb(initp, negll, lower = c(rep(-Inf,length(initp)-1),0),   x=x, y=y, offset=offset)
  if (fit$convergence>0) browser()
  return(fit)
}
