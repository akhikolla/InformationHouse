fit.robmixglm <- function(x,y,family,offset,gh=norm.gauss.hermite(21),notrials,EMTol,cores,verbose) {

 out <- switch(family,
          gaussian = gaussian.fit.robmixglm(x,y,offset,gh,notrials,EMTol,calcHessian=TRUE,cores,verbose,starting.values=NULL),
          binomial = binomial.fit.robmixglm(x,y,offset,gh,notrials,EMTol,calcHessian=TRUE,cores,verbose,starting.values=NULL),
          poisson = poisson.fit.robmixglm(x,y,offset,gh,notrials,EMTol,calcHessian=TRUE,cores,verbose,starting.values=NULL),
         gamma = gamma.fit.robmixglm(x,y,offset,gh,notrials,EMTol,calcHessian=TRUE,cores,verbose,starting.values=NULL),
         truncpoisson = truncpoisson.fit.robmixglm(x,y,offset,gh,notrials,EMTol,calcHessian=TRUE,cores,verbose,starting.values=NULL),
         nbinom = nbinom.fit.robmixglm(x,y,offset,gh,notrials,EMTol,calcHessian=TRUE,cores,verbose,starting.values=NULL)
 )
  out
}
