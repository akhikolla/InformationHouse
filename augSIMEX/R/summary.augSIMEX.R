summary.augSIMEX <-function (object, dispersion=NULL,...) 
{   
  ans <- list()
  est.disp <- FALSE
  df.r <- object$df.residual
  if (is.null(dispersion))	# calculate dispersion if needed
  {dispersion <-if(object$family$family %in% c("poisson", "binomial"))  1}
  else if(df.r > 0) {
    est.disp <- TRUE
    if(any(object$weights==0))
      warning("observations with zero weight not used for calculating dispersion")
    sum((object$weights*object$residuals^2)[object$weights > 0])/ df.r
  } else {
    est.disp <- TRUE
    NaN
  }
  
  ans$aliased <- is.na(object$coefficients)
  
  p.names <- names(object$coefficients)
  est <- object$coefficients
  est.table <- list()
  
  se <- sqrt(diag(object$vcov))
  tvalue <-est/se
  pvalue<- 2 * pnorm(-abs(tvalue))
  est.table <- cbind(est, se, tvalue, pvalue)
  dimnames(est.table) <- list(p.names, c("Estimate","Std. Error","t value","Pr(>|t|)"))
  
  keep <- match(c("call","family","deviance", "aic",
                  "df.residual","null.deviance","df.null",
                  "na.action"), names(object), 0L)
  
  ans$coefficients <- est.table
  ans$extrapolation <- object$extrapolation
  ans$err.var <- object$err.var
  ans$mis.var <- object$mis.var
  ans$err.true <- object$err.true
  ans$mis.true <- object$mis.true
  ans$B<-object$B
  ans$nBoot<-object$nBoot
  ans$dispersion<-dispersion
  ans<-c(ans,object[keep])
  
  class(ans) <- "summary.augSIMEX"
  
  ans
  
}
