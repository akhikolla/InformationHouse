summary.robmixglm <- function(object, ...) {
  if (!inherits(object, "robmixglm"))
    stop("Use only with 'robmixglm' objects.\n")
  out <- list()
  out$logLik <- logLik(object)
  out$AIC <- AIC(object)
  out$BIC <- BIC(object)
  coefficients <- matrix(object$fit@coef,ncol=1)
  coefficients[,1] <- ifelse(!(object$coef.names=="Outlier p."),coefficients[,1],1.0/(1.0+exp(-coefficients)))
  lastnames <- 2
  if (object$family %in% c("gaussian","gamma","nbinom")) lastnames <- 3
  coef.se <- sqrt(diag(object$fit@vcov))
  coef.se[(length(coefficients[,1])-lastnames+1):length(coefficients[,1])] <- NA
   coefficients <- cbind(coefficients,coef.se)
  
  coefficients <- cbind(coefficients,coefficients[,1]/coefficients[,2])
  coefficients <- cbind(coefficients,2*(1-pnorm(abs(coefficients[,3]))))
  dimnames(coefficients)[[1]] <- object$coef.names
  dimnames(coefficients)[[2]] <- c("Estimate","Std. Error","z value","Pr(>|z|)")
  out$coefficients <- coefficients
  class(out) <- "summary.robmixglm"
  out
}