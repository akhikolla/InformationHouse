logLik.robmixglm <- function(object, ...)
{
     val <- object$logLik
    attr(val, "df") <- object$np
    attr(val, "nobs") <- object$nobs
    class(val) <- "logLik"
    val
}

BIC.robmixglm <-
function (object, ...) 
{
    if (!is.element("robmixglm", class(object))) 
        stop("Argument 'object' must be an object of class \"robmixglm\".")
    BIC(logLik(object))
}

AIC.robmixglm <-
function(object,...,k=2) {
    if (!inherits(object, "robmixglm"))
        stop("Use only with 'robmixglm' objects.\n")
    AIC(logLik(object),...,k)
}

extractAIC.robmixglm <- function(fit, scale = 0, k = 2, ...)
{
  if (!inherits(fit, "robmixglm"))
    stop("Use only with 'robmixglm' objects.\n")
  n <- fit$nobs
  edf <- n  - fit$np # assumes dispersion is known
  aic <- AIC(fit)
  c(edf, aic + (k-2) * edf)
}
