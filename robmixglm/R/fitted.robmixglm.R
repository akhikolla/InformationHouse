fitted.robmixglm <- function(object, ...) {
  if (!inherits(object, "robmixglm"))
    stop("Use only with 'robmixglm' objects.\n")
    lp <- object$X %*% matrix(coef(object),ncol=1) + object$offset
  out <- switch (
    object$family,
    gaussian = lp,
    binomial = 1.0/(1.0+exp(-lp)),
    poisson = exp(lp),
    gamma = exp(lp),
    truncpoisson = exp(lp),
    nbinom = exp(lp)
  )
  return(out)
}
