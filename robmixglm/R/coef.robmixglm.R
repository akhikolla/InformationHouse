coef.robmixglm <- function(object, ...) {
  if (!inherits(object, "robmixglm"))
    stop("Use only with 'robmixglm' objects.\n")
  res <- switch (
    object$family,
    gaussian = object$fit@coef[1:(length(object$fit@coef)-3)],
    binomial = object$fit@coef[1:(length(object$fit@coef)-2)],
    poisson = object$fit@coef[1:(length(object$fit@coef)-2)],
    gamma = object$fit@coef[1:(length(object$fit@coef)-3)],
    truncpoisson = object$fit@coef[1:(length(object$fit@coef)-2)],
    nbinom = object$fit@coef[1:(length(object$fit@coef)-3)]
    )
    names(res) <- object$coef.names[1:length(res)]
    res
}