baseline.default <- function(object, ...) {

  match.call()

  if(!(class(object) %in% c("alacoxIC", "unpencoxIC"))) stop("Please input the object returned by 'alacoxIC' or 'unpencoxIC'.")

  lambda <- object$lambda
  nlambda <- length(lambda)
  lambda.set <- object$lambda.set

  result <- data.frame(lambda.set, lambda, cumsum(lambda))
  colnames(result) <- c("lower.set", "upper.set", "lambda", "clambda")

  class(result) <- "baseline"

  return(result)

}
