outlierProbs <-
  ## Short form for generic function 
  function(object) UseMethod("outlierProbs")

outlierProbs.robmixglm <- function(object) {
  if (!inherits(object, "robmixglm"))
    stop("Use only with 'robmixglm' objects.\n")
  outliers <- object$prop[,2]
  class(outliers) <- "outlierProbs"
  return(outliers)
}