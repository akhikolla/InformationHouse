#' Variable Selection Frequencies
#'
#' Extract variable selection frequencies from `randomForest` and `ranger` model objects
#'
#' @param x a `randomForest` or `ranger` object
#' @return `tibble` of variable selection frequencies
#'
#' @export
#' @importFrom tibble tibble
#' @examples
#' library(randomForest)
#' data(iris)
#' iris.rf <- randomForest(iris[,-5], iris[,5], forest = TRUE)
#'
#' iris.freqs <- selection_freqs(iris.rf)
#' print(iris.freqs)


selection_freqs <- function(x) {
  if (is.rf(x) == FALSE &
      is.ranger(x) == FALSE & is.parsnip(x) == FALSE) {
    stop(
      deparse(substitute(x)),
      " is not a valid randomForest, ranger or parsnip object",
      call. = FALSE
    )
  }


  if (is.parsnip(x) == FALSE) {
    if (is.null(x$forest)) {
      stop(deparse(substitute(x)), " has no forest", call. = FALSE)
    }
  }



  if(any(class(x) == "randomForest")){
    var <- numeric(length(x$forest$ncat))
    var_freqs <- table(x$forest$bestvar[x$forest$bestvar > 0])
    var[as.numeric(names(var_freqs))] <- var_freqs
    var_df <- tibble(variable = names(x$forest$xlevels), freq = var)
  }


  if(any(class(x) == "ranger")){
    var <- numeric(length(x$forest$independent.variable.names))
    var_freqs <-
      table(unlist(x$forest$split.varIDs)[unlist(x$forest$split.varIDs) > 0])
    var[as.numeric(names(var_freqs))] <- var_freqs
    var_df <-
      tibble(variable = x$forest$independent.variable.names,
             freq = var)
  }



  if (class(x)[1] == '_ranger' & class(x)[2] == 'model_fit') {
    var <- numeric(length(x$fit$forest$independent.variable.names))
    var_freqs <-
      table(unlist(x$fit$forest$split.varIDs)[unlist(x$fit$forest$split.varIDs) > 0])
    var[as.numeric(names(var_freqs))] <- var_freqs
    var_df <-
      tibble(variable = x$fit$forest$independent.variable.names,
             freq = var)
  }


  if (class(x)[1] == '_randomForest' & class(x)[2] == 'model_fit') {
    var <- numeric(length(x$fit$forest$ncat))
    var_freqs <-
      table(x$fit$forest$bestvar[x$fit$forest$bestvar > 0])
    var[as.numeric(names(var_freqs))] <- var_freqs
    var_df <-
      tibble(variable = names(x$fit$forest$xlevels),
             freq = var)
  }


  return(var_df)
}
