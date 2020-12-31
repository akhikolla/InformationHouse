
#' @keywords internals
is.rf <- function(x)inherits(x, 'randomForest')


#' @keywords internals
is.ranger <- function(x)inherits(x, 'ranger')

#' @keywords internals

is.parsnip <- function(x)inherits(x, 'model_fit')
