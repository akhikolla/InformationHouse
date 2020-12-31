#' Extract forest parameters
#'
#' For a `randomForest` or `ranger` classification object, extract the parameters needed to calculate an approximate selection frequency threshold
#'
#' @param x a `randomForest`, `ranger` or `parsnip` object
#' @return a list of four elements
#' * __Fn__ The number of features considered at each internal node (mtry)
#' * __Ft__ The total number of features in the data set
#' * __K__ The average number of binary tests/internal nodes across the enitre forest
#' * __Tr__ The total number of trees in the forest
#'
#'
#' @author Tom Wilson \email{tpw2@@aber.ac.uk}
#' @export
#' @examples
#' library(randomForest)
#' data(iris)
#' iris.rf <- randomForest(iris[,-5], iris[,5], forest = TRUE)
#'
#' iris.params <- extract_params(iris.rf)
#' print(iris.params)

extract_params <- function(x)
  {

  if(is.rf(x) == FALSE & is.ranger(x) == FALSE & is.parsnip(x) == FALSE){
    stop(deparse(substitute(x)), " is not a valid randomForest, ranger or parsnip object", call. = FALSE)
  }


  if(is.parsnip(x) == FALSE) {
    if (is.null(x$forest)) {
      stop(deparse(substitute(x)), " has no forest", call. = FALSE)
    }
  }



  if(any(class(x) == "randomForest")){

          Fn <- x$mtry
          Fe <- nrow(x$importance)
          K <- round(mean(apply(x$forest$nodestatus,2,function(x)(length(which(x == 1))))), digits = 0)
          Tr <- x$forest$ntree

  }


  if(any(class(x) == "ranger")){

          Fn <- x$mtry
          Fe <- x$num.independent.variables
          K <- round(mean(sapply(x$forest$split.varIDs, function(x)(length(which(x != 0))))), digits = 0)
          Tr <- x$num.trees
  }

  if(class(x)[1] == '_ranger' & class(x)[2] == 'model_fit'){
    Fn <- x$fit$mtry
    Fe <- x$fit$num.independent.variables
    K <- round(mean(sapply(x$fit$forest$split.varIDs, function(x)(length(which(x != 0))))), digits = 0)
    Tr <- x$fit$num.trees
  }


  if(class(x)[1] == '_randomForest' & class(x)[2] == 'model_fit'){
    Fn <- x$fit$mtry
    Fe <- nrow(x$fit$importance)
    K <- round(mean(apply(x$fit$forest$nodestatus,2,function(x)(length(which(x == 1))))), digits = 0)
    Tr <- x$fit$forest$ntree
  }

  return(list(Fn = Fn,Ft = Fe, K = K, Tr = Tr))
  }
