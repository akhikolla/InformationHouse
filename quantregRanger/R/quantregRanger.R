##' Creates a quantile regression forest like described in Meinshausen, 2006.
##' @title Quantile Regression with Ranger
##' 
##' @param formula Object of class \code{formula} or \code{character} describing the model to fit.
##' @param data Training data of class \code{data.frame}, \code{matrix} or \code{gwaa.data} (GenABEL).
##' @param params.ranger List of further parameters that should be passed to ranger. 
##'   See \code{\link[ranger]{ranger}} for possible parameters.
##' @author Philipp Probst
##' @references
##' Meinshausen, Nicolai. "Quantile regression forests." The Journal of Machine Learning Research 7 (2006): 983-999.
##' @seealso \code{\link{predict.quantregRanger}}
##' @useDynLib quantregRanger
##' @importFrom Rcpp evalCpp sourceCpp
##' @examples
##' y = rnorm(150)
##' x = cbind(y + rnorm(150), rnorm(150))
##' data = data.frame(x,y)
##' mod = quantregRanger(y ~ ., data = data, params.ranger = list(mtry = 2))
##' predict(mod, data = data[1:5, ], quantiles = c(0.1, 0.5, 0.9))
##' @export
quantregRanger = function(formula = NULL, data = NULL, params.ranger = NULL) {
  cl = match.call()
  cl[[1]] = as.name("quantregRanger")
  if (is.null(params.ranger))
    params.ranger = list()
  qrf = do.call(ranger::ranger, c(list(formula = formula) , list(data = data), params.ranger))
  class(qrf) = c("quantregRanger","ranger")
  qrf[["call"]] = cl
  qrf[["origNodes"]] = getnodes(qrf, data)
  qrf[["origObs"]] = stats::model.frame(formula, data)[[1]]
  importance = FALSE
  quantiles = c(0.1, 0.5, 0.9)
  if(importance == TRUE){
    qrf[["importance"]] = predict.imp(qrf, quantiles = quantiles, formula = formula, data = data)
    qrf[["quantiles"]] = quantiles
  }
  qrf
}
