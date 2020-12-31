##' Predicts quantiles for a quantile regression forest trained with quantregRanger. 
##' @title quantregRanger prediction
##'
##' @param object \code{quantregRanger} object.
##' @param data New test data of class \code{data.frame}
##' @param quantiles Numeric vector of quantiles that should be estimated
##' @param all A logical value. all=TRUE uses all observations for prediction. 
##' all=FALSE uses only a certain number of observations per node for prediction 
##' (set with argument obs). The default is all=TRUE
##' @param obs An integer number. Determines the maximal number of observations per node 
##' @param ... Currently ignored.
##' to use for prediction. The input is ignored for all=TRUE. The default is obs=1
##' @return A matrix. The first column contains the conditional quantile 
##' estimates for the first entry in the vector quantiles. The second 
##' column contains the estimates for the second entry of quantiles and so on.
##' @export
predict.quantregRanger = function(object, data = NULL, quantiles = c(0.1, 0.5, 0.9), all = TRUE, obs = 1, ...) {
  #check if all is logical
  if (!is.logical(all)) {
    stop("all has to be logical")
  }
  if (all==FALSE){#don't use all observations for prediction
    #check if obs is integer
    if(!{obs%%1==0}) {
      stop("obs has to be an integer")
    }
    predict.fast(object, data, quantiles, obs)
  }else{
    predict.all(object, data, quantiles)
  }
}