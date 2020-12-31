#' Predictions from a fitted metacart object
#'
#' Returns a data frame of predicted effect sizes and moderators from a fitted metacart object
#'
#' @param object fitted model object of class "FEmrt".
#' @param newdata data frame containing the values at which predictions are required.
#' @param ... Arguments that pass to other methods.
#' @return  A data frame containing the predicted effect size, the moderators, and the corresponding node lables in the fitted tree.
#' @importFrom stats as.formula delete.response model.frame model.response terms predict
#' @export
predict.FEmrt <- function(object, newdata, ...){
  if (!inherits(object, "FEmrt"))
    warning("calling predict.FEmrt(<fake-FEmrt-object>) ...")
  if (length(object$n) < 2) {
    warning("No tree was detected, all effect sizes are predicted as overall effect size")
    data.frame(g = rep(object$g, nrow(newdata)))
  } else {
    mf <- as.formula(object$call)
    tt <- terms(mf)
    ms <- model.frame(delete.response(tt), newdata)
    tree <- object$tree
    pred.efk <- predict(tree, newdata, type = "vector", ...)
    inx <- match(pred.efk, object$g)
    ci.lb <- object$ci.lb[inx]
    ci.ub <- object$ci.ub[inx]
    res <- data.frame(g = pred.efk, ci.lb = ci.lb, ci.ub = ci.ub, ms)
    row.names(res) <- NULL
    res
  }
  
}

#' Predictions from a fitted metacart object
#'
#' Returns a data frame of predicted effect sizes and moderators from a fitted metacart object
#'
#' @param object fitted model object of class "REmrt".
#' @param newdata data frame containing the values at which predictions are required.
#' @param ... Arguments that pass to other methods.
#' @return  A data frame containing the predicted effect size, the moderators, and the corresponding node lables in the fitted tree.
#' @importFrom stats as.formula delete.response model.frame model.response terms predict
#' @export
predict.REmrt <- function(object, newdata, ...){
  if (length(object$n) < 2) {
    warning("No tree was detected, all effect sizes are predicted as overall effect size")
    data.frame(g = rep(object$g, nrow(newdata)))
  } else {
    allNodes <- prednode_cpp(object, newdata)
    TNodes <- allNodes[, ncol(allNodes)] 
    pred.g <- object$g[as.character(TNodes)]
    ci.lb <- object$ci.lb[as.character(TNodes)]
    ci.ub <- object$ci.ub[as.character(TNodes)]
    res <- data.frame(g = pred.g, ci.lb = ci.lb, ci.ub = ci.ub, newdata)
    row.names(res) <- NULL
    res
    
  }

  }
