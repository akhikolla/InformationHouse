#' Prune a tree
#'
#' Prune an initial rpart tree by "c-standard-error" rule.
#' @inheritParams prune.rpart
#' @param tree A initial tree fitted by rpart, needs to an rpart object.
#' @param c A scalar to prune the  tree by selecting the tree with minum cross-validation error plus the standard error multiplied by c.
#' @param ... Additional arguments passed to prune.rpart().
#' @return The pruned tree
#' @importFrom rpart prune.rpart
#' @importFrom methods is
#' @keywords internal
treepruner <- function(tree, c, ...){
  prune <- NULL
  if (!is(tree, "rpart")) stop("The pruned tree should be an rpart object")
  else {
    if (!is.numeric(c) | length(c) != 1 | c < 0) {
      stop("The pruning parameter c should be a positive constant number")
    } else {
      tree <- tree
      c <- c
      mindex <- which.min(tree$cptable[,4])  # find the row of the minimum x-error
      cp.minse <- tree$cptable[mindex,4] + c*tree$cptable[mindex,5]  # the minimum x-error + c*SE
      cp.row <- min(which(tree$cptable[,4]<= cp.minse))  # find the smallest tree within the minimum x-error + c*SE
      cp.take <- tree$cptable[cp.row, 1]  # get the cp value for the smallest tree
      prune(tree, cp=cp.take, ...)  # prune the tree
    }
  }
}

