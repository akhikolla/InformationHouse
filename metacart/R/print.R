#' Print function for FEmrt
#'
#' Print the results of a FEmrt object
#'
#' @param x fitted tree of class \code{FEmrt}.
#' @param \dots additional arguments to be passed.
#' @details
#' The function returns the objects concerning the analysis results.
#' @export
print.FEmrt<- function(x, ...){
  if (length(x$n) < 2) {
    cat("\n")
    cat("Fixed Effects Meta-Tree (K = ", sum(x$n), " studies); ",
        sep = "")
    cat("\n")
    print(x$call)
    cat("\n")
    cat("No moderator effect was detected" )
    cat("\n")
    cat("use summary() to see the meta-analysis results")


  } else {
    cat("\n")
    cat("Fixed Effects Meta-tree (K = ", sum(x$n), " studies); ",
        sep = "")
    cat("\n")
    print(x$call)
    cat("\n")
    cat("A tree with ", length(x$n), " terminal nodes was detected", sep="" )
    cat("\n")
    cat("The moderators are ", paste(as.character(x$moderators), collapse = ", "), sep = "")
    cat("\n")
    cat("Use summary() and plot() to inspect the moderator analysis results and the tree structure.")
    cat("\n")
    print(x$tree)
  }
}


#' Print function for REmrt
#'
#' Print the results of a REmrt object
#'
#' @param x fitted tree of class \code{FEmrt}.
#' @param \dots additional arguments to be passed.
#' @details
#' The function returns the results (e.g., the value of the Q-between) after each split of the tree.
#' @export
print.REmrt<- function(x, ...){
  if (length(x$n) < 2) {
    cat("\n")
    cat("Random Effects Meta-Tree (K = ", sum(x$n), " studies); ",
        sep = "")
    cat("\n")
    print(x$call)
    cat("\n")
    cat("No moderator effect was detected" )
    cat("\n")
    cat("Use summary() to inspect the meta-analysis results")


  } else {
    cat("\n")
    cat("Random Effects Meta-tree (K = ", sum(x$n), " studies); ",
        sep = "")
    cat("\n")
    print(x$call)
    cat("\n")
    cat("A tree with ", length(x$n), " terminal nodes was detected", sep="" )
    cat("\n")
    cat("The moderators are ", paste(as.character(x$moderators), collapse = ", "), sep = "")
    cat("\n")
    cat("use summary() and plot() to see the moderator analysis results and the tree structure")
    cat("\n")
    print(x$tree)
  }
}
