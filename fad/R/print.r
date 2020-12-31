#' Print the Output of Factor Analysis
#' @rdname print.fad
#' @description Prints the output of the \code{fad}.
#' @param x an object of class \code{fad}.
#' @param digits number of decimal places to use in printing uniquenesses and loadings.
#' @param \dots further arguments to \code{print}.
#' @return None.
#' @export
print.fad <- function(x, digits = 3, ...)
{
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  cat("Uniquenesses:\n")
  print(round(x$uniquenesses, digits), ...)
  print(x$loadings, digits = digits, ...)
  if (!is.null(x$rotmat)){
    
    tmat <- solve(x$rotmat)
    R <- tmat %*% t(tmat)
    factors <- x$factors
    rownames(R) <- colnames(R) <- paste0("Factor", 1:factors)
    
    if (TRUE != all.equal(c(R), c(diag(factors)))){
      cat("\nFactor Correlations:\n")
      print(R, digits=digits, ...)
    }
    
    
  }
  
  if(!is.na(x$BIC))    cat("\nThe BIC is: ",x$BIC)
  
  invisible(x)
}
