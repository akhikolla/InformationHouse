#' Screen variable before penalized regression
#'
#' Expands a contingency table to a data frame where each observation in the table becomes a single observation in the data frame with corresponding information for each for each combination of the table dimensions.
#'
#' Note that no standardization is done (not necessary?)
#'
#' @param x A table or matrix
#' @param y A vector of outcomes
#' @param lambda a vector of positive values used for the penalization parameter.
#' @param method a string giving the method used for screening. Two possibilities are "global-strong" and "global-DPP"
#' @references Hastie, Tibshirani and Wainwright (2015). "Statistical Learning with Sparsity". CRC Press.
#' @return A list with three elements: lambda which contains the lambda values, selected which contains the indices of the selected variables, and method a string listing the method used.
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @keywords manip
#' @examples
#'
#' x <- matrix(rnorm(50*100), nrow=50)
#' y <- rnorm(50, mean=x[,1])
#' screen_variables(x, y, lambda=c(.1, 1, 2))
#'
#' @export
screen_variables <- function(x, y, lambda=.1, method=c("global-strong", "global-DPP")) {

    ## Sanity checks
    if (!any(c("matrix", "table") %in% class(x)))
        stop("needs matrix or table as input")

    if (any(lambda<=0))
        stop("lambda must be positive")

    method <- match.arg(method)

    ## Compute  |x^t y|
    res <- abs(as.vector(crossprod(x, y)))

    lambdamax <- max(res)

    ## Criterion
    if (method=="global-strong") {
        criterion <- 2*lambda-lambdamax
    } else {
        criterion <- lambdamax - norm(x, 2)*sqrt(sum(y^2))*(lambdamax - lambda)/lambda
    }

    discard <- outer(res, criterion, "<")

    # The code below might be modified and fast but requires data.table

    list(lambda=lambda,
         selected=unname(apply(discard, 2, which)),
         method=method)
}
