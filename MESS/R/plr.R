#' Fast computation of several simple linear regressions
#'
#' Fast computation of several simple linear regression, where the outcome is analyzed with several 
#' marginal analyses, or where several outcome are analyzed separately, or a combination of both.
#'
#' @param y either a vector (of length N) or a matrix (with N rows)
#' @param x a matrix with N rows
#' @param addintercept boolean. Should the intercept be included in the model by default (TRUE)
#' @return a data frame (if Y is a vector) or list of data frames (if Y is a matrix)
#' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk}
#' @seealso \code{mfastLmCpp}
#' @keywords datagen
#' @examples
#' 
#' N <- 1000  # Number of observations
#' Nx <- 20   # Number of independent variables
#' Ny <- 80   # Number of dependent variables
#'
#' # Simulate outcomes that are all standard Gaussians
#' Y <- matrix(rnorm(N*Ny), ncol=Ny)  
#' X <- matrix(rnorm(N*Nx), ncol=Nx)
#'
#' plr(Y, X)
#'
#' @rdname plr
#' @export
plr <- function (y, x, addintercept=TRUE)
{
  UseMethod("plr")
}


#' @rdname plr
#' @export
plr.numeric <- function (y, x, addintercept=TRUE)
{
    mfastLmCpp(y=y, x=x, addintercept=addintercept)
}

#' @rdname plr
#' @export
plr.matrix <- function (y, x, addintercept=TRUE)
{
    apply(y, 2, FUN=function(Y) { mfastLmCpp(y=Y, x=x, addintercept=addintercept) })
}

