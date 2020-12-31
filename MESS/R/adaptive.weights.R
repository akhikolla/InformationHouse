#' Compute weights for use with adaptive lasso.
#' 
#' Fast computation of weights needed for adaptive lasso based on Gaussian
#' family data.
#' 
#' The weights returned are 1/abs(beta_hat)^nu where the beta-parameters are
#' estimated from the corresponding linear model (either multivariate or
#' univariate).
#' 
#' @param x input matrix, of dimension nobs x nvars; each row is an observation
#' vector.
#' @param y response variable.
#' @param nu non-negative tuning parameter
#' @param weight.method Should the weights be computed for multivariate
#' regression model (only possible when the number of observations is larger
#' than the number of parameters) or by individual marginal/"univariate"
#' regression coefficients.
#' @return Returns a list with two elements: \item{weights }{the computed
#' weights} \item{nu }{the value of nu used for the computations}
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{glmnet}
#' @references Xou, H (2006). The Adaptive Lasso and Its Oracle Properties.
#' JASA, Vol. 101.
#' @keywords manip
#' @examples
#' 
#' set.seed(1)
#' x <- matrix(rnorm(50000), nrow=50)
#' y <- rnorm(50, mean=x[,1])
#' weights <- adaptive.weights(x, y)
#'
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'     res <- glmnet::glmnet(x, y, penalty.factor=weights$weights)
#'     head(res)
#' }
#' 
#' @export adaptive.weights
adaptive.weights <- function(x, y, nu=1, weight.method=c("multivariate", "univariate")) {

    weight.method <- match.arg(weight.method)

    nobs <- nrow(x)
    if (nobs+1 <= ncol(x)) {
        warning("using univariate weight method since p>n")
        weight.method <- "univariate"
    }

    # Get OLS fits
    weights <- switch(weight.method,
                      "univariate" = apply(x, 2, function(xi){lm.fit(x=cbind(1,xi), y=y)$coefficients[2]}),
                      "multivariate" =  (lm.fit(cbind(1, x), y)$coefficients)[-1]
                      )
    weights <- 1/abs(weights)^nu
    list(weights=weights, nu=nu)
}


