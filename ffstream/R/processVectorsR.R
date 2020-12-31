#' Quick computation of FFF mean of a given vector
#'
#' Given a vector \code{x} and a value \code{lambda} for a fixed forgetting
#' factor, returns the value of the fixed forgetting factor mean
#' \eqn{\bar{x}_{N, \lambda}}, where \eqn{N} is the length of \code{x}.
#' Algorithm is implemented in C++.
#'
#' @param x Vector of numeric values values. Default is \code{c(0)}, 
#'          a vector of one element (zero)
#'
#' @param lambda Value for the fixed forgetting factor in \eqn{[0,1]}.
#'               Default is \code{lambda=0.99}.
#'
#' @seealso \code{\link{computeAFFMean}}
#'
#'
#' @section Author:
#' Dean Bodenham
#'
#'
#' @section References:
#' D. A. Bodenham and N. M. Adams (2016) 
#' \emph{Continuous monitoring for changepoints in data 
#' streams using adaptive estimation}. 
#' Statistics and Computing  
#' doi:10.1007/s11222-016-9684-8
#'
#'
#' @export
computeFFFMean <- function(x=c(0), lambda=0.99){

    LAMBDAMIN <- 0
    LAMBDAMAX <- 1

    #some checks
    #check none are NA
    #(vectors cannot contain NULL values, so no need to check)
    if(anyNA(x)){
        stop("x contains NA values.")
    }

    #check all numeric
    if (is.numeric(x)){
        #great!
        #now check if any are infinite
        if ( any(is.infinite(x)) ){
            #could give index?
            stop("x contains non-finite values.")
        }
    } else {
        stop("x contains non-numeric values.")
    }

    if (is.null(lambda)){
        stop("lambda is NULL.")
    }
    if (is.finite(lambda)){
        #great!
    } else {
        stop("lambda is not a finite numeric value.")
    }

    if ( !(isInBounds(lambda, LAMBDAMIN, LAMBDAMAX)) ){ 
        message <- paste0("lambda not in range [", LAMBDAMIN, ", ", LAMBDAMAX, "]")
        stop(message)
    }

    xbar <- cpp_computeFFFMean(x, lambda)
    return(xbar)
}


#' Quick computation of AFF mean of a given vector
#'
#' Given a vector \code{x} and a value \code{eta} for step-size
#' in the stochastic gradient descent for the adaptive forgetting
#' factor, this returns the value of the fixed forgetting factor mean
#' \eqn{\bar{x}_{N, \overrightarrow{\lambda} }}, where \eqn{N} is the 
#' length of \code{x}. Algorithm is implemented in C++.
#'
#' @param x Vector of numeric values values. Default is \code{c(0)}, 
#'          a vector of one element (zero)
#'
#' @param eta Value for the step size in the gradient descent step. 
#'        Default is \code{eta=0.01}.
#'
#' @seealso \code{\link{computeFFFMean}}
#'
#'
#' @section Author:
#' Dean Bodenham
#'
#'
#' @section References:
#' D. A. Bodenham and N. M. Adams (2016) 
#' \emph{Continuous monitoring for changepoints in data 
#' streams using adaptive estimation}. 
#' Statistics and Computing  
#' doi:10.1007/s11222-016-9684-8
#'
#'
#' @export
computeAFFMean <- function(x=c(0), eta=0.01){

    ETAMIN <- 0
    ETAMAX <- 1
    #some checks
    #check none are NA
    if(anyNA(x)){
        stop("x contains NA values.")
    }

    #check all numeric
    if (is.numeric(x)){
        #great!
        #now check if any are infinite
        if ( any(is.infinite(x)) ){
            #could give index?
            stop("x contains non-finite values.")
        }
    } else {
        stop("x contains non-numeric values.")
    }

    if (is.null(eta)){
        stop("eta is NULL.")
    }
    if (is.finite(eta)){
        #great!
    } else {
        stop("eta is not a finite numeric value.")
    }

    if ( !(isInBounds(eta, ETAMIN, ETAMAX)) ){ 
        message <- paste0("eta not in range [", ETAMIN, ", ", ETAMAX, "]")
        stop(message)
    }

    xbar <- cpp_computeAFFMean(x, eta)
    return(xbar)
}


