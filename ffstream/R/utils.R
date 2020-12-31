## Some utility functions in R

#-------------------------------------------------------------------------#
#' Check if an argument is in the closed set [lower, upper]
#'
#' Function to check that a value \code{x} is between two values.
#'
#'
#' @param x The value to check.
#'
#' @param lower The lower bound.
#'
#' @param upper The upper bound.
#'
#'
#' @keywords internal
isInBounds <- function(x, lower, upper){
    inBounds <- TRUE
    if (is.finite(x)){
        if ( (x < lower) || (x > upper) )
            inBounds <- FALSE
    } else {
        inBounds <- FALSE
    }
    return (inBounds)
}



#-------------------------------------------------------------------------#
#' Check if an argument is in the open set (lower, upper)
#'
#' Function to check that a value \code{x} is between strictly two values.
#'
#'
#' @param x The value to check.
#'
#' @param lower The lower bound.
#'
#' @param upper The upper bound.
#'
#'
#' @keywords internal
isInBoundsStrictly <- function(x, lower, upper){
    inBounds <- FALSE
    if (is.finite(x)){
        if ( (x > lower) && (x < upper) )
            inBounds <- TRUE
    } 
    return (inBounds)
}


#-------------------------------------------------------------------------#
#' Check if an argument above a certain value
#'
#' Function to check that a value \code{x} is above a certain lower bound.
#'
#'
#' @param x The value to check.
#'
#' @param lower The lower bound.
#'
#'
#' @keywords internal
isAboveBound <- function(x, lower){
    aboveBound <- TRUE
    if (is.finite(x)){
        if (x < lower) 
            aboveBound <- FALSE
    } else {
        aboveBound <- FALSE
    }
    return (aboveBound)
}



#-------------------------------------------------------------------------#
#' Check if a value is an integer
#'
#' Function to check that a value \code{x} is an integer.
#'
#'
#' @param x The value to check.
#'
#'
#' @keywords internal
isInteger <- function(x){
    isInteger <- FALSE
    if (is.finite(x)){
        if(x%%1==0)
            isInteger <- TRUE
    }
    return(isInteger)
}


#-------------------------------------------------------------------------#
#' Compute a one-sided p-value from a quantile value
#'
#' Wrapper function to compute a one-sided p-value from a statistics
#' or quantile value, with the aid of lower and upper bounds.
#'
#'
#' @param x The two-sided p-value to be converted.
#'
#' @param a The lower bound, e.g. 0.
#'
#' @param b The upper bound, e.g. 2.
#'
#'
#' @keywords internal
computeOneSidedPvalueR <- function(x, a, b){
    return(computeOneSidedPvalueRcpp(x, a, b))
}



#-------------------------------------------------------------------------#
#' Convert a 'p-value' to the correct side.
#'
#' Wrapper function to convery a two-sided p-value to a one-sided
#' 'p-value' (values close to 0 and close to 1 both become close to 0).
#' For example, both 0.01 and 0.99 are mapped to 0.01.
#'
#'
#' @param p The two-sided p-value to be converted.
#'
#'
#' @keywords internal
convertPvalueToCorrectSideR <- function(p){
    return(convertPvalueToCorrectSideRcpp(p))
}



#-------------------------------------------------------------------------#
#' Combine two p-values into a single p-value
#'
#' Wrapper function to combine two p-values into a single p-value.
#' Uses a type of 'reverse Bonferroni'.
#'
#' @param p1 The first p-value.
#'
#' @param p2 The second p-value.
#'
#'
#' @keywords internal
combineTwoOneSidedPvaluesR <- function(p1, p2){
    return(combineTwoOneSidedPvaluesRcpp(p1, p2))
}


#-------------------------------------------------------------------------#
#' Compute a two-sided p-value from a quantile value
#'
#' Wrapper function to compute a two-sided p-value from a statistics
#' or quantile value, with the aid of lower and upper bounds.
#'
#'
#' @param x The two-sided p-value to be converted.
#'
#' @param a The lower bound, e.g. 0.
#'
#' @param b The upper bound, e.g. 2.
#'
#'
#' @keywords internal
computeTwoSidedPvalueR <- function(x, a, b){
    return(computeTwoSidedPvalueRcpp(x, a, b))
}



## a simple function to create a stream with multiple changepoints
#numChanges is the number of changes in the stream
#regimeLength is the number of observations per regime (also gives location of changepoint)
#seednum is an optional argument in order to set the seed of the random number generator
#-------------------------------------------------------------------------#
#' Create a stream with a certain number of changes
#'
#' A function, used almost exclusively in testing, to generate a stream
#' of observations with a certain number of changepoints.
#'
#'
#' @param x numChanges The number of changepoints. Default is \code{3}.
#'
#' @param regimeLength The number of observations between changepoints
#'                     (before and after). Default is \code{150}.
#'
#' @param seednum The seed for the random number generator. Default is
#'                \code{NULL}, in which case no seed will be used.
#'
#' @param mu0 The initial mean of the process.
#'
#' @param sigma0 The variance of each stream.
#'
#'
#' @details Creates a stream of observations, where after each changepoint
#'          the mean increases by one standard deviation.
#'
#'
#' @return A vector of observations.
#'
#'
#' @keywords internal
makeStreamMeanChangeR <- function(numChanges=3, regimeLength=150, seednum=NULL, mu0=0, sigma0=1){
    if (is.numeric(seednum)){
        set.seed(seednum)
    }
    #need numChanges plus one, because need to "bookend" the changepoints
    N <- (regimeLength) * (numChanges + 1)
    stepChanges <- rep(c(0:numChanges), each=regimeLength)
    x <- rnorm(N, mean=mu0, sd=sigma0)  +  stepChanges
    return(x)
}
