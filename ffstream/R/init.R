#will initialise FFF, AFF, FFFCD, AFFCD, CUSUM AND EWMA

#functions used from utils.R:
#   isInBounds()
#   isAboveBound()



#' Initialisation of FFF
#'
#' This function makes it simple to initalise an FFF mean object.
#'
#'
#' @param lambda The value of the fixed forgetting factor. 
#'               Default is lambda=1.
#'
#' @examples
#' library(Rcpp)
#' fff1 <- initFFFMean()     # initialises with lambda=1
#'
#' fff2 <- initFFFMean(0.9)  # initialises with lambda=0.9
#'
#' @export 
initFFFMean <- function(lambda=1){
    checkFFFargs(lambda=lambda, functionName="initFFFMean")
    return( new (FFF, lambda) )
}



#' Initialisation of AFF mean
#'
#' This function makes it simple to initalise an AFF object.
#'
#'
#' @param eta The value of the step size in the gradient descent. 
#'            Default is eta=0.01.
#'
#' @examples
#' library(Rcpp)
#' aff1 <- initAFFMean()     # initialises with eta=0.01
#'
#' aff2 <- initAFFMean(0.1)  # initialises with eta=0.1
#'
#' @export 
initAFFMean <- function(eta=0.1){
    checkAFFargs(eta=eta, functionName="initAFFMean")
    return( new (AFF, eta) )
}



#' Initialisation of FFF mean change detector
#'
#' This function makes it simple to initalise an FFF object.
#'              
#' @param alpha The value of the significance level. 
#'              Default value is \code{0.01}, although it is
#'              recommended that the user set this parameter.
#'              
#' @param lambda The value of the fixed forgetting factor. 
#'               Default value is \code{lambda=1}.
#'              
#' @param BL The length of the burn-in region. Default value is \code{BL=50}.
#'              
#' @examples
#' library(Rcpp)
#' fffcd1 <- initFFFMeanCD()               # initialises with alpha=0.01
#'                                         
#' fffcd2 <- initFFFMeanCD(0.05, 0.95)     # initialises with alpha=0.05
#'                                         # and lambda=0.95 (and BL=50)
#' @export 
initFFFMeanCD <- function(alpha=0.01, lambda=1, BL=50){
    BL <- checkFFFMeanCDargs(alpha, lambda, BL, "initFFFMeanCD")

    #if pass all these tests, then can initialise:
    return( new(FFFChangeDetector, lambda, alpha, BL))
}



#' Initialisation of AFF change detector
#'
#' This function makes it simple to initalise an FFF object.
#'
#'
#' @param alpha The value of the significance level. 
#'              Default value is \code{0.01}, although it is
#'              recommended that the user set this parameter.
#'
#' @param eta The value of the step-size in the gradient descent. 
#'            Default is \code{eta=0.01}.
#'
#' @param BL The length of the burn-in region. Default value is \code{BL=50}. 
#'           Must be at least greater than or equal to \code{2}. No maximum.
#'           However, there is an exception: \code{BL=0} also works, but in
#'           this case the user needs to specify the \code{streamEstMean} and
#'           \code{streamEstSigma}; see 
#'
#' @examples
#' library(Rcpp)
#' affmeancd1 <- initAFFMeanCD()              # initialises with alpha=0.01, 
#'                                            # eta=0.01 and BL=50
#' 
#' affmeancd2 <- initAFFMeanCD(alpha=0.005, eta=0.1,  BL=100) 
#'
#' 
#' affmeancd3 <- initAFFMeanCD(alpha=0.005, eta=0.1,  BL=0)     # Example 3
#' affmeancd3$streamEstMean <- 0
#' affmeancd3$streamEstSigma <- 1
#'
#'
#' @export 
initAFFMeanCD <- function(alpha=0.01, eta=0.01, BL=50){
    BL <- checkAFFMeanCDargs(alpha, eta, BL, "initAFFMeanCD")
    #if pass all these tests, then can initialise:
    return( new(AFFChangeDetector, alpha, eta, BL))
}



#' Initialisation of CUSUM
#'
#' This function makes it simple to initalise a CUSUM object.
#'
#'
#' @param k One of the CUSUM control parameters. 
#'          Default value is \code{k=0.25}.
#'
#' @param h One of the CUSUM control parameters. 
#'          Default value is \code{h=8.00}.
#'
#' @param BL The burn-in length to be used with a CUSUM change detector. 
#'           Default value is \code{BL=50}.
#'
#' @examples
#' library(Rcpp)
#' c1 <- initCUSUMMeanCD()     # initialises with k=0.25, h=8.00, BL=50
#'
#' c2 <- initCUSUMMeanCD(k=0.5, h=4.00, BL=30)  
#' 
#' @export 
initCUSUMMeanCD <- function(k = 0.25, h = 8.00, BL=50){
    BL <- checkCUSUMMeanCDargs(k, h, BL, "initCUSUMMeanCD")

    #if pass all these tests, then can initialise:
    return( new(CusumChangeDetector, k, h, BL) )
}



#' Initialisation of EWMA
#'
#' This function makes it simple to initialise a EWMA object.
#'
#'
#' @param r One of the EWMA control parameters. 
#'          Default value is \code{r=0.20}.
#'
#' @param L One of the EWMA control parameters. 
#'          Default value is \code{L=3.00}.
#'
#' @param BL The burn-in length to be used with a EWMA change detector. 
#'           Default value is \code{BL=50}.
#'
# @examples
#' library(Rcpp)
#' e1 <- initEWMAMeanCD()                   #initialises with r=0.20, L=3.00
#'
#' e1 <- initEWMAMeanCD(r=0.05, L=0.275)    #initialises with r=0.20, L=3.00
#'
#' @export
initEWMAMeanCD <- function(r=0.2, L=3.0, BL=50){
    BL <- checkEWMAMeanCDargs(r, L, BL, "initEWMAMeanCD")

    #if pass all these tests, then can initialise:
    return( new(EwmaChangeDetector, r, L, BL) )
}
