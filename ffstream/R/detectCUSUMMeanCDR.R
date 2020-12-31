#' Detect a change/changes in a vector using CUSUM method
#'
#' Given a vector \code{x}, use the CUSUM method to sequentially detect 
#' changes (or a single change) in the MEAN of the vector. 
#'
#'
#' @param x The vector (stream) in which to detect change(s).
#'
#' @param k control parameter for CUSUM. Default is \code{0.25}.
#'
#' @param h control parameter for CUSUM. Defqult is \code{8.00}.
#'
#' @param BL The burn-in length. Default is \code{BL=50}.
#'
#' @param multiple Boolean to use to decide whether to detect multiple changes
#'                or only a single change. Default is \code{TRUE} (i.e. detect 
#'                multiple changes).
#'
#' @param single Boolean to use to decide whether to detect only a single 
#'               change or multiple changes. Set to \code{!multiple}, i.e.
#'               default is FALSE. If both \code{single} and \code{multiple}
#'               are set to \code{TRUE}, then only a single change will be
#'               detected; if both set to \code{FALSE} then multiple changes
#'               will be detected (i.e. \code{single} dominates).
#'
#' @param usePrechange Boolean indicating whether prechange parameters
#'                       (mean and variance) are known and will be used 
#'                       (or not). Default is
#'                       \code{FALSE}. If \code{TRUE}, then prechange mean
#'                       and standard deviation variance must be specified.
#'                       See parameters \code{prechangeMean},
#'                       \code{prechangeSigma} and \code{prechangeVar}.
#'
#' @param prechangeMean Value to be used for the prechange mean. 
#'                      Default is \code{NULL}. 
#'                      If \code{prechangeKnown = TRUE} and value is
#'                      \code{NULL}, this will result in an error.
#'
#' @param prechangeSigma Value to be used for the prechange standard 
#'                       deviation. Default is \code{NULL}. 
#'                       If \code{prechangeKnown = TRUE} and value is
#'                       \code{NULL}, this will result in an error, unless
#'                       \code{prechangeVar} is not \code{NULL}.
#'
#' @param prechangeVar Value to be used for the prechange variance. 
#'                     Default is \code{NULL}. 
#'                     If \code{prechangeKnown = TRUE} and value is
#'                     \code{NULL}, this will result in an error, unless
#'                     \code{prechangeSigma} is not \code{NULL}.
#'                     \code{prechangeVar} is set to 
#'                     \code{sqrt(prechangeSigma)}.
#'
#' @param skipCheck A boolean which allows the function to skip the check 
#'                   of the stream. Default is \code{FALSE}.
#'
#' 
#' @details CUSUM updates via: 
#'          \deqn{S_{j} = \max{0, S_{j-1} + (x_{j} - \mu)/ \sigma - k}}
#'          and 
#'          \deqn{T_{j} = \max{0, S_{j-1} - (x_{j} - \mu)/ \sigma - k}}
#'          where \eqn{\mu} and \eqn{\sigma} are, respectively, the mean 
#'          and variance of the in-control stream, 
#'          \eqn{x_j} is the observation at time \eqn{j}
#'          and \eqn{k} 
#'          is a control parameter for CUSUM. Then, a change is signalled
#'          if \eqn{S_j > h} or \eqn{T_j > h},
#'          where \eqn{h} is the other control parameter. This is the 
#'          formulation for using CUSUM to detect an increase or decrease
#'          in the mean.
#'
#'
#' @return A list with the following elements:
#'         \describe{
#'             \item{\code{tauhat}}{A vector of the changepoints found.}
#'          }
#'
#'
#' @section Author:
#' Dean Bodenham
#'
#' @section References:
#' E. S. Page (1954) \emph{Continuous inspection schemes}. 
#' Biometrika, 41(1/2), 100-115
#'
#'
#' @examples
#' # create a stream with three changepoints
#' set.seed(8)
#' x <- rnorm(400, 5, 1) + rep(c(0:3), each=100) # mean is 5 and s.d. is 1
#'
#' # multiple changepoints
#' list_cusum <- detectCUSUMMean(x, k=0.25, h=8.00, BL=50, multiple=TRUE)
#'
#' # now only a single (the first) changepoint
#' list_cusum2 <- detectCUSUMMean(x, k=0.25, h=8.00, BL=50, single=TRUE)
#'
#' # now only a single (the first) changepoint, but with the prechange 
#' # mean and variance known
#' list_cusum3 <- detectCUSUMMean(x, k=0.25, h=8.00, BL=50, single=TRUE,
#'                                prechangeMean=5, prechangeSigma=1)
#'
#'
#'
#'
#' @export
detectCUSUMMean <- function(x, k=0.25, h=8.00, BL=50, multiple=TRUE, 
                        single=!multiple, usePrechange=FALSE, prechangeMean=NULL,
                        prechangeSigma=NULL, prechangeVar=NULL, skipCheck=FALSE){

    
    #check FFF parameters are properly set
    checkCUSUMMeanCDargs(k, h, BL, functionName="detectCUSUMMean")

    #check booleans
    checkBooleans(multiple, single, usePrechange, skipCheck, "detectCUSUMMean")

    if (skipCheck==FALSE){
        #check stream
        checkStream(x, "detectCUSUMMean")
    }

    if (single){
        if (usePrechange){
            #first check prechange
            checkPrechange(prechangeMean, prechangeSigma, prechangeVar, "detectCUSUMMean")
            return ( cpp_detectCUSUMMeanSinglePrechange(x, k, h, 
                                                   prechangeMean, prechangeSigma) )
        } else {
            return ( cpp_detectCUSUMMeanSingle(x, k, h, BL) )
        }
    }
    #last option...
    return ( cpp_detectCUSUMMeanMultiple(x, k, h, BL) );
}

