#' Detect a change/changes in a vector using FFF method
#'
#' Given a vector \code{x}, use the FFF method to sequentially detect changes 
#' (or a single change) in the MEAN of the vector. 
#'
#'
#' @param x The vector (stream) in which to detect change(s).
#'
#' @param lambda The value for the forgetting factor. 
#'              Default is \code{lambda=0.95}.
#'
#' @param alpha The value for the threshold. Default is \code{alpha=0.01}.
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
#' @return A list with the following elements:
#'         \describe{
#'             \item{\code{tauhat}}{A vector of the changepoints found.}
#'          }
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
#' @examples
#' # create a stream with three changepoints
#' set.seed(8)
#' x <- rnorm(400, 5, 1) + rep(c(0:3), each=100) # mean is 5 and s.d. is 1
#'
#' # multiple changepoints
#' list_fff <- detectFFFMean(x, alpha=0.01, lambda=0.95, BL=50, multiple=TRUE)
#'
#' # now only a single (the first) changepoint
#' list_fff2 <- detectFFFMean(x, alpha=0.01, lambda=0.95, BL=50, single=TRUE)
#'
#' # now only a single (the first) changepoint, but with the prechange 
#' # mean and variance known
#' list_fff3 <- detectFFFMean(x, alpha=0.01, lambda=0.95, BL=50, single=TRUE,
#'                           prechangeMean=5, prechangeSigma=1)
#'
#'
#' @export
detectFFFMean <- function(x, lambda=0.95, alpha=0.01, BL=50, multiple=TRUE, 
                        single=!multiple, usePrechange=FALSE, prechangeMean=NULL,
                        prechangeSigma=NULL, prechangeVar=NULL, skipCheck=FALSE){


    #check FFF parameters are properly set
    checkFFFMeanCDargs(alpha, lambda, BL, functionName="detectFFFMean")

    #check booleans
    checkBooleans(multiple, single, usePrechange, skipCheck, "detectFFFMean")

    if (skipCheck==FALSE){
        #check stream
        checkStream(x, "detectFFFMean")
    }

    if (single){
        if (usePrechange){
            #first check prechange
            checkPrechange(prechangeMean, prechangeSigma, prechangeVar, "detectFFFMean")
            return ( cpp_detectFFFMeanSinglePrechange(x, lambda, alpha, 
                                                   prechangeMean, prechangeSigma) )
        } else {
            return ( cpp_detectFFFMeanSingle(x, lambda, alpha, BL) )
        }
    }
    #last option - do multiple
    return ( cpp_detectFFFMeanMultiple(x, lambda, alpha, BL) );
}
