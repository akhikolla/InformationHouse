#' Detect a change/changes in a vector using EWMA method
#'
#' Given a vector \code{x}, use the EWMA method to sequentially detect changes 
#' (or a single change) in the MEAN of the vector. 
#'
#'
#' @param x The vector (stream) in which to detect change(s).
#'
#' @param r Control parameter for EWMA. Must be in range \eqn{[0,1]}.
#'          Default is \code{r=0.25}.
#'
#' @param L Control parameter for EWMA. Default is \code{L=3.00}
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
#' @details EWMA updates via: 
#'          \deqn{Z_{j} = (1-r) Z_{j-1} + r x_{j}}
#'          where \eqn{\mu} is the mean of the in-control stream, 
#'          \eqn{x_j} is the observation at time \eqn{j} and \eqn{r} 
#'          is a control parameter for EWMA. Then, a change is signalled
#'          if \deqn{|Z_j - \mu|  > L \sigma_{Z_j}}, 
#'          where \eqn{L} is the other control parameter, and 
#'          \eqn{\sigma_{Z_j}} is a scaled version of the in-control
#'          variance \eqn{\sigma}.
#'          This is the formulation for using EWMA to detect an increase or
#'          decrease in the mean.
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
#' S. W. Roberts (1959) \emph{Control chart tests based 
#' on geometric moving averages}. Technometrics, 1(3), 239-250
#'
#'
#' @examples
#' # create a stream with three changepoints
#' set.seed(8)
#' x <- rnorm(400, 5, 1) + rep(c(0:3), each=100) # mean is 5 and s.d. is 1
#'
#' # multiple changepoints
#' list_ewma <- detectEWMAMean(x, r=0.25, L=3.023, BL=50, multiple=TRUE)
#'
#' # now only a single (the first) changepoint
#' list_ewma2 <- detectEWMAMean(x, r=0.25, L=3.023, BL=50, single=TRUE)
#'
#' # now only a single (the first) changepoint, but with the prechange 
#' # mean and variance known
#' list_ewma3 <- detectEWMAMean(x, r=0.25, L=3.023, BL=50, single=TRUE,
#'                              prechangeMean=5, prechangeSigma=1)
#'
#'
#' @export
detectEWMAMean <- function(x, r=0.25, L=3.00, BL=50, multiple=TRUE, 
                        single=!multiple, usePrechange=FALSE, prechangeMean=NULL,
                        prechangeSigma=NULL, prechangeVar=NULL, skipCheck=FALSE){


    #check FFF parameters are properly set
    checkEWMAMeanCDargs(r, L, BL, functionName="detectEWMAMean")

    #check booleans
    checkBooleans(multiple, single, usePrechange, skipCheck, "detectEWMAMean")

    if (skipCheck==FALSE){
        #check stream
        checkStream(x, "detectEWMAMean")
    }
    if (single){
        if (usePrechange){
            #first check prechange
            checkPrechange(prechangeMean, prechangeSigma, prechangeVar, "detectEWMAMean")
            return ( cpp_detectEWMAMeanSinglePrechange(x, r, L, 
                                                   prechangeMean, prechangeSigma) )
        } else {
            return ( cpp_detectEWMAMeanSingle(x, r, L, BL) )
        }
    }
    #last option...
    return ( cpp_detectEWMAMeanMultiple(x, r, L, BL) );
}

