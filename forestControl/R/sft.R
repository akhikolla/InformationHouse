#' Selection Frequency Threshold
#'
#' Determine the selecton frequency threshold of a model at a specified false positive rate
#'
#' @param x a `randomForest` or `ranger` object
#' @param alpha a false positive rate (ie, 0.01)
#' @return a list of two elements
#' * __sft__ Tthe selection frequency threshold
#' * __probs_atsft__ The esimated false positive rate
#'
#' @author Tom Wilson \email{tpw2@aber.ac.uk}
#' @export
#' @examples
#' library(randomForest)
#' data(iris)
#' iris.rf <- randomForest(iris[,-5], iris[,5], forest = TRUE)
#'
#' # For a false positive rate of 1%
#' iris.sft <- sft(iris.rf, 0.01)
#' print(iris.sft)
#'
#' # To iterate through a range of alpha values
#'
#' alpha <- c(0.01,0.05, 0.1,0.15,0.2, 0.25)
#' threshold <- NULL
#' for(i in seq_along(alpha)){
#'     threshold[i] <- sft(iris.rf, alpha[i])$sft
#' }
#'
#' plot(alpha, threshold, type = 'b')
#'

sft <- function(x, alpha)
  {
  params <- extract_params(x)
  sft_res <- sft_calc(Ft = params$Ft,Fn = params$Fn, K = params$K, Tr = params$Tr, alpha = alpha)
  sft_res <- list(sft = round(sft_res[1]), prob = sft_res[2])
  return(sft_res)
  }
