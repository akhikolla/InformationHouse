#' Rank histogram for ensemble forecasts
#'
#' @description Calculate the rank histogram for an archive of ensemble forecasts and their corresponding verifying observations.
#' @param ens matrix of dimension (N,K). An archive of K-member ensemble forecasts for N time instances.
#' @param obs vector of length N. The corresponding verifying observations.
#' @param reduce.bins number of adjacent bins that will be merged into one bin; has to be a divisor of K+1
#' @param handle.na how should missing values in ensemble and observation data be handled; possible values are 'na.fail' (fails if any data is missing) and 'use.complete' (only uses times where all ensemble members and obs are available); default: 'na.fail'
#' @return a vector of length (K+1)/reduce.bins containing the rank counts
#' @examples 
#' data(eurotempforecast)
#' rh <- Rankhist(ens, obs)
#' @references Anderson J.L. (1996). A Method for Producing and Evaluating Probabilistic Forecasts from Ensemble Model Integrations. J. Climate, 9, 1518--1530. 
#' Hammill T.M. (2001). Interpretation of Rank Histograms for Verifying Ensemble Forecasts. Mon. Wea. Rev., 129, 550--560. 
#' @seealso PlotRankhist, TestRankhist
#' @export

Rankhist <- function(ens, obs, reduce.bins=1, handle.na="na.fail") {

# @references Anderson J.L. (1996). A Method for Producing and Evaluating Probabilistic Forecasts from Ensemble Model Integrations. J. Climate, 9, 1518--1530. \doi{10.1175/1520-0442(1996)009%3C1518:AMFPAE%3E2.0.CO;2}
# Hammill T.M. (2001). Interpretation of Rank Histograms for Verifying Ensemble Forecasts. Mon. Wea. Rev., 129, 550--560. \doi{10.1175/1520-0493(2001)129%3C0550:IORHFV%3E2.0.CO;2}


  ## sanity checks
  stopifnot(nrow(ens) == length(obs))
  stopifnot(any(!is.na(ens)))
  stopifnot(any(!is.na(obs)))

  ## handle NA's

  # remove missing ensemble members (columns with all NA)
  nna.cols <- apply(ens, 2, function(z) !all(is.na(z)))
  ens <- ens[, nna.cols, drop=FALSE]

  # check handle.na argument
  if (handle.na == "na.fail") {
    if (any(is.na(c(ens, obs)))) {
      stop("missing values")
    }
  } else if (handle.na == "use.complete") {
    nna <- !is.na(obs) & !is.na(rowSums(ens))
    if (all(nna == FALSE)) {
      stop("there are missing values at all times")
    }
    obs <- obs[nna]
    ens <- ens[nna, ]
  } else {
    stop("unknown 'handle.na' argument")
  }


  N <- dim(ens)[1]
  K <- dim(ens)[2]

  if ((K+1) %% reduce.bins != 0) {
    stop("number of histogram bins is not a multiple of reduce.bins")
  }

  ranks <- apply(cbind(obs, ens), 1, rank, ties.method="random")[1, ]
  rank.hist <- hist(ranks, breaks=seq(0.5, K+1.5, 
                    reduce.bins), plot=FALSE)[["counts"]]
  return(rank.hist)
}

