#' Calculate difference between areas under the ROC curve (AUC) between a forecast and a reference forecast for the same observation, and estimate the variance of the AUC difference
#'
#' @param fcst vector of forecasts
#' @param fcst.ref vector of reference forecasts
#' @param obs vector of binary observations (0 for non-occurrence, 1 for occurrence of the event)
#' @param handle.na how should missing values in forecasts and observations be handled; possible values are 'na.fail' and 'only.complete.triplets'; default: 'na.fail'
#' @param use_fn the function used for the calculation: 'C++' (default) for the fast C++ implementation, or 'R' for the slow (but more readable) R implementation 
#' @return vector with AUC difference, and estimated standard deviation
#' @examples
#' data(eurotempforecast)
#' AucDiff(rowMeans(ens.bin), ens.bin[, 1], obs.bin)
#' @seealso Auc
#' @references DeLong et al (1988): Comparing the Areas under Two or More Correlated Receiver Operating Characteristic Curves: A Nonparametric Approach. Biometrics.  \url{https://www.jstor.org/stable/2531595}
#' Sun and Xu (2014): Fast Implementation of DeLong's Algorithm for Comparing the Areas Under Correlated Receiver Operating Characteristic Curves. IEEE Sign Proc Let 21(11). \doi{10.1109/LSP.2014.2337313}
#' @export
AucDiff <- function(fcst, fcst.ref, obs, handle.na=c('na.fail', 'only.complete.triplets'), 
                    use_fn=c('C++','R')) {


  ## sanity checks
  stopifnot(length(fcst) == length(obs))
  stopifnot(length(fcst.ref) == length(obs))


  ## handle NA's
  handle.na = match.arg(handle.na)
  if (handle.na == "na.fail") {
    if (any(is.na(c(fcst, fcst.ref, obs)))) {
      stop("missing values")
    }
  } 
  if (handle.na == "only.complete.triplets") {
    nna <- !is.na(fcst) & !is.na(fcst.ref) & !is.na(obs)
    if (all(nna == FALSE)) {
      stop("there are no complete triplets of forecasts and observations")
    }
    fcst <- fcst[nna]
    fcst.ref <- fcst.ref[nna]
    obs <- obs[nna]
  }


  ## after removing any NA's, check if observations are either 0 or 1
  stopifnot(all(obs %in% c(0,1)))
  L = length(obs)
  m = sum(obs)
  n = L - m
  if (m == 0 | n == 0) {
    stop("need at least one event and one non-event")
  }

  
  # decide whether to use the fast C function or the slow but more readable R function 

  use_fn = match.arg(use_fn)

  if (use_fn == 'C++') {

    ans = aucdiff_cpp(fcst, fcst.ref, obs)

  }

  if (use_fn == 'R') {

  
    ## calculate sets of forecasts with events (X) and forecasts with non-events (Y)
    X <- cbind(fcst[obs == 1], fcst.ref[obs == 1])
    Y <- cbind(fcst[obs == 0], fcst.ref[obs == 0])
  
  
    ## Delong's Psi function as matrix (Psi.mat[i, j] = Psi(X[i], Y[j]))
    psi.fun <- function(x, y) {
      return(1 * (x > y) + 0.5 * (x == y))
    }
    Psi <- outer(X[, 1], Y[, 1], psi.fun)
    Psi.ref <- outer(X[, 2], Y[, 2], psi.fun)
  
    
    ## AUC estimates
    auc <- mean(Psi)
    auc.ref <- mean(Psi.ref)
  
  
    ## AUC difference
    auc.diff <- auc - auc.ref
  
  
    ## Delong (1988) variance estimation
    V <- rowMeans(Psi)
    V.ref <- rowMeans(Psi.ref)
    W <- colMeans(Psi)
    W.ref <- colMeans(Psi.ref)
  
    v.11 <- var(V)
    v.22 <- var(V.ref)
    v.12 <- cov(V, V.ref)
  
    w.11 <- var(W)
    w.22 <- var(W.ref)
    w.12 <- cov(W, W.ref)

    sd.auc = sqrt(v.11 / m + w.11 / n)
    sd.auc.ref = sqrt(v.22 / m + w.22 / n)
    sd.auc.diff = sqrt((v.11 + v.22 - 2 * v.12) / m +
                       (w.11 + w.22 - 2 * w.12) / n)

    ans = c(auc, sd.auc, auc.ref, sd.auc.ref, auc.diff, sd.auc.diff)
  
  }


  ## return
  names(ans) = c('auc', 'auc_sd', 'auc_ref', 'auc_ref_sd', 'auc_diff', 'auc_diff_sd')
  return(ans)

}

