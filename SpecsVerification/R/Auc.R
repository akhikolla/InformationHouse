#' Calculate area under the ROC curve (AUC) for a forecast and its verifying binary observation, and estimate the variance of the AUC 
#'
#' @param fcst vector of forecasts
#' @param obs vector of binary observations (0 for non-occurrence, 1 for occurrence of the event)
#' @param handle.na how should missing values in forecasts and observations be handled; possible values are 'na.fail' and 'only.complete.pairs'; default: 'na.fail'
#' @param use_fn the function used for the calculation: 'C++' (default) for the fast C++ implementation, or 'R' for the slow (but more readable) R implementation 
#' @return vector containing AUC and its estimated sampling standard deviation
#' @examples
#' data(eurotempforecast)
#' Auc(rowMeans(ens.bin), obs.bin)
#' @seealso AucDiff
#' @references DeLong et al (1988): Comparing the Areas under Two or More Correlated Receiver Operating Characteristic Curves: A Nonparametric Approach. Biometrics. \url{https://www.jstor.org/stable/2531595}
#' Sun and Xu (2014): Fast Implementation of DeLong's Algorithm for Comparing the Areas Under Correlated Receiver Operating Characteristic Curves. IEEE Sign Proc Let 21(11). \doi{10.1109/LSP.2014.2337313}
#' @export
Auc = function(fcst, obs, handle.na=c("na.fail", "only.complete.pairs"), use_fn=c('C++','R')) {


  ## sanity checks
  stopifnot(length(fcst) == length(obs))


  ## handle NA's
  handle.na = match.arg(handle.na)
  if (handle.na == "na.fail") {
    if (any(is.na(c(fcst, obs)))) {
      stop("missing values")
    }
  } else if (handle.na == "only.complete.pairs") {
    nna = !is.na(fcst) & !is.na(obs)
    if (all(nna == FALSE)) {
      stop("there are no complete sets of forecasts and observations")
    }
    fcst = fcst[nna]
    obs = obs[nna]
  } 


  ## after removing any NA's, check if observations are either 0 or 1
  stopifnot(all(obs %in% c(0,1)))
  n = length(obs)
  n1 = sum(obs)
  n0 = n - n1 
  if (n0 == 0 | n1 == 0) {
    stop("need at least one event and one non-event")
  }

  
  ## decide whether to use the fast C function or the slow but more readable R function

  use_fn = match.arg(use_fn)

  if (use_fn == 'C++') {

    ans = auc_cpp(fcst, obs)

  } 


  if (use_fn == 'R') {

    ## calculate sets of forecasts with events (X) and forecasts with non-events (Y)
    X = fcst[obs == 1]
    Y = fcst[obs == 0]
  
  
    ## Delong's Psi function as matrix (Psi_mat[i, j] = Psi(X[i], Y[j]))
    psi_fun = function(x, y) {
      return(1 * (x > y) + 0.5 * (x == y))
    }
    Psi = outer(X, Y, psi_fun)
  
    
    ## AUC estimates
    auc = mean(Psi)
  
  
    ## Delong et al (1988) variance estimate
    V = rowMeans(Psi)
    W = colMeans(Psi)
    v = var(V)
    w = var(W)
    auc_sd = sqrt(v / n1 + w / n0)
  
    ans = c(auc, auc_sd)

  }


  ## return
  names(ans) = c('auc', 'auc_sd')
  return(ans)

}

