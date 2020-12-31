#' The normal reference bandwidth selection for weighted data
#'
#' @description This function computes the data-driven bandwidth for smoothing the ROC (or distribution) function using the NR method of Beyene and El Ghouch (2020). This is an extension of the classical (unweighted) normal reference bandwith selection method to the case of weighted data.
#' 
#' @param X The numeric data vector.
#' @param wt The non-negative weight vector.
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". By default, the "\code{normal}" kernel is used.
#' @return Returns the computed value for the bandwith parameter.
#' @details See Beyene and El Ghouch (2020) for details.
#' @author
#' Kassu Mehari Beyene, Catholic University of Louvain. \code{<kasu.beyene@uclouvain.be>}
#'
#' Anouar El Ghouch, Catholic University of Louvain. \code{<anouar.elghouch@uclouvain.be>}
#' @references Beyene, K. M. and El Ghouch A. (2020). Smoothed time-dependent ROC curves for right-censored survival data. \emph{submitted}.
#' @examples library(cenROC)
#'
#' X <- rnorm(100) # random data vector
#' wt <- runif(100) # weight vector
#'
#' ## Normal reference bandwidth selection
#' NR(X = X, wt = wt)$bw
#'
#' @export

NR <- function (X, wt, ktype="normal") {
  nx <- length(X)
  mul <- (nx * sum(wt * wt)) / ((sum(wt)) ^ 2)
  stdv <- sqrt(wvar(X = X, wt = wt))
  IQR <- wIQR(X = X, wt = wt)
  sigma <- min(stdv, IQR / 1.349)
  c <- (4 * sqrt(pi) * (muro(ktype)$ro) / ((muro(ktype = ktype)$mu2) ^ 2)) ^ (1 / 3)
  wbw <- (c * sigma) * (mul ^ (1 / 3)) * nx ^ (-1 / 3)

  return(list(bw = wbw))
}


#' The plug-in bandwidth selection for weighted data
#'
#' @description This function computes the data-driven bandwidth for smoothing the ROC (or distribution) function using the PI method of Beyene and El Ghouch (2020). This is an extension of the classical (unweighted) direct plug-in bandwith selection method to the case of weighted data.
#'
#' @param X The numeric vector of random variable.
#' @param wt The non-negative weight vector.
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". By default, the "\code{normal}" kernel is used.
#' @return Returns the computed value for the bandwith parameter.
#' @details See Beyene and El Ghouch (2020) for details.
#' @author
#' Kassu Mehari Beyene, Catholic University of Louvain. \code{<kasu.beyene@uclouvain.be>}
#'
#' Anouar El Ghouch, Catholic University of Louvain. \code{<anouar.elghouch@uclouvain.be>}
#' @references Beyene, K. M. and El Ghouch A. (2020). Smoothed time-dependent ROC curves for right-censored survival data. \emph{submitted}.
#' @examples library(cenROC)
#'
#' X <- rnorm(100) # random data vector
#' wt <- runif(100) # weight vector
#'
#' ## Plug-in bandwidth selection
#' PI(X = X, wt = wt)$bw
#' 
#' @export

PI <- function(X, wt, ktype="normal")
{
  ### bandwidth estimation ###################
  band <- function (X, wt, psi, ktype) {
    n <- length(X)
    mul <- (n * sum(wt * wt)) / ((sum(wt)) ^ 2)
    #### Estimates of E(wt^2)/((E(wt))^2)
    co <- ((muro(ktype = ktype)$ro) / ((muro(ktype = ktype)$mu2) ^ 2)) ^ (1 / 3)
    wbww <- (co) * (mul ^ (1 / 3)) * (n ^ (-1 / 3)) * ((-psi) ^ (-1 / 3))

    return(wbww)
  }

  ####### Estimation of psi  ################
  psi <- function(r, g) {
    n <- length(X)
    w <- (n * wt) / (sum(wt)) #### Estimate for wt/E(wt)
    ww <- outer(w, w, "*")
    aux <- outer(X, X, "-") / g
    aux <- (ww) * (dnorkernel(r, aux))
    result <- (sum(aux)) * (((g) ^ (-r - 1)) * (n ^ (-2)))
    return(result)
  }
  g0 <- NR(X, wt, ktype = ktype)$bw ### Intial bandwidt estimatioh using normal reference
  psi2 <- psi(2, g0) ### Plug-in estimate
  PIbw <- band(X, wt, psi2, ktype = ktype)

  return(list(bw = PIbw))
}


#############################################################################
## This function is modified from kerdiest R package CVbw function ##########
##  of Quintela-del-Rio and Estevez-Perez (2015)                   ##########
#############################################################################

#' The cross-validation bandwidth selection for weighted data
#'
#' @description This function computes the data-driven bandwidth for smoothing the ROC (or distribution) function using the CV method of Beyene and El Ghouch (2020). This is an extension of the classical (unweighted) cross-validation bandwith selection method to the case of weighted data.
#'
#' @param X The numeric data vector.
#' @param wt The non-negative weight vector.
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". By default, the "\code{normal}" kernel is used.
#' @return Returns the computed value for the bandwith parameter.
#' @details Bowman et al (1998) proposed the cross-validation bandwidth selection method for unweighted kernal smoothed distribution function. This method is implemented in the \code{R} package \code{kerdiest}.
#' We adapted this for the case of weighted data by incorporating the weight variable into the cross-validation function of Bowman's method. See Beyene and El Ghouch (2020) for details.
#'
#' @author
#' Kassu Mehari Beyene, Catholic University of Louvain. \code{<kasu.beyene@uclouvain.be>}
#'
#' Anouar El Ghouch, Catholic University of Louvain. \code{<anouar.elghouch@uclouvain.be>}
#' @references Beyene, K. M. and El Ghouch A. (2020). Smoothed time-dependent ROC curves for right-censored survival data. \emph{submitted}.
#' @references Bowman A., Hall P. and Trvan T.(1998). Bandwidth selection for the smoothing of distribution functions. \emph{Biometrika} 85:799-808.
#' @references Quintela-del-Rio, A. and Estevez-Perez, G. (2015). \code{kerdiest:} Nonparametric kernel estimation of the distribution function, bandwidth selection and estimation of related functions. \code{R} package version 1.2.
#' @examples 
#' \dontrun{library(cenROC)
#' 
#' X <- rnorm(100) # random data vector
#' wt <- runif(100) # weight vector
#'
#' ## Cross-validation bandwidth selection
#' CV(X = X, wt = wt)$bw
#'
#' }
#' @export

CV <- function(X, wt, ktype = "normal")
{
  mul <- length(wt) / sum(wt)
  prob_quantile <- 0
  ss <- quantile(X, c(prob_quantile, 1 - prob_quantile))
  y <- seq(ss[1], ss[2], length.out = 100)
  seq_bws = seq((max(X) - min(X)) / 200, (max(X) - min(X)) / 2, length.out = 50)
  n_bws <- length(seq_bws)
  CVfunction <- numeric(length = n_bws)
  for (i in 1:n_bws)
  {
    integrand <- apply((t((mul * wt) * (t(outer(y, X, "-") >= 0))) - t(ker_dis_i(ktype = ktype, y = y, X = X, wt = wt, bw = seq_bws[i]))) ^ 2, 1, mean)
    CVfunction[i] <- integ(x = y, fx = integrand, method = "simps")
  }
  i0 <- which.min(CVfunction)
  CVbw_val <- seq_bws[i0]
  wbw <- CVbw_val
  return(list(bw = wbw))
}


#' Function to select the bandwidth parameter needed for smoothing the time-dependent ROC curve.
#'
#' @description This function computes the data-driven bandwidth value for smoothing the ROC curve.
#'              It contains three methods: the normal refrence, the plug-in and the cross-validation methods.
#'
#' @param X The numeric data vector.
#' @param wt The non-negative weight vector.
#' @param bw A character string specifying the bandwidth selection method. The possible options are "\code{NR}" for the normal reference, the plug-in "\code{PI}" and cross-validation "\code{CV}".
#' @param ktype A character string indicating the type of kernel function: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". Default is "\code{normal}" kernel.
#' @return Returns the estimated value for the bandwith parameter.
#' @author
#' Kassu Mehari Beyene, Catholic University of Louvain. \code{<kasu.beyene@uclouvain.be>}
#'
#' Anouar El Ghouch, Catholic University of Louvain. \code{<anouar.elghouch@uclouvain.be>}
#' @references Beyene, K. M. and El Ghouch A. (2020). Smoothed time-dependent ROC curves for right-censored survival data. \emph{submitted}.
#' @keywords internal

wbw <- function(X, wt, bw = "NR", ktype = "normal")
{
  if (is.numeric(bw)) {
    bwv <- bw
  } else if (bw == "CV") {
    bwv <- CV(X = X, wt = wt, ktype = ktype)$bw
  } else if (bw == "NR") {
    bwv <- NR(X = X, wt = wt, ktype = ktype)$bw
  } else if (bw == "PI") {
    bwv <- PI(X = X, wt = wt, ktype = ktype)$bw
  } else{
    print("Please check your bandwidth options")
  }
  return(list(bw = bwv))
}
