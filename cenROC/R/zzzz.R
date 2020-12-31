#' Weighted variance estimation
#'
#' @param X The numeric data vector.
#' @param wt The non-negative weight vector.
#' @param na.rm The character indicator wether to consider missing value(s) or not. The defult is FALSE.
#' @keywords internal

wvar <- function(X, wt, na.rm = FALSE) {
  if (na.rm) {
    wt <- wt[i <- !is.na(X)]
    X <- X[i]
  }
  wsum <- sum(wt)
  wmean = sum(wt * X) / wsum
  varr = sum(wt * (X - wmean) ^ 2) / (wsum)
  return(varr)
}

#' Weighted quartile estimation
#'
#' @param X The numeric data vector.
#' @param wt The non-negative weight vector.
#' @param p The percentile value. The defult is 0.5.
#' @keywords internal

wquantile <- function(X, wt, p = 0.5)
{
  if (!is.numeric(wt) || length(X) != length(wt))
    stop("X and wt must be numeric and equal-length vectors")
  if (!is.numeric(p) || any(p < 0 | p > 1))
    stop("Quartiles must be 0<=p<=1")
  if (min(wt) < 0)
    stop("Weights must be non-negative numbers")
  ord <- order(X)
  X <- X[ord]
  cusumw <- cumsum(wt[ord])
  sumW <- sum(wt)
  plist <- cusumw / sumW
  qua <- withCallingHandlers(approx(plist, X, p)$y, warning=function(w){invokeRestart("muffleWarning")})

  return(qua)
}

#' Weighted inter-quartile range estimation
#'
#' @param X The numeric data vector.
#' @param wt The non-negative weight vector.
#' @keywords internal

wIQR <- function(X, wt) {
  (wquantile(X = X, wt = wt, p = 0.75) - wquantile(X = X, wt = wt, p = 0.25))
}



#' Numerical Integral function using Simpson's rule
#'
#' @param x The numeric data vector.
#' @param fx The function.
#' @param n.pts Number of points.
#' @param method The character string specifying method of numerical integration. The possible options are \code{trap} for trapezoidal rule and \code{simps} for simpson'r rule.
#' @keywords internal

integ <- function(x, fx, method, n.pts = 256) {
  n = length(x)
  if (method == "simps") {
    if (class(fx) == "function")
      fx = fx(x)
    if (n != length(fx))
      stop("Unequal input vector lengths")
    if (n.pts < 64)
      n.pts = 64
    ap = approx(x, fx, n = 2 * n.pts + 1)
    h = diff(ap$x)[1]
    integral = h * (ap$y[2 * (1:n.pts) - 1] + 4 * ap$y[2 * (1:n.pts)] + ap$y[2 * (1:n.pts) + 1]) / 3
    value = sum(integral)
  }
  if (method == "trap") {
    if (!is.numeric(x) | !is.numeric(fx))
    {
      stop('The variable of integration "x" or "fx" is not numeric.')
    }
    if (length(x) != length(fx))
    {
      stop("The lengths of the variable of integration and the integrand do not match.")
    }
    # integrate using the trapezoidal rule
    integral <- 0.5 * sum((x[2:(n)] - x[1:(n - 1)]) * (fx[1:(n - 1)] + fx[2:n]))
    value <- integral
  }
  return(value)
}

#' Derivative of normal distribution
#'
#' @param X The numeric data vector.
#' @param ord The order of derivative.
#' @keywords internal

dnorkernel <- function(ord, X)
{
  if (ord == 2)
    # second derivative
    result <- (1 / (sqrt(2 * pi))) * exp(-(X ^ 2) / 2) * ((X ^ 2) - 1)
  else if (ord == 4)
    # fourth derivative
    result <- (1 / (sqrt(2 * pi))) * exp(-(X ^ 2) / 2) * (3 - (6 * (X ^ 2)) + X ^ 4)
  else if (ord == 6)
    # sixth derivative
    result <- (1 / (sqrt(2 * pi))) * exp(-(X ^ 2) / 2) * (X ^ 6 - (15 * (X ^ 4)) + (45 * (X ^ 2)) - 15)
  else if (ord == 8)
    # eighth derivative
    result <- (1 / (sqrt(2 * pi))) * exp(-(X ^ 2) / 2) * (X ^ 8 - (28 * (X ^ 6)) + (210 * (X ^ 4)) - (420 * (X ^ 2)) + 105)
  return(result)
}

#' Distribution function without the ith observation
#'
#' @param X The numeric data vector.
#' @param y The vector where the kernel estimation is computed.
#' @param wt The non-negative weight vector.
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}".
#' @param bw A numeric bandwidth value.
#' @return Returns the estimated value for the bandwith parameter.
#' @author Kassu Mehari Beyene  and Anouar El Ghouch
#' @keywords internal

ker_dis_i <- function(X, y, wt, ktype, bw)
{
  n <- length(X);
  AUX <- matrix(0, n, n);
  zero <- rep(0, n);
  ww <- outer(wt, zero, "-");
  diag(ww) <- 0;
  den <- apply(ww, 2, sum);
  resu <- matrix(0, n, length(y));
  for (j in 1:length(y))
  {
    AUX <- matrix(rep.int(outer(y[j], X, "-"), n), nrow = n, byrow = TRUE) / bw;
    aux <- kfunc(ktype = ktype, difmat = AUX );
    aux1 <- t(wt * t(aux));
    diag(aux1) <- 0;
    resu[, j] <- (apply(aux1, 1, sum)) / den;
  }
  return(resu)
}

#' The value of squared integral x^2 k(x) dx and integral x k(x) K(x) dx
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}".
#' @keywords internal

muro <- function(ktype)
{
  if (ktype == "normal") {
    ro <- 2 * 0.28209
    mu2 <- 1
  } else if (ktype == "epanechnikov") {
    ro <- 2 * 0.12857
    mu2 <- 1 / 5
  } else if (ktype == "biweight") {
    ro <- 2 * 0.10823
    mu2 <- 1 / 7
  } else if (ktype == "triweight") {
    ro <- 2 * 0.095183
    mu2 <- 1 / 9
  }

  return(list(ro = ro, mu2 = mu2))
}

#' Kernel distribution function
#'
#' @param X A numeric vector of sample data.
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}".
#' @return Returns a vector resulting from evaluating X.
#' @keywords internal

kfunction <- function(ktype, X) {
  if (ktype == "normal") {
    result <- pnorm(X)

  }
  else if (ktype == "epanechnikov") {
    result <- (0.75 * X * (1 - (X ^ 2) / 3) + 0.5)

  }
  else if (ktype == "biweight") {
    result <- ((15 / 16) * X - (5 / 8) * X ^ 3 + (3 / 16) * X ^ 5 + 0.5)

  }
  else if (ktype == "triweight") {
    result <- ((35 / 32) * X - (35 / 32) * X ^ 3 + (21 / 32) * X ^ 5 - (5 / 32) * X ^ 7 + 0.5)
  }
  return(result)
}

#' Function to evaluate the matrix of data vector minus the grid points divided by the bandwidth value.
#'
#' @param difmat A numeric matrix of sample data (X) minus evaluation points (x0) divided by bandwidth value (bw).
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". By default, the "\code{normal}" kernel is used.
#' @return Returns the matrix resulting from evaluating \code{difmat}.
#' @keywords internal

kfunc <- function(ktype = "normal", difmat)
{
  if (ktype == "normal")
  {
    estim <- kfunction(ktype = "normal", X = difmat)
  }
  else if (ktype == "epanechnikov")
  {
    estim <- difmat
    low <- (difmat <= -1)
    up <- (difmat >= 1)
    btwn <- (difmat > -1 & difmat < 1)
    estim[low] <- 0
    estim[up] <- 1
    value <- estim[btwn]
    estim[btwn] <- kfunction(ktype = "epanechnikov", X = value)
  }
  else if (ktype == "biweight")
  {
    estim <- difmat
    low <- (difmat <= -1)
    up <- (difmat >= 1)
    btwn <- (difmat > -1 & difmat < 1)
    estim[low] <- 0
    estim[up] <- 1
    value <- estim[btwn]
    estim[btwn] <- kfunction(ktype = "biweight", X = value)
  }
  else if (ktype == "triweight")
  {
    estim <- difmat
    low <- (difmat <= -1)
    up <- (difmat >= 1)
    btwn <- (difmat > -1 & difmat < 1)
    estim[low] <- 0
    estim[up] <- 1
    value <- estim[btwn]
    estim[btwn] <- kfunction(ktype = "triweight", X = value)
  }
  return(estim)
}

#'  ROC estimation function
#'
#' @param U The vector of grid points where the ROC curve is estimated.
#' @param D The event indicator.
#' @param M The numeric vector of marker values for which the time-dependent ROC curves is computed.
#' @param bw The bandwidth parameter for smoothing the ROC function. The possible options are \code{NR} normal reference method; \code{PI} plug-in method and \code{CV} cross-validation method. The default is the \code{NR} normal reference method.
#' @param method is the method of ROC curve estimation. The possible options are \code{emp} emperical metod; \code{untra} smooth without boundary correction and \code{tra} is smooth ROC curve estimation with boundary correction.
#' @param ktype A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}".
#' @author Beyene K. Mehari and El Ghouch Anouar
#' @references Beyene, K. M. and El Ghouch A. (2019). Smoothed time-dependent ROC curves for right-censored survival data. <\url{https://dial.uclouvain.be/pr/boreal/object/boreal:219643}>.
#' @keywords internal

RocFun <- function(U, D, M, bw = "NR", method, ktype) {
  oM <- order(M)
  D <- (D[oM])
  nD <- length(D)
  sumD <- sum(D)
  Z <- 1 - cumsum(1 - D) / (nD - sumD)
  AUC <- sum(D * Z) / sumD
  if (method == "emp") {
    difmat <- (outer(U, Z, "-"))
    resul <- (difmat >= 0)
    roc1 <- sweep(resul, 2, D, "*")
    roc <- apply(roc1, 1, sum) / sumD
    bw1 <- NA
  }
  else if (method == "untra") {
    Zt <- Z
    Ut <- U
    Ztt <- Zt[D != 0]
    wt <- D[D != 0]
    bw1 <- wbw(X = Ztt, wt = wt, bw = bw, ktype = ktype)$bw
    difmat <- (outer(Ut, Ztt, "-")) / bw1
    resul <- kfunc(ktype = ktype, difmat = difmat)
    w <- wt / sum(wt)
    roc1 <- sweep(resul, 2, w, "*")
    roc <- apply(roc1, 1, sum)
  }
  else if (method == "tra") {
    mul <- nD / (nD + 1)
    Zt <- qnorm(mul * Z + (1 / nD ^ 2))
    Ut <- qnorm(mul * U + (1 / nD ^ 2))
    Ztt <- Zt[D != 0]
    wt <- D[D != 0]
    bw1 <- wbw(X = Ztt, wt = wt, bw = bw, ktype = ktype)$bw
    difmat <- (outer(Ut, Ztt, "-")) / bw1
    resul <- kfunc(ktype = ktype, difmat = difmat)
    w <- wt / sum(wt)
    roc1 <- sweep(resul, 2, w, "*")
    roc <- apply(roc1, 1, sum)
  }
  else{
   stop("The specified method is not correct.")
  }
  return(list(roc = roc, auc = AUC, bw = bw1))
}


#' Survival probability conditional to the observed data estimation for right censored data.
#'
#' @param Y The numeric vector of event-times or observed times.
#' @param M The numeric vector of marker values for which we want to compute the time-dependent ROC curves.
#' @param censor The censoring indicator, \code{1} if event, \code{0} otherwise.
#' @param t A scaler time point at which we want to compute the time-dependent ROC curve.
#' @param h A scaler for the bandwidth of Beran's weight calculaions. The defualt is using the method of Sheather and Jones (1991).
#' @param kernel A character string giving the type kernel to be used: "\code{normal}", "\code{epanechnikov}", , "\code{tricube}", "\code{boxcar}", "\code{triangular}", or "\code{quartic}". The defaults is "\code{normal}" kernel density.
#' @return Return a vectors:
#' @return \code{positive    }    \code{P(T<t|Y,censor,M)}.
#' @return \code{negative    }     \code{P(T>t|Y,censor,M)}.
#' @references Beyene, K. M. and El Ghouch A. (2019). Smoothed time-dependent ROC curves for right-censored survival data. <\url{https://dial.uclouvain.be/pr/boreal/object/boreal:219643}>.
#' @references Li, Liang, Bo Hu and Tom Greene (2018).  A simple method to estimate the time-dependent receiver operating characteristic curve and the area under the curve with right censored data, Statistical Methods in Medical Research, 27(8): 2264-2278.
#' @references Pablo Martínez-Camblor and Gustavo F. Bayón and Sonia Pérez-Fernández (2016). Cumulative/dynamic roc curve estimation, Journal of Statistical Computation and Simulation, 86(17): 3582-3594.
#' @keywords internal

Csurv <- function(Y, M, censor, t, h = NULL, kernel="normal") {
  if (is.null(h)) {
    h <- bw.SJ(M, method = "dpi")
  }
  if(kernel=="normal"){
    kernel <- "gaussian"
  }
  n <- length(M)
  positive <- rep(NA, n)
  for (i in 1:n) {
    if (Y[i] > t) {
      positive[i] <- 0

    } else {
      if (censor[i] == 1) {
        positive[i] <- 1

      } else {
        St <- Beran(time = Y, status = censor, covariate = M, x = M[i], y = t, kernel = kernel, bw = h)
        Sy <- Beran(time = Y, status = censor, covariate = M, x = M[i], y = Y[i], kernel = kernel, bw = h)
        if (Sy == 0) {
          positive[i] <- 1

        } else {
          positive[i] <- 1 - St / Sy

        }
      }
    }
  }
  negative <- 1 - positive

  return(list(positive = positive, negative = negative))

}


# Function to compute the knots.
# This functions are based on the R package intcensROC.
.knotT <- function(U, V, delta, dim) {
  size <- dim - 2
  knot_pre_t = c(U[delta != 1], U[delta != 2], V[delta != 2], V[delta != 3])
  qt = rep(0, size + 1)
  qt[1] = 0
  for (i in 2:(size + 1)) qt[i] = quantile(knot_pre_t, (i - 1)/size, name = F,
                                           na.rm = TRUE)
  knots = c(qt[1], qt[1], qt, qt[size + 1], qt[size + 1])
}

.knotM <- function(marker, dim) {
  size <- dim - 2
  knot_pre_m = marker
  qt = rep(0, size + 1)
  qt[1] = 0
  for (i in 2:(size + 1)) qt[i] = quantile(knot_pre_m, (i - 1)/size, name = F,
                                           na.rm = TRUE)
  qt[size + 1] = max(marker + 0.1)
  knots = c(qt[1], qt[1], qt, qt[size + 1], qt[size + 1])
}


#'  Compute the conditional survival function for Interval Censored Survival Data
#'
#' @description  A method to compute the survival function for the
#' interval censored survival data based on a spline function based constrained
#' maximum likelihood estimator. The maximization process of likelihood is
#' carried out by generalized gradient projection method.
#' @usage condS(L, R, M, Delta, t, m)
#' @param L The numericvector of left limit of observed time. For left censored observations \code{L == 0}.
#' @param R The numericvector of right limit of observed time. For right censored observation \code{R == inf}.
#' @param M  An array contains marker levels for the samples.
#' @param Delta An array of indicator for the censored type, use 1, 2, 3 for
#' event happened before the left bound time, within the defined time range, and
#' after.
#' @param t A scalar indicates the predict time.
#' @param m A scalar for the cutoff of the marker variable.
#' @references Wu, Yuan; Zhang, Ying. Partially monotone tensor spline estimation
#' of the joint distribution function with bivariate current status data.
#' Ann. Statist. 40, 2012, 1609-1636 <doi:10.1214/12-AOS1016>
#' @keywords internal

condS <- function(L, R, M, Delta=NULL, t, m) {
  n <- length(L)
  U <- L
  V <- R
  Marker <- M
  PredictTime <- t
  ind <- (U<=0)
  ind1 <- (V==Inf)
  if (any(ind1 == TRUE)){
    V[ind1] <- 10000000
  }
  if(is.null(Delta)){
    Delta <- rep(2, n)
    Delta[ind] <- 1
    Delta[ind1] <- 3
  }
  if (any(Marker < 0))
    stop(paste0("Negative marker value found!"))

  # detemine the dimension of spline function
  size <- length(U)
  cadSize <- size^(1/3)
  if (cadSize - floor(cadSize) < ceiling(cadSize) - cadSize) {
    Dim <- floor(cadSize) + 2
  } else {
    Dim <- ceiling(cadSize) + 2
  }

  # compute the knots for time and marker
  knotT <- .knotT(U, V, Delta, Dim)
  knotM <- .knotM(Marker, Dim)

  # compute the thetas
  theta = .Call("_cenROC_sieve", PACKAGE = "cenROC", U, V, Marker, Delta,
                knotT, knotM, Dim)

  m2 <- m - 0.0001
  m  <- m + 0.0001
  Fm = .Call("_cenROC_surva", PACKAGE = "cenROC", theta, m,
             m2, PredictTime, knotT, knotM)

  Fest <- 1 - Fm

  return(Fest)
}



#' Survival probability conditional on the observed data estimation for interval censored data
#'
#' @param L The numericvector of left limit of observed time. For left censored observations \code{L == 0}.
#' @param R The numericvector of right limit of observed time. For right censored observation \code{R == inf}.
#' @param M The numeric vector of marker value.
#' @param t A scaler time point used to calculate the the ROC curve
#' @param method A character indication type of modeling. This include nonparametric \code{"np"},parmetric \code{"pa"} and semiparametric \code{"sp"}.
#' @param dist A character incating the type of distribution for parametric model. This includes are \code{"exponential"}, \code{"weibull"}, \code{"gamma"}, \code{"lnorm"}, \code{"loglogistic"} and \code{"generalgamma"}.
#' @return Return a vectors:
#' @return \code{positive    }    \code{P(T<t|L,R,M)}.
#' @return \code{negative    }     \code{P(T>t|L,R,M)}.
#' @references Beyene, K. M. and El Ghouch A. (2020). Time-dependent ROC curves estimator for interval-censored survival data.
#' @keywords internal

ICsur <- function( L, R, M,  t,  method, dist) {
  data <- data.frame(L=L, R=R, M=M)
  n <- length(M) ;
  positive <- rep(NA, n);
  for (i in 1:n) {
    if (R[i] <= t) {
      positive[i] <- 1;
    } else {
      if (L[i] >= t) {
        positive[i] <- 0;
      } else {
        if (method=="np"){
          tmp1 <- condS(L=L, R=R, M=M, t=t, m=M[i])
          tmp2 <- condS(L=L, R=R, M=M, t=L[i], m=M[i])
          tmp3 <- condS(L=L, R=R, M=M, t=R[i], m=M[i])
          tmp <- c(tmp1, tmp2, tmp3)
        } else if (method=="pa"){
          formula <- Surv(time=L, time2=R, type="interval2") ~ M
          fit <- ic_par(formula, model = "aft", dist = dist, data=data, weights = NULL)
          newdat <- data.frame(M=c(M[i]));
          tmp <- 1 - (getFitEsts(fit, newdat, q=c(t, L[i], R[i])));
        } else if (method=="sp"){
          formula <- Surv(time=L, time2=R, type="interval2") ~ M
          fit <- ic_sp(formula, model = "ph", data=data, weights = NULL)
          newdat <- data.frame(M=c(M[i]));
          tmp <- 1 - (getFitEsts(fit, newdat, q=c(t, L[i], R[i])));
        }
        positive[i] <- ifelse(R[i]==Inf, 1-(tmp[1]/tmp[2]), ifelse(L[i]==0, (1-tmp[1])/(1-tmp[3]), (tmp[2]-tmp[1])/(tmp[2]-tmp[3])))
      }
    }
  }
  negative <- 1 - positive;
  return(list(positive = positive, negative = negative));
}
