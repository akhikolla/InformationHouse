# Generic function for computing rq-residuals for count data from a CDF.
rqRes = function(y, pFun, ...) {
  b = pFun(y, ...)
  a = b*0
  a = pFun(y - 1, ...)
  if (is.null(dim(y))) {
    res = numeric(length(y)) + NA
  } else {
    res = array(data = NA, dim = dim(y))
  }
  isnna = which(!is.na(y) & !is.na(a) & !is.na(b))
  res[isnna] = rqRes0(a[isnna], b[isnna])
  # if (any(is.infinite(res))) {
  #   warning("Some residuals infinite.")
  # }
  res
}

rqRes0 = function(a, b) {
  stopifnot(length(a) == length(b))
  stats::qnorm(runif(length(a), a, b))
}

#' Plot residuals against fitted values
#'
#' Plots randomized-quantile residuals for binomial N-mixture models against fitted values.
#'
#' @param umFit An object from a model fitted using \link[unmarked]{pcount}.
#' @param type The type of randomized quantile residual to plot. One of 'marginal', 'site-sum' or 'observation'.
#' @param ... Plot arguments.
#'
#' @export
#'
#' @examples
#' library(unmarked)
#' umf = unmarkedFramePCount(y = shoveler$y, obsCovs = shoveler$obs, siteCovs = shoveler$site)
#' fmP = pcount(~scale(date) + scale(reedcover) ~ scale(log(water)) + scale(latitude), 
#'       data = umf, K = 80)
#' residfit(fmP, "marginal")
#' residfit(fmP, "site-sum")
#' residfit(fmP, "observation")
residfit = function(umFit, type = "marginal", ...) {
  if (!inherits(umFit,"unmarkedFitPCount")) {
    stop("Argument needs to be a model fitted by the function pcount in the package unmarked.")
  }
  type = match.arg(type, c("marginal", "site-sum", "observation"))
  if (identical(type, "marginal")) {
    rqr = rqresiduals(umFit, type = "marginal")
    fitted = fitted(umFit, K = umFit@K)
  } else if (identical(type, "site-sum")) {
    rqr = rqresiduals(umFit, type = "site-sum")
    fitted = apply(fitted(umFit, K = umFit@K), 1, sum, na.rm = TRUE)
  } else {
    rqr = rqresiduals(umFit, type = "observation")
    fitted = fitted(umFit, K = umFit@K)
  }
  plotArgs = list(fitted, rqr, ...)
  if (!("ylab" %in% names(plotArgs)))
      plotArgs = c(plotArgs, list(ylab = paste(type, "rq-residuals")))
  if (!("xlab" %in% names(plotArgs)))
    plotArgs = c(plotArgs, list(xlab = "fitted values"))
  do.call(plot, plotArgs)
}

#' Qq plot of randomized quantile residuals against standard normal quantiles
#'
#' @param umFit An object of class \link[unmarked]{unmarkedFit} from a model fitted using \link[unmarked]{pcount}.
#' @param type The type of randomized quantile residual to plot. One of 'site-sum' or 'observation'.
#' @param main Plot label.
#' @param plotLine If true, the identity line is added to the plot.
#' @param ... Further arguments passed to qqnorm.
#'
#' @return A list with x and y coordinates of the qq plot, see \link[stats]{qqnorm}.
#' @export
#'
#' @examples
#' library(unmarked)
#' umf = unmarkedFramePCount(y = shoveler$y, obsCovs = shoveler$obs, siteCovs = shoveler$site)
#' fmP = pcount(~scale(date) + scale(reedcover) ~ scale(log(water)) + scale(latitude), 
#'       data = umf, K = 80)
#' residqq(fmP, "site-sum")
#' residqq(fmP, "observation")
residqq =  function(umFit, type = "site-sum", main = "Residual qq plot", plotLine = TRUE, ...) {
    if (!inherits(umFit,"unmarkedFitPCount")) {
      stop("Argument needs to be a model fitted by the function pcount in the package unmarked.")
    }
    type = match.arg(type, c("marginal", "site-sum", "observation"))
    if (identical(type, "marginal")) {
      stop("Marginal qq-plot not implemented")
    } else if (identical(type, "site-sum")) {
      rqr = rqresiduals(umFit, type = "site-sum")
    } else {
      rqr = rqresiduals(umFit, type = "observation")
    }
    qqArgs = list(rqr ,main = main, ...)
    if (!("ylab" %in% names(qqArgs)))
      qqArgs = c(qqArgs , ylab = "residual quantile")
    qq = do.call(stats::qqnorm, qqArgs)
    if (plotLine)
      abline(a = 0, b = 1)
    invisible(qq)
}

#' Randomized quantile resiudals for binomial N-mixture models.
#'
#' Computes three types of randomized quantile residuals for binomial N-mixture models.
#' @param umFit An object of class \link[unmarked]{unmarkedFit} from a model fitted using \link[unmarked]{pcount}.
#' @param type The type of rq residuals to compute, one of 'marginal', 'site-sum' or 'observation'.
#'
#' @return A matrix (if \code{type} is 'marginal' or 'site-sum') or vector (for )  con.
#' @export
#'
#' @examples
#' library(unmarked)
#' umf = unmarkedFramePCount(y = shoveler$y, obsCovs = shoveler$obs, siteCovs = shoveler$site)
#' fmP = pcount(~scale(date) + scale(reedcover) ~ scale(log(water)) + scale(latitude), 
#'       data = umf, K = 80)
#' qqnorm(rqresiduals(fmP, "s"))
#' qqnorm(rqresiduals(fmP, "o"))
#' par(mfcol = c(3,4))
#' invisible(apply(rqresiduals(fmP, "m"), 2, qqnorm))
rqresiduals = function(umFit, type = "marginal") {
  if (!inherits(umFit,"unmarkedFitPCount")) {
    stop("Argument needs to be a model fitted by the function pcount in the package unmarked.")
  }
  type = match.arg(type, c("marginal", "site-sum", "observation"))
  if (identical(type, "marginal"))
    res = rqResM(umFit)
  else if (identical(type, "site-sum"))
    res = rqResS(umFit)
  else
    res = rqResO(umFit)
  res
}


rqResM = function(umFit) {
  fitval = unmarked::fitted(umFit, K = umFit@K) # A bug in older versions of unmarked (fixed now) may cause incorrect fitted values for NB models unless K is supplied.
  y = umFit@data@y
  if (identical(umFit@mixture, "P")) {
    rqr = rqRes(y, pFun = stats::ppois, lambda = fitval)
  } else if (identical(umFit@mixture, "NB")) {
    if (!identical(umFit@estimates["alpha"]@invlink, "exp"))
      stop("Unknown link function.")
    size = exp(umFit@estimates["alpha"]@estimates)
    rqr = rqRes(y, pFun = stats::pnbinom, mu = fitval, size = size)
  } else if (identical(umFit@mixture, "ZIP")) {
    if (!identical(umFit@estimates["psi"]@invlink, "logistic"))
      stop("Unknown link function.")
    pZIP = function(y, lambda, psi) {
      psi*(y >= 0) + (1 - psi) * stats::ppois(y, lambda/(1 - psi))
    }
    psi = stats::plogis(umFit@estimates["psi"]@estimates)
    rqr = rqRes(y, pFun = pZIP, lambda = fitval, psi = psi)
  } else {stop("Mixture not recognized.")}
  rqr
}

rqResO = function(umFit) {
  rN = integer(nrow(umFit@data@y)) + NA
  if (length(umFit@sitesRemoved) > 0)
    rN[-umFit@sitesRemoved] = apply(unmarked::ranef(umFit)@post[,,1], 1 , sample, x = umFit@K + 1, size = 1, replace = FALSE) - 1
  else
    rN = apply(unmarked::ranef(umFit)@post[,,1], 1 , sample, x = umFit@K + 1, size = 1, replace = FALSE) - 1
  p = unmarked::getP(umFit, na.rm = FALSE)
  res = rqRes(umFit@data@y, pFun = stats::pbinom, size = kronecker(rN, t(rep(1, ncol(p)))), prob = p)
  # if (any(is.infinite(res))) {
  #   warning(paste(sum(is.infinite(res)), " residuals infinite."))
  # }
  res
}


rqResS = function(umFit) {
  lam = unmarked::predict(umFit, type = "state", na.rm = FALSE)[,1]
  p = unmarked::getP(umFit, na.rm = FALSE)
  res = numeric(nrow(umFit@data@y)) + NA
  if (length(umFit@sitesRemoved) > 0) {
    y = umFit@data@y[-umFit@sitesRemoved,]
    p = p[-umFit@sitesRemoved,]
    lam = lam[-umFit@sitesRemoved]
  }
  else
    y = umFit@data@y
  cumProb = matrix(0, nrow = nrow(y), ncol = 2)
  dfun = switch(umFit@mixture,
                P = function(N) {stats::dpois(N, lam)},
                NB = function(N) {stats::dnbinom(N, mu = lam, size = exp(unmarked::coef(umFit, type = "alpha")))},
                ZIP = function(N) {psi = stats::plogis(unmarked::coef(umFit, type = "psi"))
                                   (1 - psi)*stats::dpois(N, lam/(1 - psi)) + psi*(N == 0)}
         )
  for (N in 0:umFit@K) {
    cumProb = cumProb + kronecker(dfun(N), t(c(1,1))) * pbinsum(y, rep(N, nrow(y)), p)[,2:3]
  }
  if (length(umFit@sitesRemoved) > 0)
    res[-umFit@sitesRemoved] = rqRes0(cumProb[,1], cumProb[,2])
  else
    res = rqRes0(cumProb[,1], cumProb[,2])
  # if (any(is.infinite(res))) {
  #   warning(paste(sum(is.infinite(res)), " residuals infinite."))
  # }
  res
}



#' Overdispersion metrics for binomial N-mixture models.
#'
#' Computes various types of overdispersion metrics, based on Pearson residuals, for binomial N-mixture models.
#' @param umFit An object of class \link[unmarked]{unmarkedFit} from a model fitted using \link[unmarked]{pcount}.
#' @param type The type of metric to compute, one of 'marginal', 'site-sum' or 'observation'.
#'
#' @return An estimate of overdispersion relative to the fitted model.
#' @export
#'
#' @examples
#' library(unmarked)
#' data(mallard)
#' fm.mallard <- pcount(~ 1 ~ 1, unmarkedFramePCount(y = mallard.y), K=100)
#' chat(fm.mallard, "m")
#' chat(fm.mallard, "s")
#' chat(fm.mallard, "o")
chat = function(umFit, type = "marginal") {
  if (!inherits(umFit,"unmarkedFitPCount")) {
    stop("Argument needs to be a model fitted by the function pcount in the package unmarked.")
  }
  type = match.arg(type, c("marginal", "site-sum", "observation"))
  if (identical(type, "marginal"))
    chat = chatM(umFit)
  else if (identical(type, "site-sum"))
    chat = chatS(umFit)
  else
    chat = chatO(umFit)
  chat
}

chatS = function(umFit) {
  p = unmarked::getP(umFit, na.rm = FALSE)
  lam = unmarked::predict(umFit, "state", na.rm = FALSE)[,1]
  if (length(umFit@sitesRemoved) > 0) {
    y = umFit@data@y[-umFit@sitesRemoved,]
    expected = unmarked::fitted(umFit, K = umFit@K)[-umFit@sitesRemoved,]
    lam = lam[-umFit@sitesRemoved]
    p = p[-umFit@sitesRemoved, ]
  } else {
    y = umFit@data@y
    expected = unmarked::fitted(umFit, K = umFit@K)
  }
  naMat = (is.finite(y + expected))
  naMat[which(naMat != 1)] = NA
  obs.site = apply(y * naMat, 1, sum, na.rm = TRUE)
  exp.site = apply(expected * naMat, 1, sum, na.rm = TRUE)
  odPar = switch(umFit@mixture, P = NA, ZIP = unmarked::coef(umFit)["psi(psi)"], NB  = unmarked::coef(umFit)["alpha(alpha)"])
  nVar = switch(umFit@mixture, P = lam, ZIP = lam*(1 + lam * stats::plogis(odPar)/(1 - stats::plogis(odPar))), NB = lam + lam^2 * exp(-odPar))
    var.site = lam * apply(p * (1 - p) * naMat, 1, sum, na.rm = TRUE) +
    nVar * apply(p * naMat, 1, function(row) {sum(outer(row, row), na.rm = TRUE)})
  chi2 = sum((obs.site - exp.site)^2/var.site)
  chi2/(nrow(y) - length(unmarked::coef(umFit)))
}


chatO = function(umFit) {
  p = unmarked::getP(umFit, na.rm = FALSE)
  if (length(umFit@sitesRemoved) > 0) {
    y = umFit@data@y[-umFit@sitesRemoved,]
    p = p[-umFit@sitesRemoved, ]
  } else {
    y = umFit@data@y
  }
  rN = apply(unmarked::ranef(umFit)@post[,,1], 1 , sample, x = umFit@K + 1, size = 1, replace = FALSE) - 1
  naMat = (is.finite(y + p))
  naMat[which(naMat != 1)] = NA
  obs.site = apply(y * naMat, 1, sum, na.rm = TRUE)
  exp.site = apply(p * naMat, 1, sum, na.rm = TRUE) * rN
  var.site = apply(p * (1 - p) * naMat, 1, sum, na.rm = TRUE) * rN
  posY = which(obs.site > 0)
  chi2 = sum((obs.site[posY] - exp.site[posY])^2/var.site[posY])
  chi2/(length(posY) - length(unmarked::coef(umFit)))
}

chatM = function(umFit) {
  odPar = switch(umFit@mixture, P = NA, ZIP = unmarked::coef(umFit)["psi(psi)"], NB  = unmarked::coef(umFit)["alpha(alpha)"])
  if (length(umFit@sitesRemoved) > 0) {
    y = umFit@data@y[-umFit@sitesRemoved,]
    expected = unmarked::fitted(umFit, K = umFit@K)[-umFit@sitesRemoved,]
  } else {
    y = umFit@data@y
    expected = unmarked::fitted(umFit, K = umFit@K)
  }
  ## For ZIP, zero inflation has been included in expected.
  variance = switch(umFit@mixture, P = expected, ZIP = expected*(1 + expected * stats::plogis(odPar)/(1 - stats::plogis(odPar))), NB = expected + expected^2 * exp(-odPar))
  res = (y - expected)/sqrt(variance)
  chi2 = sum(res^2, na.rm = TRUE)
  chi2/(sum(!is.na(res)) - length(unmarked::coef(umFit)))
}


#' Plot residuals against covariates
#'
#'
#' A convenience function to plot rq residuals against all untransformed numeric covariates.
#' Site-sum randomized quantile residuals are used for site covariates while marginal residuals are used for observation covariates.
#' The same random residual draws are reused for different covariates.
#'
#' @param umFit An object of class \link[unmarked]{unmarkedFit} from a model fitted using \link[unmarked]{pcount}.
#' @param ... Plot arguments.
#'
#' @export
#'
#' @examples
#' library(unmarked)
#' umf = unmarkedFramePCount(y = shoveler$y, obsCovs = shoveler$obs, siteCovs = shoveler$site)
#' fmP = pcount(~scale(date) + scale(reedcover) ~ scale(log(water)) + scale(latitude), 
#'       data = umf, K = 80)
#' residcov(fmP)
residcov = function(umFit, ...) {
  data = unmarked::getData(umFit)
  if (!is.null(data@siteCovs)) {
    rqS = rqresiduals(umFit, type = "site-sum")
    for (i in 1:ncol(data@siteCovs)) {
      if (is.numeric(data@siteCovs[[i]])) {
        plot(data@siteCovs[,i], rqS, xlab = colnames(data@siteCovs)[i],ylab = "site-sum rq resiudal", ...)
      }
    }
  }
  if (!is.null(data@obsCovs)) {
    rqM = rqresiduals(umFit, type = "marginal")
    for (i in 1:ncol(data@obsCovs)) {
      if (is.numeric(data@obsCovs[[i]])) {
        plot(data@obsCovs[,i], rqM, xlab = colnames(data@obsCovs)[i],ylab = "marginal rq resiudal", ...)
      }
    }
  }
}

