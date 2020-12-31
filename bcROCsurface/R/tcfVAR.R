####==========================================================================####
## The R code for the asymptotic covariance matrix of true class fractions.     ##
##                                                                              ##
##                                                                              ##
####==========================================================================####
#'
#' @title Asymptotic variance-covariance estimation for True Class Fractions (TCFs) at the cut point \eqn{(c_1, c_2)}
#'
#' @description
#' \code{asyCovTCF} computes the asymptotic variance-covariance matrix of full data (FULL) and bias-corrected estimators (i.e. full imputation, mean score imputation, inverse probability weighting, semiparametric efficient and K nearest neighbor) of TCFs.
#'
#'
#' @param obj_tcf  a result of a call to \code{\link{ROCs.tcf}}.
#' @param T  a numeric vector containing the diagnostic test values. \code{NA} values of \code{T} are not accepted.
#' @param Dvec a n * 3  binary matrix with three columns, corresponding to the three classes of the disease status. In row i, 1 in column j indicates that the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param V  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param rhoEst  a result of a call to \code{\link{rhoMLogit}} of \code{\link{rhoKNN}} to fit the disease model.
#' @param piEst  a result of a call to \code{\link{psglm}} to fit the verification model.
#' @param BOOT  a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a bootstrap resampling is employed to estimate the asymptotic variance-covariance matrix of bias-corrected TCFs.
#' @param nR  the number of bootstrap replicates, used when \code{BOOT = TRUE} or for FULL estimator. Usually this will be a single positive integer. Default 250.
#' @param parallel  a logical value. If \code{TRUE}, a parallel computing is employed in the bootstrap resampling process.
#' @param ncpus  number of processes to be used in parallel computing. Default is half of available cores.
#'
#' @details For bias-corrected estimators of TCFs, the asymptotic variance-covariance matrix at a fixed cut point is estimated by using the Delta method. The function \code{asyCovTCF} implements the explicit forms presented in To Duc et al. (2016a, 2016b). In addition, the bootstrap procedure is also available.
#'
#' For FULL estimator, the asymptotic variance-covariance matrix is computed via bootstrap only.
#'
#' @return This function returns an estimated asymptotic variance-covariance matrix for FULL estimator and bias-corrected estimators of TCFs at a fixed cut point.
#'
#' @references
#'
#' To Duc, K., Chiogna, M. and Adimari, G. (2016a)
#' Bias-corrected methods for estimating the receiver operating characteristic surface of continuous diagnostic tests.
#' \emph{Electronic Journal of Statistics}, \bold{10}, 3063-3113.
#'
#' To Duc, K., Chiogna, M. and Adimari, G. (2018)
#' Nonparametric estimation of ROC surfaces in presence of verification bias.
#' \emph{REVSTAT Statistical Journal}. Accepted.
#'
#'
#' @examples
#' data(EOC)
#'
#' # FULL data estimator
#' Dfull <- preDATA(EOC$D.full, EOC$CA125)
#' Dvec.full <- Dfull$Dvec
#'
#' full.tcf <- ROCs.tcf("full", T = EOC$CA125, Dvec = Dvec.full, cps = c(2, 4))
#' full.var <- asyCovTCF(full.tcf, T = EOC$CA125, Dvec = Dvec.full)
#'
#' # Preparing the missing disease status
#' Dna <- preDATA(EOC$D, EOC$CA125)
#' Dfact.na <- Dna$D
#' Dvec.na <- Dna$Dvec
#'
#' rho.out <- rhoMLogit(Dfact.na ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#'
#' ## FI estimator
#' fi.tcf <- ROCs.tcf("fi", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V,
#'                    rhoEst = rho.out, cps = c(2,4))
#' fi.var <- asyCovTCF(fi.tcf, T = EOC$CA125, Dvec = Dvec.na, V = EOC$V,
#'                     rhoEst = rho.out)
#'
#' ## MSI estimator
#' msi.tcf <- ROCs.tcf("msi", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V,
#'                     rhoEst = rho.out, cps = c(2,4))
#' msi.var <- asyCovTCF(msi.tcf, T = EOC$CA125, Dvec = Dvec.na, V = EOC$V,
#'                      rhoEst = rho.out)
#'
#' ## IPW estimator
#' pi.out <- psglm(V ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#'
#' ipw.tcf <- ROCs.tcf("ipw", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V,
#'                     piEst = pi.out, cps = c(2,4))
#' ipw.var <- asyCovTCF(ipw.tcf, T = EOC$CA125, Dvec = Dvec.na, V = EOC$V,
#'                      piEst = pi.out)
#'
#' ## SPE estimator
#' spe.tcf <- ROCs.tcf("spe", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V,
#'                    rhoEst = rho.out, piEst = pi.out, cps = c(2,4))
#' spe.var <- asyCovTCF(spe.tcf, T = EOC$CA125, Dvec = Dvec.na, V = EOC$V,
#'                      rhoEst = rho.out, piEst = pi.out)
#'
#' ## KNN estimators
#' XX <- cbind(EOC$CA125, EOC$CA153, EOC$Age)
#' rho.1nn <- rhoKNN(X = XX, Dvec = Dvec.na, V = EOC$V, K = 1, type = "mahala")
#' knn.tcf <- ROCs.tcf("knn", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V,
#'                    rhoEst = rho.1nn, cps = c(2,4))
#' knn.var <- asyCovTCF(knn.tcf, T = EOC$CA125, Dvec = Dvec.na, V = EOC$V,
#'                      rhoEst = rho.1nn)
#'
#'
#' @import parallel
#' @export
asyCovTCF <- function(obj_tcf, T, Dvec, V = NULL, rhoEst = NULL, piEst = NULL,
                      BOOT = FALSE, nR = 250, parallel = FALSE,
                      ncpus = ifelse(parallel, detectCores()/2, NULL)){
  if(is.list(obj_tcf) & length(obj_tcf) > 1) {
    check.tcf <- function(x){
      class(x) == "tcfs"
    }
    if(all(sapply(obj_tcf, check.tcf))) stop("The \"obj_tcf\" is a list containing results of TCFs at more than 1 cut point! \n Please, choose one result of \"obj_tcf\" \n")
    else stop("The argument \"obj_tcf\" is not a result of ROC.tcf()")
  }
  if(class(obj_tcf) != "tcfs") stop("The argument \"obj_tcf\" is not a result of a call to ROC.tcf()")
  ## checking the argument T
  if(missing(T)) stop("argument \"T\" is missing \n")
  if(!inherits(T, "numeric") | any(is.na(T))) stop("variable \"T\" must be a numeric vector and not include NA values")
  ## checking Dvec
  if(missing(Dvec)) stop("argument \"Dvec\" is missing \n")
  if(!inherits(Dvec, "matrix")| ncol(Dvec) != 3 | !all(is.element(na.omit(Dvec), c(0,1)))) stop("variable \"Dvec\" must be a binary matrix with 3 columns")
  if(length(T) != nrow(Dvec)) stop(gettextf("arguments imply differing number of observation: %d", length(T)), gettextf(", %d", nrow(Dvec)), domain = NA)
  ##
  method <- tolower(attr(obj_tcf, "name"))
  if(!is.element(method, c("full", "fi")) & is.null(V)) stop("argument \"V\" is missing \n")
  if(is.element(method, c("fi", "msi", "knn", "spe")) & is.null(rhoEst)) stop("argument \"rhoEst\" is missing \n")
  if(is.element(method, c("ipw", "spe")) & is.null(piEst)) stop("argument \"piEst\" is missing \n")
  ##
  if(!is.null(rhoEst)){
    if(!is.element(class(rhoEst),  c("prob_dise", "prob_dise_knn")))
      stop("\"rhoEst\" not a result of rhoMlogit or rhoKNN \n")
    if(is.element(class(rhoEst),  c("prob_dise_knn"))){
      if(is.null(rhoEst$K)) stop("The \"rhoEst\" is a list and contains results of rhoKNN corresponding to many K! \n Please, choose one result of them \n")
    }
  }
  if(!is.null(piEst)){
    if(!is.element(class(piEst),  c("prob_veri")))
      stop("\"piEst\" not a result of psglm \n")
  }
  ##
  Dvec.temp <- Dvec
  if(any(is.na(Dvec))) Dvec.temp[is.na(Dvec)] <- 99
  tcf.orgi <- obj_tcf
  tcf.thet <- attr(tcf.orgi, "theta")
  tcf.bet <- attr(tcf.orgi, "beta")
  cp <- attr(tcf.orgi, "cp")
  if(method == "full"){
    bst.full <- function(dt, inds, cp){
      dat <- dt[inds, ]
      ROCs.tcf(method = "full", T = dat[,1], Dvec = as.matrix(dat[, c(2:4)]),
               cps = cp)
    }
    data <- data.frame(T, Dvec)
    if(parallel){
      res.bst <- boot(data, bst.full, R = nR, cp = cp, parallel = "snow",
                      ncpus = ncpus)
    }
    else res.bst <- boot(data, bst.full, R = nR, cp = cp)
    ans <- var(res.bst$t)
  }
  else if(method == "knn"){
    if(!BOOT){
      XX <- rhoEst$X
      rho.temp <- rhoKNN(XX, Dvec, V, K = 2, type = rhoEst$type)$values
      pi.est <- psknn(XX, V, rhoEst$type)
      ans <- asy.Cov.KNN(T, tcf.thet, tcf.bet, cp, pi.est, rho.temp, rhoEst$K)
    }
    else{
      bst.knn <- function(dt, inds, cp, k, type){
        dat <- dt[inds,]
        XX <- as.matrix(dat[, -c(1:5)])
        rho.knn <- rhoKNN(XX, as.matrix(dat[, c(2:4)]), dat[, 5], K = k,
                          type = type)
        ROCs.tcf(method = "knn", T = dat[,1], Dvec = as.matrix(dat[, c(2:4)]),
                 V = dat[, 5], rhoEst = rho.knn, cps = cp)
      }
      data <- data.frame(T, Dvec, V, rhoEst$X)
      if(parallel){
        res.bst <- boot(data, bst.knn, R = nR, cp = cp, k = rhoEst$K,
                        type = rhoEst$type, parallel = "snow", ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.knn, R = nR, cp = cp, k = rhoEst$K,
                           type = rhoEst$type)
      ans <- var(res.bst$t)
    }
  }
  else if(method == "fi"){
    if(!BOOT){
      term1 <- EstFuncIE_deriv(V, T, rhoEst$X, cp, rhoEst$coeff, rhoEst$values,
                               rhoEst$Hess, tcf.thet, tcf.bet, m = 0)
      term2 <- EstFuncIE(Dvec.temp, V, T, rhoEst$X, cp, rhoEst$coeff,
                         rhoEst$values, tcf.thet, tcf.bet, m = 0)
      term1_inv <- solve(term1)
      Sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
      hh <- h_deriv(tcf.thet, tcf.bet, length(rhoEst$coeff))
      ans <- hh %*% Sig %*% t(hh)
    }
    else{
      bst.fi <- function(dt, inds, formula, cp){
        dat <- dt[inds,]
        out <- rhoMLogit(formula, data = dat)
        ROCs.tcf(method = "fi", T = dat[,1],
                 Dvec = as.matrix(dat[, c("D1", "D2", "D3")]),
                 rhoEst = out, cps = cp)
      }
      name.var <- as.character(attributes(terms(rhoEst$formula))$variables)[-1]
      data <- data.frame(T, rhoEst$D, rhoEst$X[,name.var[-1]], Dvec)
      names(data) <- c("Temp", name.var, "D1", "D2", "D3")
      if(parallel){
        res.bst <- boot(data, bst.fi, R = nR, formula = rhoEst$formula, cp = cp,
                        parallel = "snow", ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.fi, R = nR, formula = rhoEst$formula, cp = cp)
      ans <- var(res.bst$t)
    }
  }
  else if(method == "msi"){
    if(!BOOT){
      term1 <- EstFuncIE_deriv(V, T, rhoEst$X, cp, rhoEst$coeff, rhoEst$values,
                               rhoEst$Hess, tcf.thet, tcf.bet, m = 1)
      term2 <- EstFuncIE(Dvec.temp, V, T, rhoEst$X, cp, rhoEst$coeff,
                         rhoEst$values, tcf.thet, tcf.bet, m = 1)
      term1_inv <- solve(term1)
      Sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
      hh <- h_deriv(tcf.thet, tcf.bet, length(rhoEst$coeff))
      ans <- hh %*% Sig %*% t(hh)
    }
    else{
      bst.msi <- function(dt, inds, formula, cp){
        dat <- dt[inds,]
        out <- rhoMLogit(formula, data = dat)
        ROCs.tcf(method = "msi", T = dat[,1],
                 Dvec = as.matrix(dat[, c("D1", "D2", "D3")]), V = dat[, "V"],
                 rhoEst = out, cps = cp)
      }
      name.var <- as.character(attributes(terms(rhoEst$formula))$variables)[-1]
      data <- data.frame(T, rhoEst$D, rhoEst$X[,name.var[-1]], V, Dvec)
      names(data) <- c("Temp", name.var, "V", "D1", "D2", "D3")
      if(parallel){
        res.bst <- boot(data, bst.msi, R = nR, formula = rhoEst$formula, cp = cp,
                        parallel = "snow", ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.msi, R = nR, formula = rhoEst$formula, cp = cp)
      ans <- var(res.bst$t)
    }
  }
  else if(method == "ipw"){
    if(!BOOT){
      term1 <- EstFuncIPW_deriv(Dvec.temp, V, T, piEst$X, cp, piEst$coeff,
                                piEst$values, piEst$Hess, tcf.thet, tcf.bet,
                                piEst$model)
      term2 <- EstFuncIPW(Dvec.temp, V, T, piEst$X, cp, piEst$coeff,
                          piEst$values, tcf.thet, tcf.bet, piEst$model)
      term1_inv <- solve(term1)
      Sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
      hh <- h_deriv(tcf.thet, tcf.bet, length(piEst$coeff))
      ans <- hh %*% Sig %*% t(hh)
    }
    else{
      bst.ipw <- function(dt, inds, formula, cp){
        dat <- dt[inds,]
        out <- psglm(formula, data = dat, test = FALSE, trace = FALSE)
        ROCs.tcf(method = "ipw", T = dat[,1],
                 Dvec = as.matrix(dat[, c("D1", "D2", "D3")]), V = dat[, 2],
                 piEst = out, cps = cp)
      }
      name.var <- as.character(attributes(terms(piEst$formula))$variables)[-1]
      data <- data.frame(T, V, piEst$X[,name.var[-1]], Dvec)
      names(data) <- c("Temp", name.var, "D1", "D2", "D3")
      if(parallel){
        res.bst <- boot(data, bst.ipw, R = nR, formula = piEst$formula, cp = cp,
                        parallel = "snow", ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.ipw, R = nR, formula = piEst$formula, cp = cp)
      ans <- var(res.bst$t)
    }
  }
  else{
    if(!BOOT){
      term1 <- EstFuncSPE_deriv(Dvec.temp, V, T, rhoEst$X, piEst$X, cp,
                                rhoEst$coeff, rhoEst$values, rhoEst$Hess,
                                piEst$coeff, piEst$values, piEst$Hess, tcf.thet,
                                tcf.bet, piEst$model)
      term2 <- EstFuncSPE(Dvec.temp, V, T, rhoEst$X, piEst$X, cp, rhoEst$coeff,
                          rhoEst$values, piEst$coeff, piEst$values, tcf.thet,
                          tcf.bet, piEst$model)
      term1_inv <- solve(term1)
      Sig <- term1_inv %*% (t(term2) %*% term2) %*% t(term1_inv)
      hh <- h_deriv(tcf.thet, tcf.bet, length(rhoEst$coeff) + length(piEst$coeff))
      ans <- hh %*% Sig %*% t(hh)
    }
    else{
      bst.spe <- function(dt, inds, formula.rho, formula.pi, cp){
        dat <- dt[inds,]
        out.rho <- rhoMLogit(formula.rho, data = dat)
        out.pi <- psglm(formula.pi, data = dat, test = FALSE, trace = FALSE)
        ROCs.tcf(method = "spe", T = dat[,1],
                 Dvec = as.matrix(dat[, c("D1", "D2", "D3")]), V = dat[, "V"],
                 rhoEst = out.rho, piEst = out.pi, cps = cp)
      }
      name.var.rho <- as.character(attributes(terms(rhoEst$formula))$variables)[-1]
      name.var.pi <- as.character(attributes(terms(piEst$formula))$variables)[-1]
      data.rho <- data.frame(rhoEst$X[,name.var.rho[-1]])
      data.pi <- data.frame(piEst$X[,name.var.pi[-1]])
      data <- data.frame(T, V, rhoEst$D, union(data.rho, data.pi), Dvec)
      names(data) <- c("Temp", name.var.pi[1], name.var.rho[1],
                       union(name.var.rho[-1], name.var.pi[-1]), "D1", "D2", "D3")
      if(parallel){
        res.bst <- boot(data, bst.spe, R = nR, formula.rho = rhoEst$formula,
                        formula.pi = piEst$formula, cp = cp, parallel = "snow",
                        ncpus = ncpus)
      }
      else res.bst <- boot(data, bst.spe, R = nR, formula.rho = rhoEst$formula,
                           formula.pi = piEst$formula, cp = cp)
      ans <- var(res.bst$t)
    }
  }
  rownames(ans) <- colnames(ans) <- paste("TCF", c(1:3), sep = "")
  return(ans)
}

