####==========================================================================####
## This file consists of some functions that are related to the computation VUS ##
## by using the VUS formula.                                                    ##
## Date: 04/07/2016																															##
####==========================================================================####
##
#' @title Estimation methods for volume under ROC surface (VUS)
#'
#' @description \code{vus} computes bias-corrected estimates of the volume under the ROC surface for evaluating the accuracy of a continuous diagnostic test.
#'
#' @param method  name of bias-corrected estimation method to be used for estimating the VUS in presence of verification bias. See \code{\link{ROCs}} for more details.
#' @param T  a numeric vector containing the diagnostic test values. \code{NA} values are not admitted.
#' @param Dvec  a n * 3  binary matrix with the three columns, corresponding to three classes of the disease status. In row i, 1 in column j indicates that the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param V  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param rhoEst  a result of a call to \code{\link{rhoMLogit}} of \code{\link{rhoKNN}} to fit the disease model.
#' @param piEst  a result of a call to \code{\link{psglm}} to fit the verification model.
#' @param ci a logical value. If TRUE (default), computes an confidence interval of VUS and tests the null hypothesis H0: VUS = 1/6.
#' @param ci.level  an confidence level to be used for constructing the confidence interval; default 0.95.
#' @param BOOT a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a bootstrap resampling is employed to estimate the asymptotic variance of bias-corrected VUS estimates. See \code{\link{asyVarVUS}}.
#' @param nR  the number of bootstrap replicates, which is used for FULL or KNN estimator, or option \code{BOOT = TRUE}. Usually this will be a single positive integer.
#' @param parallel  a logical value. If \code{TRUE}, a parallel computing is employed to the bootstrap resampling process.
#' @param ncpus  number of processes to be used in parallel computing. Default is a half of available cores.
#' @param trace a logical value. If \code{TRUE}, tracing information on the progress of the estimation is produced.
#'
#' @details
#' The function implements five bias-corrected estimation methods in To Duc et al. (2016, 2018) for estimating VUS of a three-class continuous diagnostic test in presence of verification bias. The estimators are full imputation (FI), mean score imputation (MSI), inverse probability weighted (IPW), semiparametric efficient (SPE) and K nearest-neighbor (KNN), see \code{\link{ROCs}}. These esitmators work under MAR assumption.
#'
#' The standard error of the estimates are obtained through the function \code{\link{asyVarVUS}}. In particular, the standard error of the FULL estimate is computed by bootrap resampling method or by Jackknife approach proposed in Guangming et al. (2013). For the bias-corrected estimates, the standard errors are computed by using asymptotic theory (with respect to FI, MSI, IPW and SPE estimator) or bootstrap resampling method (with respect to KNN estimator). A confidence interval of VUS also is given. A logit transformation is also applied for obtaining the confidence interval.
#'
#' The default value of the number of bootstrap replicates is 250.
#'
#' Note that, before apply the functions \code{vus}, the use of \code{\link{preDATA}} might be needed to check the monotone ordering disease classes and to create the matrix format for disease status.
#'
#' @return \code{vus} returns an object of class inheriting from "vus" class.
#'
#' The function \code{\link{print.vus}} can be used to print a summary of the results.
#'
#' An object of class "vus" is a list containing at least the following components:
#'
#' \item{vus.fit}{the estimate of VUS.}
#' \item{std}{the standard error, obtained by using asymptotic theory or bootstrap resampling method.}
#' \item{call}{the matched call.}
#' \item{t.stat}{t-statistic.}
#' \item{p.val_norm}{p-value correspond to normal-test.}
#' \item{ci.norm}{the confidence interval of VUS by using normal approximation.}
#' \item{ci.logit}{the confidence interval of VUS via logit transform.}
#' \item{ci.level}{the confidence level used.}
#' \item{BOOT}{the value of \code{BOOT}.}
#' \item{nR}{the number of bootstrap replicates used.}
#'
#' In addition, the name of method used to estimate VUS also is given as the attribute of \code{vus.fit}.
#'
#' @references
#' To Duc, K., Chiogna, M. and Adimari, G. (2018)
#' Nonparametric estimation of ROC surfaces in presence of verification bias.
#' \emph{REVSTAT Statistical Journal}. Accepted.
#'
#' To Duc, K., Chiogna, M. and Adimari, G. (2016)
#' Bias-corrected methods for estimating the receiver operating characteristic surface of continuous diagnostic tests.
#' \emph{Electronic Journal of Statistics}, \bold{10}, 3063-3113.
#'
#' Guangming, P., Xiping, W. and Wang, Z. (2013)
#' Non-parameteric statistical inference for $P(X < Y < Z)$.
#' \emph{Sankhya A}, \bold{75}, 1, 118-138.
#'
#' @examples
#' data(EOC)
#' head(EOC)
#'
#'
#' \dontrun{
#' # FULL data estimator
#' Dfull <- preDATA(EOC$D.full, EOC$CA125)
#' Dvec.full <- Dfull$Dvec
#' vus("full", T = EOC$CA125, Dvec = Dvec.full)
#' }
#'
#' # Preparing the missing disease status
#' Dna <- preDATA(EOC$D, EOC$CA125)
#' Dfact.na <- Dna$D
#' Dvec.na <- Dna$Dvec
#' # FI estimator
#' rho.out <- rhoMLogit(Dfact.na ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#'
#' vus("fi", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.out)
#'
#' \dontrun{
#' # MSI estimator
#' vus("msi", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.out)
#'
#' # IPW estimator
#' pi.out <- psglm(V ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#' vus("ipw", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, piEst = pi.out)
#'
#' # SPE estimator
#' vus("spe", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.out, piEst = pi.out)
#'
#' # KNN estimator, K = 1, Mahalanobis distance
#' XX <- cbind(EOC$CA125, EOC$CA153, EOC$Age)
#' rho.maha.1nn <- rhoKNN(X = XX, Dvec = Dvec.na, V = EOC$V, K = 1, type = "mahala")
#' vus("knn", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.maha.1nn)
#' }
#'
#' @import utils
#' @import parallel
#' @importFrom Rcpp evalCpp
#' @useDynLib bcROCsurface, .registration = TRUE
#' @export
vus <- function(method = "full", T, Dvec, V, rhoEst = NULL, piEst = NULL,
                ci = TRUE, ci.level = ifelse(ci, 0.95, NULL), BOOT = FALSE,
                nR = ifelse(ci, 250, NULL), parallel = FALSE,
                ncpus = ifelse(parallel, detectCores()/2, NULL), trace = TRUE){
  ## return the function call
  call <- match.call()
  ## checking the method
  methodtemp <- substitute(me, list(me = method))
  okMethod <- c("full", "fi", "msi", "ipw", "spe", "knn")
  if(length(methodtemp) > 1) stop(gettextf("Please, choose one method from %s", paste(sQuote(okMethod), collapse = ", ")), domain = NA)
  if (!is.character(methodtemp)) methodtemp <- deparse(methodtemp)
  if (!is.element(methodtemp, okMethod)){
    stop(gettextf("the required method \"%s\" is not available; it should be one of %s", methodtemp, paste(sQuote(okMethod), collapse = ", ")),
         domain = NA)
  }
  ## checking the argument T
  if(missing(T)) stop("argument \"T\" is missing \n")
  if(!inherits(T, "numeric") | any(is.na(T))) stop("variable \"T\" must be a numeric vector and not include NA values")
  name_diagnostic <- substitute(T)
  if (!is.character(name_diagnostic))  {
    name_diagnostic <- deparse(name_diagnostic)
  }
  name_diagnostic <- unlist(strsplit(name_diagnostic, NULL))
  if(any(name_diagnostic %in% c("$"))){
    id.name <- which(name_diagnostic %in% c("$"))
    name_diagnostic <- paste(name_diagnostic[(id.name[1] + 1) : length(name_diagnostic)],
                             collapse = "")
  }
  else name_diagnostic <- paste(name_diagnostic, collapse = "")
  method_name <- toupper(method)
  ## checking Dvec
  if(missing(Dvec)) stop("argument \"Dvec\" is missing \n")
  if(!inherits(Dvec, "matrix") | ncol(Dvec) != 3 | !all(is.element(na.omit(Dvec), c(0,1)))) stop("variable \"Dvec\" must be a binary matrix with 3 columns")
  if(length(T) != nrow(Dvec)) stop(gettextf("arguments imply differing number of observation: %d", length(T)), gettextf(", %d", nrow(Dvec)), domain = NA)
  Dvec.flag <- any(is.na(Dvec))
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
  ## checking V
  if(missing(V)){
    if(methodtemp != "full" | Dvec.flag) stop("argument \"V\" is missing, in addition, the method is not \"full\" or argument \"Dvec\" includes NA")
    if(trace){
      cat("Hmm, look likes the full data.\n")
      cat("The verification status is not available.\n")
      cat("You are working on FULL or Complete Case approach.\n")
      cat("Number of observation:", length(T), "\n")
      cat("The diagnostic test:", name_diagnostic, "\n")
      cat("Processing .... \n")
      flush.console()
    }
    V <- NULL
  }
  else{
    if(all(V == 0) | !all(is.element(V, c(0,1))) ) stop("There are mistakes in \"V\". Please, check your input and see whether it is correct or not")
    if(nrow(Dvec) != length(V) | length(T) != length(V)) stop(gettextf("arguments imply differing number of observation: %d", nrow(Dvec)), gettextf(", %d", length(T)), gettextf(", %d", length(V)), domain = NA)
    if(all(V == 1)){
      if(methodtemp != "full" | Dvec.flag) stop("Please, check your inputs and see whether they are correct or not.\n If you want to estimate Complete Case approach, please, remove the missing values in the \n data set and try again with the option of \"full\" method.")
      if(trace){
        cat("Hmm, look likes the full data\n")
        cat("Number of observation:", length(T), "\n")
        cat("All subjects underwent the verification process\n")
        cat("You are working on FULL or Complete Case approach\n")
        cat("The diagnostic test:", name_diagnostic, "\n")
        cat("Processing .... \n")
        flush.console()
      }
    }
    else{
      rv <- mean(V)
      if(!Dvec.flag){
        if(trace){
          cat("Warning: There are no NA values in variable Dvec, while",
              paste(round(rv*100), "%", sep = ""), "of the subjects receive disease verification. \n")
          cat("BE CAREFULL OF YOUR INPUT AND RESULTS \n")
          cat("Number of observation:", length(T), "\n")
          cat("You required estimate VUS using", method_name, "approach \n")
          cat("The diagnostic test:", name_diagnostic, "\n")
          cat("Processing .... \n")
          flush.console()
        }
      }
      else{
        if(trace){
          cat("Hmm, look likes the incomplete data\n")
          cat("Number of observation:", length(T), "\n")
          cat(paste(round(rv*100), "%", sep = ""), "of the subjects receive disease verification. \n")
          cat("You required estimate VUS using", method_name, "approach \n")
          cat("The diagnostic test:", name_diagnostic, "\n")
          cat("Processing .... \n")
          flush.console()
        }
      }
    }
  }
  ## Main body
  if(methodtemp == "full"){
    if(Dvec.flag){
      cat("Opp! Look likes wrong method \n")
      ques <- readline("Do you want use Complete Case (CC) approach? [y/n]: ")
      if(ques %in% c("y", "n")){
        if(ques == "y"){
          if(trace){
            cat("Number of observation:", length(T), "\n")
            cat("We are estimating VUS by using CC method \n")
            cat("BE CAREFULL OF THE RESULTS. THIS CAN MAKE THE DISTORTED INFERENCE IN VUS \n")
            cat("Processing .... \n")
            flush.console()
          }
          Dvec <- Dvec[V == 1,]
          T <- T[V == 1]
          ans <- vusC(T, Dvec)
        }
        else if(ques == "n") stop("Sorry! The FULL method can not compute with NA value(s) of disease status Dvec.")
      }
      else stop("The answer was wrong. Please, choose y for Yes, or n for No!")
    }
    ans <- vusC(T, Dvec)
  }
  else if(methodtemp == "fi"){
    if(is.null(rhoEst)) stop("The input of argument \"rhoEst\" is needed for ", method_name, " estimator")
    if(length(T) != nrow(rhoEst$values)) stop(gettextf("arguments imply differing number of observation: %d", length(T), ", %d", nrow(rhoEst$values)), domain = NA)
    ans <- vusC(T, rhoEst$values)
  }
  else if(methodtemp %in% c("msi", "knn")){
    if(is.null(rhoEst)) stop("argument \"rhoEst\" is needed for", method_name, "estimator")
    Dvectemp <- Dvec
    if(Dvec.flag) Dvectemp[is.na(Dvec)] <- 99
    Dmsi <- Dvectemp*V + (1 - V)*rhoEst$values
    ans <- vusC(T, Dmsi)
  }
  else if(methodtemp == "ipw"){
    if(is.null(piEst)) stop("argument \"piEst\" is needed for", method_name, "estimator")
    Dvectemp <- Dvec
    if(Dvec.flag) Dvectemp[is.na(Dvec)] <- 99
    Dipw <- Dvectemp*V/piEst$values
    ans <- vusC(T, Dipw)
  }
  else{
    if(is.null(rhoEst) | is.null(piEst)) stop("arguments \"rhoEst\" and \"pi.est\" are needed for", method_name, "estimator")
    Dvectemp <- Dvec
    if(Dvec.flag) Dvectemp[is.na(Dvec)] <- 99
    Dspe <- Dvectemp*V/piEst$values - (V/piEst$values - 1)*rhoEst$values
    ans <- vusC(T, Dspe)
  }
  attr(ans, "name") <- method_name
  res <- list(vus.fit = ans, call = call)
  class(res) <- "vus"
  if(ci){
    var.ans <- asyVarVUS(res, T = T, Dvec = Dvec, V = V, rhoEst = rhoEst,
                         piEst = piEst, BOOT = BOOT, nR = nR, parallel = parallel,
                         ncpus = ncpus)
    cf.ans <- ans + c(-1, 1)*qnorm((1 + ci.level)/2)*sqrt(var.ans)
    sd.log <- sqrt(var.ans)/(ans*(1 - ans))
    log.cf <- log(ans/(1 - ans)) + c(-1, 1)*qnorm((1 + ci.level)/2)*sd.log
    cf.ans.tran <- exp(log.cf)/(1 + exp(log.cf))
    W <- (ans - 1/6)/sqrt(var.ans)
    p.val_norm <- 1 - pnorm(W)
    res <- list(vus.fit = ans, std = sqrt(var.ans), t.stat = W,
                p.val_norm = p.val_norm, ci.norm = cf.ans, ci.logit = cf.ans.tran,
                call = call, ci.level = ci.level, BOOT = BOOT, nR = nR)
		class(res) <- "vus"
  }
  if(trace) cat("DONE\n")
  res
}


## The function print.vus
#' @title Print summary results of VUS
#'
#' @description \code{print.vus} prints the results for the output of function \code{\link{vus}}.
#'
#' @method print vus
#' @param x an object of class "vus", a result of a call to \code{\link{vus}}.
#' @param digits minimal number of significant digits, see \code{\link{print.default}}.
#' @param ... further arguments passed to \code{\link{print}} method.
#'
#' @details \code{print.vus} shows a nice format of the summary table for the VUS estimate results. Some information on the diagnostic test, the fitted values of VUS, and confidence intervals are shown.
#'
#' @seealso \code{\link{vus}}
#'
#' @export
print.vus <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\n")
  cat("CALL: ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n \n", sep = "")
  if(!is.null(x$std)){
    res.tab <- c(x$ci.level*100, x$ci.norm, x$ci.logit)
    res.tab <- format(round(res.tab, digits = digits))
    res.tab[1L] <- paste("\n",x$ci.level*100,"% ",sep = "")
    res.tab[2*(1L:2L)] <- paste(" (",res.tab[2*(1L:2L)],",",sep = "")
    res.tab[2*(1L:2L) + 1L] <- paste(res.tab[2*(1L:2L) + 1L],") ")
    p.val <- x$p.val_norm
    stat <- x$t.stat
    test.tab <- cbind(stat, p.val)
    colnames(test.tab) <- c("Test Statistic", "P-value")
    rownames(test.tab) <- "Normal-test"
    ci.name <- c("       Normal        ", "       Logit        ")
    cat("Estimate of VUS:", format(round(x$vus.fit, digits = digits)), "\n")
    cat("Standard error:", format(round(x$std, digits = digits)), "\n")
    cat("\nIntervals:")
    cat("\nLevel", ci.name)
    cat(res.tab)
    if(attr(x$vus, "name") %in% c("FULL", "KNN")){
      if(x$BOOT) cat("\nEstimation of Standard Error and Intervals are based on Bootstrap with", x$nR, "replicates\n")
      else cat("\nEstimation of Standard Error and Intervals are based on Jackknife approach \n")
    }
    else{
      cat("\nEstimation of Standard Error and Intervals are based on Asymptotic Theory \n")
    }
    cat("\n")
    cat("Testing the null hypothesis H0: VUS = 1/6 \n")
    printCoefmat(test.tab, has.Pvalue = TRUE)
  }
  else{
    cat("Estimate of VUS:", round(x$vus.fit, digits = digits), "\n")
  }
  invisible(x)
}
