####==========================================================================####
## The R code for the true class fractions.                                     ##
##                                                                              ##
####==========================================================================####
#' @name ROCsurface
#'
#' @title Receiver operating characteristics surface for a continuous diagnostic test
#' @aliases ROCs
#' @aliases ROCs.tcf
#'
#'
#' @description
#' \code{ROCs.tcf} is used to obtain bias-corrected estimates of the true class fractions (TCFs) for evaluating the accuracy of a continuous diagnostic test for a given cut point \eqn{(c_1, c_2)}, with \eqn{c_1 < c_2}.
#'
#' \code{ROCs} provides bias-corrected estimates of the ROC surfaces of the continuous diagnostic test by using TCF.
#'
#'
#' @param method  a estimation method to be used for estimating the true class fractions in presence of verification bias. See 'Details'.
#' @param T  a numeric vector containing the diagnostic test values. \code{NA} values are not allowed.
#' @param Dvec  a n * 3  binary matrix with the three columns, corresponding to three classes of the disease status. In row i, 1 in column j indicates that the i-th subject belongs to class j, with j = 1, 2, 3. A row of \code{NA} values indicates a non-verified subject.
#' @param V  a binary vector containing the verification status (1 verified, 0 not verified).
#' @param A  a vector/matrix of dimension n * q containing the values of the covariate(s). If the method is \code{"knn"} and \code{ellipsoid = TRUE}, \code{A} is needed to compute the asymptotic covariance of TCFs at a fixed cut point. The default \code{NULL} is suitable for the remaining methods.
#' @param rhoEst  a result of a call to \code{\link{rhoMLogit}} of \code{\link{rhoKNN}} to fit the disease model.
#' @param piEst  a result of a call to \code{\link{psglm}} to fit the verification model.
#' @param ncp  the dimension of cut point grid. It is used to determine the cut points (see 'Details'). Default 100.
#' @param plot  if \code{TRUE}(the default), a 3D plot of ROC surface is produced.
#' @param ellipsoid  a logical value. If TRUE, adds an ellipsoidal confidence region for TCFs at a specified cut point to current plot of ROC surface.
#' @param cpst  a specified cut point, which used to construct the ellipsoid confidence region. If \code{m} ellipsoid confidence regions are required, \code{cpst} must be matrix with \code{m} rows and 2 columns. Default \code{NULL}.
#' @param level  an confidence level to be used for constructing the ellipsoid confidence region; default 0.95.
#' @param sur.col  color to be used for plotting ROC surface and ellipsoid.
#' @param ...  optional arguments to be passed to \code{\link[rgl]{plot3d}}, \code{\link[rgl]{surface3d}}.
#' @param cps  a cut point \eqn{(c_1, c_2)}, with \eqn{c_1 < c_2}, which used to estimate TCFs. If \code{m} estimates of TCFs are required, \code{cps} must be matrix with \code{m} rows and 2 columns.
#' @param BOOT a logical value. Default = \code{FALSE}. If set to \code{TRUE}, a bootstrap resampling is employed to estimate the asymptotic variance-covariance matrix of TCFs at the cut point \code{cpst}. See more details in \code{\link{asyCovTCF}}.
#' @param nR  the number of bootstrap replicates, which is used for FULL estimator, or option \code{BOOT = TRUE}. Usually this will be a single positive integer. Default 250.
#' @param parallel  a logical value. If \code{TRUE}, a parallel computing is employed to the bootstrap resampling process.
#' @param ncpus  number of processes to be used in parallel computing. Default is half of of available cores.
#'
#'
#' @details In a three-class diagnostic problem, quantities used to evaluate the accuracy of a diagnostic test are the true class fractions (TCFs). For a given pair of cut points \eqn{(c_1, c_2)} such that \eqn{c_1 < c_2}, subjects are classified into class 1 (\eqn{D_1}) if \eqn{T < c_1}; class 2 (\eqn{D_2}) if \eqn{c_1 \le T < c_2}; class 3 (\eqn{D_3}) otherwise. The true class fractions of the test \eqn{T} at \eqn{(c_1, c_2)} are defined as
#' \deqn{TCF_1(c_1) = P(T < c_1| D_1 = 1) = 1 - P(T \ge c_1| D_1 = 1),}
#' \deqn{TCF_2(c_1, c_2) = P(c_1 \le T < c_2| D_2 = 1) = P(T \ge c_1| D_2 = 1) - P(T \ge c_2| D_2 = 1),}
#' \deqn{TCF_3(c_2) = P(T > c_2| D_3 = 1) = P(T \ge c_2| D_3 = 1). }
#'
#' The receiver operating characteristic (ROC) surface is the plot of \eqn{TCF_1}, \eqn{TCF_2} and \eqn{TCF_3} by varying the cut point \eqn{(c_1, c_2)} in the domain of the diagnostic test. The cut points \eqn{(c_1, c_2)} are produced by designing a cut point grid with \code{ncp} dimension. In this grid, the points satisfying \eqn{c_1 < c_2} are selected as the cut points. The number of the cut points are obtained as \eqn{ncp(ncp - 1)/2}, for example, the default is 4950.
#'
#' These functions implement the bias-corrected estimators in To Duc et al (2016a, 2016b) for estimating TCF of a three-class continuous diagnostic test in presence of verification bias. The estimators work under MAR assumption. Five methods are provided, namely:
#' \itemize{
#'    \item Full imputation (FI): uses the fitted values of the disease model to replace the true disease status (both of missing and non-missing values).
#'    \item Mean score imputation (MSI): replaces only the missing values by the fitted values of the disease model.
#'    \item Inverse probability weighted (IPW): weights each observation in the verification sample by the inverse of the sampling fraction (i.e. the probability that the subject was selected for verification).
#'    \item Semiparametric efficient (SPE): replaces the true disease status by the double robust estimates.
#'    \item K nearest-neighbor (KNN): uses K nearest-neighbor imputation to obtain the missing values of the true disease status.
#' }
#'
#' The argument \code{method} must be selected from the collection of the bias-corrected methods, i.e., \code{"full"}, \code{"fi"}, \code{"msi"}, \code{"ipw"}, \code{"spe"} and \code{"knn"}.
#'
#' The ellipsoidal confidence region of TCFs at a given cut point can be constructed by using a normal approximation and plotted in the ROC surface space. The confidence level (default) is 0.95.
#'
#' Note that, before using the functions \code{ROCs} and \code{ROCs.tcf}, the use of \code{\link{preDATA}} might be needed to check the monotone ordering disease classes and to create the matrix format for disease status.
#'
#' @return \code{ROCs} returns a list, with the following components:
#' \item{vals}{the estimates of TCFs at all cut points.}
#' \item{cpoint}{the cut points are used to construct the ROC surface.}
#' \item{ncp}{dimension of the cut point grid.}
#' \item{cpst}{the cut points are used to construct the ellipsoidal confidence regions.}
#' \item{tcf}{the estimates of TCFs at the cut points \code{cpst}.}
#' \item{message}{an integer code or vector. 1 indicates the ellipsoidal confidence region is available.}
#'
#' \code{ROCs.tcf} returns a vector having estimates of TCFs at a cut point when \code{cps} is a vector with two elements, or a list of estimates of TCFs at \code{m} cut points when \code{cps} is a \code{m*2} matrix. In addition, some attributes called \code{theta}, \code{beta}, \code{cp} and \code{name} are given. Here, \code{theta} is a probability vector, with 3 element, corresponding to the disease prevalence rates of three classes. \code{beta} is also a probability vector having 4 components, which are used to compute TCFs, see To Duc el al. (2016a, 2016b) for more details. \code{cp} is the specified cut point that is used to estimate TCFs. \code{name} indicates the method used to estimate TCFs. These attributes are required to compute the asymptotic variance-covariance matrix of TCFs at the given cut point.
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
#' @seealso \code{\link{psglm}}, \code{\link{rhoMLogit}}, \code{\link[rgl]{plot3d}}.
#'
#' @examples
#' data(EOC)
#' head(EOC)
#'
#' \dontrun{
#' # FULL data estimator
#' Dfull <- preDATA(EOC$D.full, EOC$CA125)
#' Dvec.full <- Dfull$Dvec
#' ROCs("full", T = EOC$CA125, Dvec = Dvec.full, , ncp = 30, ellipsoid = TRUE,
#'      cpst = c(-0.56, 2.31))
#' }
#'
#' # Preparing the missing disease status
#' Dna <- preDATA(EOC$D, EOC$CA125)
#' Dvec.na <- Dna$Dvec
#' Dfact.na <- Dna$D
#'
#' # FI estimator
#' rho.out <- rhoMLogit(Dfact.na ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#' ROCs("fi", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.out, ncp = 30)
#'
#' \dontrun{
#' # Plot ROC surface and add ellipsoid confidence region
#' ROCs("fi", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.out, ncp = 30,
#'      ellipsoid = TRUE, cpst = c(-0.56, 2.31))
#'
#' # MSI estimator
#' ROCs("msi", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.out, ncp = 30,
#'      ellipsoid = TRUE, cpst = c(-0.56, 2.31))
#'
#' # IPW estimator
#' pi.out <- psglm(V ~ CA125 + CA153 + Age, data = EOC, test = TRUE)
#' ROCs("ipw", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, piEst = pi.out, ncp = 30,
#'      ellipsoid = TRUE, cpst = c(-0.56, 2.31))
#'
#' # SPE estimator
#' ROCs("spe", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.out, ncp = 30,
#'      piEst = pi.out, ellipsoid = TRUE, cpst = c(-0.56, 2.31))
#'
#' # 1NN estimator
#' XX <- cbind(EOC$CA125, EOC$CA153, EOC$Age)
#' K.opt <- CVknn(X = XX, Dvec = Dvec.na, V = EOC$V, type = "mahala", plot = TRUE)
#' rho.1nn <- rhoKNN(X = XX, Dvec = Dvec.na, V = EOC$V, K = K.opt, type = "mahala")
#' ROCs("knn", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.1nn, ncp = 30,
#'      ellipsoid = TRUE, cpst = c(-0.56, 2.31))
#'
#' ## Compute TCFs at three cut points
#' cutps <- rbind(c(0,0.5), c(0,1), c(0.5,1))
#' ROCs.tcf("spe", T = EOC$CA125, Dvec = Dvec.na, V = EOC$V, rhoEst = rho.out, ncp = 30,
#'          piEst = pi.out, cps = cutps)
#'}
#'

##
#' @rdname ROCs
#' @export
ROCs.tcf <- function(method = "full", T, Dvec, V = NULL, rhoEst = NULL, piEst = NULL, cps){
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
  method_name <- toupper(methodtemp)
  ## checking Dvec
  if(missing(Dvec)) stop("argument \"Dvec\" is missing \n")
  if(!inherits(Dvec, "matrix") | ncol(Dvec) != 3 | !all(is.element(na.omit(Dvec), c(0,1)))) stop("variable \"Dvec\" must be a binary matrix with 3 columns")
  if(length(T) != nrow(Dvec)) stop(gettextf("arguments imply differing number of observation: %d", length(T)), gettextf(", %d", nrow(Dvec)), domain = NA)
  ##
  if(!is.null(rhoEst)){
    if(!is.element(class(rhoEst),  c("prob_dise", "prob_dise_knn")))
      stop("\"rhoEst\" not a result of rhoMlogit or rhoKNN \n")
    if(is.element(class(rhoEst),  c("prob_dise_knn"))){
      if(is.null(rhoEst$K)) stop("The \"rhoEst\" is a list and contains results of rho corresponding to many K! \n Please, choose one result of them \n")
    }
  }
  if(!is.null(piEst)){
    if(!is.element(class(piEst),  c("prob_veri")))
      stop("\"piEst\" not a result of psglm \n")
  }
  ## Main body
  tcfCal <- function(tt, dd, cp, weight, methodd){
    thet.hat <- colSums(dd)/weight
    names(thet.hat) <- c()
    beta.func <- function(tt, dd, ctp, weight){
      bet <- numeric(4)
      bet[1] <- sum((tt >= ctp[1])*dd[,1])/weight
      bet[2] <- sum((tt >= ctp[1])*dd[,2])/weight
      bet[3] <- sum((tt >= ctp[2])*dd[,2])/weight
      bet[4] <- sum((tt >= ctp[2])*dd[,3])/weight
      return(bet)
    }
    if(!is.matrix(cp)){
      bet.hat <- beta.func(tt = tt, dd = dd, ctp = cp, weight = weight)
      res <- numeric(3)
      res[1] <- 1 - bet.hat[1]/thet.hat[1]
      res[2] <- (bet.hat[2] - bet.hat[3])/thet.hat[2]
      res[3] <- bet.hat[4]/thet.hat[3]
      attr(res, "beta") <- bet.hat
      attr(res, "theta") <- thet.hat
      attr(res, "cp") <- cp
      attr(res, "name") <- methodd
      names(res) <- c("TCF1", "TCF2", "TCF3")
      class(res) <- "tcfs"
    }
    else{
      res <- list()
      bet.hat <- t(apply(cp, 1, beta.func, tt = tt, dd = dd, weight = weight))
      temp <- matrix(0, nrow = nrow(cp), ncol = 3)
      temp[,1] <- 1 - bet.hat[,1]/thet.hat[1]
      temp[,2] <- (bet.hat[,2] - bet.hat[,3])/thet.hat[2]
      temp[,3] <- bet.hat[,4]/thet.hat[3]
      colnames(temp) <- c("TCF1", "TCF2", "TCF3")
      for(i in 1:nrow(cp)){
        temp1 <- temp[i,]
        attr(temp1, "beta") <- bet.hat[i,]
        attr(temp1, "theta") <- thet.hat
        attr(temp1, "cp") <- cp[i,]
        attr(temp1, "name") <- methodd
        res[[i]] <- temp1
        class(res[[i]]) <- "tcfs"
      }
    }
    return(res)
  }
  nsize <- length(T)
  if(method == "full"){
    if(any(is.na(Dvec))){
      stop("The", method_name, "method can not access the NA values of disease status")
    }
    else{
      ans <- tcfCal(T, Dvec, cps, weight = nsize, methodd = method_name)
    }
  }
  if(method == "fi"){
    if(is.null(rhoEst)) stop("argument \"rhoEst\" is needed for ", method_name, "estimator")
    else{
      ans <- tcfCal(T, rhoEst$values, cps, weight = nsize, methodd = method_name)
    }
  }
  if(method %in% c("msi", "knn")){
    if(is.null(rhoEst) | is.null(V)) stop("argument \"rhoEst\" is needed for ", method_name, "estimator")
    else{
      Dvectemp <- Dvec
      if(any(is.na(Dvec))) Dvectemp[is.na(Dvec)] <- 99
      Dmsi <- Dvectemp*V + (1 - V)*rhoEst$values
      ans <- tcfCal(T, Dmsi, cps, weight = nsize, methodd = method_name)
    }
  }
  if(method == "ipw"){
    if(is.null(piEst) | is.null(V)) stop("argument \"piEst\" is needed for ", method_name, "estimator")
    else{
      Dvectemp <- Dvec
      if(any(is.na(Dvec))) Dvectemp[is.na(Dvec)] <- 99
      Dipw <- Dvectemp*V/piEst$values
      ans <- tcfCal(T, Dipw, cps, weight = sum(V/piEst$values), methodd = method_name)
    }
  }
  if(method == "spe"){
    if(is.null(rhoEst) | is.null(V) | is.null(piEst)) stop("arguments \"rhoEst\" and \"piEst\" are needed for ", method_name, " estimator")
    else{
      Dvectemp <- Dvec
      if(any(is.na(Dvec))) Dvectemp[is.na(Dvec)] <- 99
      Dspe <- Dvectemp*V/piEst$values - (V/piEst$values - 1)*rhoEst$values
      ans <- tcfCal(T, Dspe, cps, weight = nsize, methodd = method_name)
    }
  }
  ans
}

## ROCs function
#' @rdname ROCs
#' @import utils
#' @import rgl
#' @import parallel
#' @export
ROCs <- function(method = "full", T, Dvec, V, A = NULL, rhoEst = NULL, piEst = NULL,
                 ncp = 100, plot = TRUE, ellipsoid = FALSE, cpst = NULL,
                 level = 0.95, sur.col = c("gray40", "green"), BOOT = FALSE,
                 nR = 250, parallel = FALSE,
                 ncpus = ifelse(parallel, detectCores()/2, NULL), ...){
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
  if(class(T) != "numeric" | any(is.na(T))) stop("variable \"T\" must be a numeric vector and not include NA values")
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
  method_name <- toupper(methodtemp)
  ## checking Dvec
  if(missing(Dvec)) stop("argument \"Dvec\" is missing \n")
  if(isFALSE("matrix" %in% class(Dvec)) | ncol(Dvec) != 3 | !all(is.element(na.omit(Dvec), c(0,1)))) stop("variable \"Dvec\" must be a binary matrix with 3 columns")
  if(length(T) != nrow(Dvec)) stop(gettextf("arguments imply differing number of observation: %d", length(T)), gettextf(", %d", nrow(Dvec)), domain = NA)
  Dvec.flag <- any(is.na(Dvec))
  ##
  if(ellipsoid & is.null(cpst)) stop("argument \"cpst\" is required to construct the ellipsoidal confidence region of TCFs")
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
    cat("Hmm, look likes the full data\n")
    cat("Number of observation:", length(T), "\n")
    cat("The verification status is not available\n")
    cat("You are working on FULL or Complete Case approach\n")
    cat("The diagnostic test:", name_diagnostic, "\n")
    cat("Processing .... \n")
    flush.console()
    V <- NULL
  }
  else{
    if(all(V == 0) | !all(is.element(V, c(0,1))) ) stop("There are mistakes in \"V\". Please, check your input and see whether it is correct or not")
    if(nrow(Dvec) != length(V) | length(T) != length(V)) stop(gettextf("arguments imply differing number of observation: %d", length(T)), gettextf(", %d", nrow(Dvec)), domain = NA)
    if(all(V == 1)){
      if(methodtemp != "full" | Dvec.flag) stop("Please, check your inputs and see whether they are correct or not.\n If you want to estimate Complete Case approach, please, remove the missing values in the \n data set and try again with the option of \"full\" method.")
      cat("Hmm, look likes the full data\n")
      cat("Number of observation:", length(T), "\n")
      cat("All subjects underwent the verification process\n")
      cat("You are working on FULL or Complete Case approach\n")
      cat("The diagnostic test:", name_diagnostic, "\n")
      cat("Processing .... \n")
      flush.console()
    }
    else{
      rv <- mean(V)
      if(!Dvec.flag){
        cat("Warning: There are no NA values in variable Dvec, while", paste(round(rv*100), "%", sep = ""), "of the subjects receive disease verification. \n")
        cat("BE CAREFULL OF YOUR INPUT AND RESULTS \n")
        cat("Number of observation:", length(T), "\n")
        cat("You required estimate ROC surface using", method_name, "approach \n")
        cat("The diagnostic test:", name_diagnostic, "\n")
        cat("Processing .... \n")
        flush.console()
      }
      else{
        cat("Hmm, look likes the incomplete data\n")
        cat("Number of observation:", length(T), "\n")
        cat(paste(round(rv*100), "%", sep = ""), "of the subjects receive disease verification. \n")
        cat("You required estimate ROC surface using", method_name, "approach \n")
        cat("The diagnostic test:", name_diagnostic, "\n")
        if(ellipsoid){
          cat("The ellipsoidal confidence region(s) of TCFs are also constructed \n")
          if(!is.matrix(cpst) & !is.vector(cpst)) stop("The cut points should be a vector, which has two elements, or matrix, which has two columns.")
          if(is.vector(cpst)){
            if(length(cpst) != 2 | cpst[1] > cpst[2]) stop("The first cut point must be less than the second one.")
          }
          if(is.matrix(cpst)){
            if(ncol(cpst) != 2 | any(cpst[,1] > cpst[,2])) stop("the first column must be less than or equal to the second column.")
          }
        }
        cat("Processing .... \n")
        flush.console()
      }
    }
  }
  cp <- c(-Inf, seq(min(T), max(T), length.out = ncp - 2), Inf)
  cp1 <- rep(cp, seq(ncp - 1, 0, by = -1))
  cp2 <- c()
  for(i in 1:(ncp - 1)){
    cp2 <- c(cp2, cp[-c(1:i)])
  }
  cpoint <- cbind(cp1, cp2)
  ROCpoint <- ROCs.tcf(method = methodtemp, T = T, Dvec = Dvec, V = V,
                       rhoEst = rhoEst, piEst = piEst, cps = cpoint)
  ROCpoint <- matrix(unlist(ROCpoint), ncol = 3, byrow = TRUE)
  colnames(ROCpoint) <- c("TCF1", "TCF2", "TCF3")
  rownames(ROCpoint) <- paste("(",round(cp1, 3)," , " ,round(cp2, 3), ")", sep = "")
  ct1 <- numeric(ncp - 1)
  for(i in 1:(ncp - 1)){
    ct1[i] <- i*ncp - i*(i + 1)/2
  }
  tcf1 <- matrix(ROCpoint[ct1, 1], ncp - 1, ncp - 1, byrow = FALSE)
  tcf3 <- matrix(ROCpoint[1:(ncp - 1), 3], ncp - 1, ncp - 1, byrow = TRUE)
  tcf2 <- matrix(0, nrow = ncp - 1, ncol = ncp - 1)
  tcf2[lower.tri(tcf2, diag = TRUE)] <- ROCpoint[, 2]
  tcf2 <- t(tcf2)
  res <- list()
  res$vals <- ROCpoint
  res$cpoint <- cpoint
  res$ncp <- ncp
  if(plot){
    open3d()
    par3d(windowRect = 50 + c(0, 0, 640, 640))
    plot3d(0, 0, 0, type = "n", box = FALSE, xlab = " ", ylab = " ", zlab = " ",
           xlim = c(0,1), ylim = c(0,1), zlim = c(0,1), ...)
    surface3d(tcf1, tcf3, tcf2, col = sur.col[1], alpha = 0.5, ...)
    if(any(c("full","fi", "msi", "ipw", "spe") %in% methodtemp)){
      title_plot <- paste(toupper(methodtemp), "estimator")
    }
    if("knn" %in% methodtemp){
      title_plot <- paste(rhoEst$K, "NN-",
                          switch(rhoEst$type, "eucli" = "Euclidean",
                                 "manha" = "Manhattan", "canber" = "Canberra",
                                 "lagran" = "Lagrange", "mahala" = "Mahalanobis"),
                          sep = "")
    }
    title3d(main = title_plot, xlab = "TCF1",
            ylab = "TCF3", zlab = "TCF2", line = 3)
    play3d(spin3d(axis = c(0, 0, 1), rpm = 12.25), duration = 2)
    play3d(spin3d(axis = c(0, 1, 0), rpm = 0.3), duration = 2)
  }
  if(ellipsoid){
    shade.ellips <- function(orgi, sig, lev){
      t1 <- sig[2, 2]
      sig[2, 2] <- sig[3, 3]
      sig[3, 3] <- t1
      t1 <- sig[1, 2]
      sig[1, 2] <- sig[1, 3]
      sig[1, 3] <- t1
      sig[lower.tri(sig)] <- sig[upper.tri(sig)]
      ellips <- ellipse3d(sig, centre = orgi[c(1,3,2)], t = sqrt(qchisq(lev, 3)))
      return(ellips)
    }
    res$cpst <- cpst
    if(is.vector(cpst)){
      tcf.orgi <- ROCs.tcf(method = methodtemp, T = T, Dvec = Dvec, V = V,
                           rhoEst = rhoEst, piEst = piEst, cps = cpst)
      tcf.sig <- asyCovTCF(tcf.orgi, T = T, Dvec = Dvec, V = V,
                           rhoEst = rhoEst, piEst = piEst, BOOT = BOOT,
                           nR = nR, parallel = parallel, ncpus = ncpus)
      sig.test <- rcond(tcf.sig)
      res$tcf <- tcf.orgi
      if(sig.test < .Machine$double.eps){
        cat("The asymptotic variance-covariance matrix of TCFs at the cut point",
            paste("(",cpst[1],", ",cpst[2],")", sep = ""), "is not singular!\n")
        cat("The ellipsoidal confidence region is not available!\n")
        if(!BOOT & methodtemp != "full") cat("Try again with bootstrap process!\n")
        plot3d(tcf.orgi[1], tcf.orgi[3], tcf.orgi[2], type = "s", col = "red",
               radius = 0.01, add = TRUE, ...)
        res$message <- 0
      }
      else{
        res$message <- 1
        ellip.tcf <- shade.ellips(tcf.orgi, tcf.sig, level)
        plot3d(ellip.tcf, box = FALSE, col = sur.col[2], alpha = 0.5, xlim = c(0,1),
               ylim = c(0, 1), zlim = c(0, 1), xlab = " ", ylab = " ", zlab = " ",
               add = TRUE, ...)
        plot3d(tcf.orgi[1], tcf.orgi[3], tcf.orgi[2], type = "s", col = "red",
               radius = 0.01, add = TRUE, ...)
      }
    }
    if(is.matrix(cpst)){
      res$tcf <- matrix(0, nrow = nrow(cpst), ncol = 3,
                        dimnames = list(paste("(", round(cpst[,1], 3), " , " ,
                                              round(cpst[,2], 3), ")", sep = ""),
                                        c("TCF1", "TCF2", "TCF3")))
      res$message <- numeric(nrow(cpst))
      for(i in 1:nrow(cpst)){
        tcf.orgi <- ROCs.tcf(method = methodtemp, T = T, Dvec = Dvec, V = V,
                             rhoEst = rhoEst, piEst = piEst, cps = cpst[i,])
        tcf.sig <- asyCovTCF(tcf.orgi, T = T, Dvec = Dvec, V = V, rhoEst = rhoEst,
                             piEst = piEst, BOOT = BOOT, nR = nR,
                             parallel = parallel, ncpus = ncpus)
        sig.test <- rcond(tcf.sig)
        res$tcf[i,] <- tcf.orgi
        if(sig.test < .Machine$double.eps){
          cat("The asymptotic variance-covariance matrix of TCFs at the cut points",
              paste("(",cpst[i,1],", ",cpst[i,2],")", sep = ""), "is not singular!\n")
          cat("The ellipsoidal confidence region is not available!\n")
          if(!BOOT & methodtemp != "full") cat("Try again with bootstrap process!\n")
          plot3d(tcf.orgi[1], tcf.orgi[3], tcf.orgi[2], type = "s", col = "red",
                 radius = 0.01, add = TRUE, ...)
          res$message[i] <- 0
        }
        else{
          res$message[i] <- 1
          ellip.tcf <- shade.ellips(tcf.orgi, tcf.sig, level)
          plot3d(ellip.tcf, box = FALSE, col = sur.col[2], alpha = 0.5, xlim = c(0, 1),
                 ylim = c(0, 1), zlim = c(0, 1), xlab = " ", ylab = " ", zlab = " ",
                 add = TRUE, ...)
          plot3d(tcf.orgi[1], tcf.orgi[3], tcf.orgi[2], type = "s", col = "red",
                 radius = 0.01, add = TRUE, ...)
        }
      }
    }
  }
  cat("DONE\n")
  cat("===============================================================\n")
  selec.ord <- order(rowSums(ROCpoint) - 1, decreasing = TRUE)
  cat("Some values of TCFs:\n")
  print(ROCpoint[selec.ord[1:6],], 3)
  cat("\n")
  if(ellipsoid){
    cat("Some information for Ellipsoidal Confidence Region(s):\n")
    cat("Confidence level:", level,"\n")
    if(is.vector(cpst)){
      cat("TCFs at", paste("(", cpst[1], ", ", cpst[2], ")", sep = ""),"are:\n")
      temp <- tcf.orgi
      attributes(temp) <- NULL
      names(temp) <- c("TCF1", "TCF2", "TCF3")
      print(temp, 3)
    }
    if(is.matrix(cpst)){
      colnames(res$tcf) <- c("TCF1", "TCF2", "TCF3")
      rownames(res$tcf) <- paste("(",round(cpst[,1], 3)," , ",round(cpst[,2], 3),")",
                                 sep = "")
      cat("TCFs at", rownames(res$tcf), "are:\n")
      print(res$tcf, 3)
    }
  }
  cat("===============================================================\n")
  invisible(res)
}

