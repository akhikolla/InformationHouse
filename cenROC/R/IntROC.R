#' Time-dependent ROC curve estimation for interval-censored survival data
#'
#' @description  This function computes the time-dependent ROC curve for interval censored survival data using the cumulative sensitivity and dynamic specificity definitions.
#'  The ROC curves can be either empirical (non-smoothed) or smoothed with/without boundary correction. It also calculates the time-dependent AUC.
#' @usage IntROC(L, R, M, t, U = NULL, method = "emp", method2 = "pa", dist = "weibull",
#'         bw = NULL, ktype = "normal", len = 151, B = 0, alpha = 0.05, plot = "TRUE")
#' @param L The numericvector of left limit of observed time. For left censored observations \code{L == 0}.
#' @param R The numericvector of right limit of observed time. For right censored observation \code{R == inf}.
#' @param M The numeric vector of marker values.
#' @param U The numeric vector of cutoff values.
#' @param t A scaler time point used to calculate the ROC curve.
#' @param len The length of the grid points for ROC estimation. Default is \code{151}.
#' @param method The method of ROC curve estimation. The possible options are "\code{emp}" empirical metod; "\code{untra}" smooth without boundary correction and "\code{tra}" is smooth ROC curve estimation with boundary correction. The default is the "\code{emp}" empirical method.
#' @param method2 A character indication type of modeling. This include nonparametric \code{"np"}, parmetric \code{"pa"} and semiparametric \code{"sp"}. The default is the "\code{np}" parametric with weibull distribution.
#' @param dist A character incating the type of distribution for parametric model. This includes are \code{"exponential"}, \code{"weibull"}, \code{"gamma"}, \code{"lnorm"}, \code{"loglogistic"} and \code{"generalgamma"}.
#' @param ktype A character string giving the type kernel distribution to be used for smoothing the ROC curve: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". By default, the "\code{normal}" kernel is used.
#' @param bw A character string specifying the bandwidth estimation method. The possible options are "\code{NR}" for the normal reference, the plug-in "\code{PI}" and the cross-validation "\code{CV}". The default is the "\code{NR}" normal reference method. It is also possible to use a numeric value.
#' @param B The number of bootstrap samples to be used for variance estimation. The default is \code{0}, no variance estimation.
#' @param alpha The significance level. The default is \code{0.05}.
#' @param plot The logigal parameter to see the ROC curve plot. Default is \code{TRUE}.
#' @details This function implments time-dependent ROC curve and the corresponding AUC using the model-band and nonparametric for the estimation of conditional survival function. The empirical (non-smoothed) ROC estimate and the smoothed ROC estimate with/without boundary correction can be obtained using this function.
#' The smoothed ROC curve estimators require selecting a bandwidth parametr for smoothing the ROC curve. To this end, three data-driven methods: the normal reference "\code{NR}", the plug-in "\code{PI}" and the cross-validation "\code{CV}" were implemented.
#' See Beyene and El Ghouch (2020) for details.
#' @return Returns the following items:
#' @return    \code{ROC     } The vector of estimated ROC values. These will be numeric numbers between zero
#' @return    \code{        } and one.
#' @return    \code{U       } The vector of grid points used.
#' @return    \code{AUC      } A data frame of dimension \eqn{1 \times 4}. The columns are: AUC, standard error of AUC, the lower
#' @return    \code{         }               and upper limits of bootstrap CI.
#' @return    \code{bw       } The computed value of bandwidth. For the empirical method this is always \code{NA}.
#' @return    \code{Dt      } The vector of estimated event status.
#' @return    \code{M       } The vector of Marker values.
#' @useDynLib cenROC
#' @importFrom stats approx dnorm pnorm qnorm dunif quantile
#' @import survival
#' @importFrom icenReg ic_par ic_sp getFitEsts
#' @importFrom Rcpp evalCpp
#' @examples library(cenROC)
#'
#' data(hds)
#'
#' IntROC(L=hds$L, R=hds$R, M=hds$M, t=2)$AUC
#'
#' @references Beyene, K. M. and El Ghouch A. (2020). Time-dependent ROC curves estimator for interval-censored survival data.
#' @export

IntROC <-  function(L, R, M, t, U = NULL, method = "emp", method2 = "pa", dist = "weibull", bw = NULL, ktype = "normal", len = 151, B = 0, alpha = 0.05, plot = "TRUE"){
  if(is.null(U)) {U <- seq(0, 1, length.out = len) }
  if(!is.vector(L, mode = "numeric") | !is.vector(M, mode = "numeric")| !is.vector(t,  mode = "numeric"))
    stop(paste0("Error! all numeric vectors L, M, and t should be specified"))

    Dt <- ICsur(L = L, R = R, M = M, t = t, method = method2, dist = dist)$positive;
    estim <- RocFun(U = U, D = Dt, M = M, method = method, bw = bw);
    ROC <- estim$roc;
    AUCc <- 1 - estim$auc;
    AUC <- data.frame(AUC = round(AUCc, 4) , sd = NA, LCL = NA, UCL = NA)
    if (B > 0){
      data <- data.frame(L = L, R = R, M = M)
      aucb <-NULL
      rocb <- matrix(NA, nrow = length(U), ncol = B)
      for (i in 1:B){
        bootsample <- sample(1:nrow(data), nrow(data), replace=TRUE)
        dat <- data[bootsample, ]
        D <- ICsur(L=dat$L, R=dat$R, M=dat$M, t=t, method=method2, dist=dist)$positive
        estim1 <- RocFun(U = U, D = D, M = dat$M, method = method, bw = bw, ktype = ktype)
        aucb[i] <- 1 - estim1$auc
        rocb[, i] <- estim1$roc
      }
      SP <- unname(quantile(aucb, p=c(alpha/2, 1-alpha/2), na.rm = TRUE))
      AUC <- data.frame(AUC = round(AUCc, 4) , sd = round(sd(aucb, na.rm = TRUE), 4), LCL = round(SP[1], 4), UCL = round(SP[2], 4))
      qroc <- unname(apply(rocb, 1, quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE))
      ROC <- data.frame(ROC = ROC, LCL = qroc[1, ], UCL = qroc[2, ])
    }


    if (plot == "TRUE") {
      if (B==0){
        plot(c(0, U, 1), c(0, (ROC), 1), type = "l", lwd = 2, col.lab = "blue", col = "blue",
             xlab = "False positive rate", ylab = "True positive rate")
        abline(coef = c(0,1))
      }else{
        plot(c(0, U, 1), c(0, (ROC[, 1]), 1), type = "l", lwd = 2, col.lab = "blue", col = "blue",
             xlab = "False positive rate", ylab = "True positive rate")
        # abline(coef = c(0,1))
        polygon(c(rev(c(0,U,1)), c(0,U,1)), c(rev(c(0,ROC[, 3],1)), c(0,ROC[, 2],1)), col = 'darkslategray1', border = NA)
        lines(c(0,U,1), c(0, (ROC[, 1]), 1), lty = 'solid', col = 'blue', lwd = 2)
        lines(c(0,U,1), c(0, (ROC[, 3]), 1), lty = 'dashed', col = 'red', lwd = 2)
        lines(c(0,U,1), c(0, (ROC[, 2]), 1), lty = 'dashed', col = 'red', lwd = 2)
      }
    }
  return(list(AUC = AUC, ROC = ROC, U = U, Dt = Dt, M = M, bw = estim$bw))
}

