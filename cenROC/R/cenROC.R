#' Estimation of the time-dependent ROC curve for right censored survival data
#'
#' @description  This function computes the time-dependent ROC curve for right censored survival data using the cumulative sensitivity and dynamic specificity definitions.
#'  The ROC curves can be either empirical (non-smoothed) or smoothed with/wtihout boundary correction. It also calculates the time-dependent area under the ROC curve (AUC).
#' @usage cenROC(Y, M, censor, t, U = NULL, h = NULL, bw = "NR", method = "tra",
#'     ktype = "normal", ktype1 = "normal", B = 0, alpha = 0.05, plot = "TRUE")
#' @param Y The numeric vector of event-times or observed times.
#' @param M The numeric vector of marker values for which the time-dependent ROC curves is computed.
#' @param censor The censoring indicator, \code{1} if event, \code{0} otherwise.
#' @param t A scaler time point at which the time-dependent ROC curve is computed.
#' @param U The vector of grid points where the ROC curve is estimated. The default is a sequence of \code{151} numbers between \code{0} and \code{1}.
#' @param bw A character string specifying the bandwidth estimation method for the ROC itself. The possible options are "\code{NR}" for the normal reference, the plug-in "\code{PI}" and the cross-validation "\code{CV}". The default is the "\code{NR}" normal reference method. The user can also introduce a numerical value.
#' @param h A scaler for the bandwidth of Beran's weight calculaions. The defualt is the value obtained by using the method of Sheather and Jones (1991).
#' @param method The method of ROC curve estimation. The possible options are "\code{emp}" emperical metod; "\code{untra}" smooth without boundary correction and "\code{tra}" is smooth ROC curve estimation with boundary correction. The default is the "\code{tra}" smooth ROC curve estimate with boundary correction.
#' @param ktype A character string giving the type kernel distribution to be used for smoothing the ROC curve: "\code{normal}", "\code{epanechnikov}", "\code{biweight}", or "\code{triweight}". By default, the "\code{normal}" kernel is used.
#' @param ktype1 A character string specifying the desired kernel needed for Beran weight calculation. The possible options are "\code{normal}", "\code{epanechnikov}", "\code{tricube}", "\code{boxcar}", "\code{triangular}", or "\code{quartic}". The defaults is "\code{normal}" kernel density.
#' @param B The number of bootstrap samples to be used for variance estimation. The default is \code{0}, no variance estimation.
#' @param alpha The significance level. The default is \code{0.05}.
#' @param plot The logical parameter to see the ROC curve plot. The default is \code{TRUE}.
#' @details The empirical (non-smoothed) ROC estimate and the smoothed ROC estimate with/without boundary correction can be obtained using this function.
#' The smoothed ROC curve estimators require selecting two bandwidth parametrs: one for Beran’s weight calculation and one for smoothing the ROC curve.
#' For the latter, three data-driven methods: the normal reference "\code{NR}", the plug-in "\code{PI}" and the cross-validation "\code{CV}" were implemented.
#' To select the bandwidth parameter needed for Beran’s weight calculation, by default, the plug-in method of Sheather and Jones (1991) is used but it is also possible introduce a numeric value.
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
#' @importFrom stats pnorm qnorm quantile approx bw.SJ integrate sd
#' @importFrom condSURV Beran
#' @importFrom graphics abline legend segments text lines polygon
#' @examples library(cenROC)
#'
#' data(mayo)
#' cenROC(Y=mayo$time, M=mayo$mayoscore5, censor=mayo$censor, t=365*6)$AUC
#'
#' @author
#' Kassu Mehari Beyene, Catholic University of Louvain. \code{<kasu.beyene@uclouvain.be>}
#'
#' Anouar El Ghouch, Catholic University of Louvain. \code{<anouar.elghouch@uclouvain.be>}
#' @references Beyene, K. M. and El Ghouch A. (2020). Smoothed time-dependent ROC curves for right-censored survival data. \emph{submitted}.
#' @references Sheather, S. J. and Jones, M. C. (1991). A Reliable data-based bandwidth selection method for kernel density estimation. \emph{Journal of the Royal Statistical Society}. Series B (Methodological) 53(3): 683–690.
#' @export

cenROC <- function(Y, M, censor, t, U = NULL, h = NULL, bw = "NR",  method = "tra", ktype = "normal", ktype1 = "normal", B = 0, alpha = 0.05, plot = "TRUE")
  {
    if (is.null(U)) {U <- seq(0, 1, length.out = 151)}
    if (!is.vector(Y, mode = "numeric") |
        !is.vector(M, mode = "numeric") |
        !is.vector(censor, mode = "numeric"))
      print("Error! all numeric vectors Y, M and censor should be specified")
    else{
      Dt <- Csurv(Y = Y, M = M, censor = censor, t = t, h = h, kernel = ktype1)$positive
      estim <- RocFun(U = U, D = Dt, M = M, method = method, bw = bw, ktype = ktype)
      ROC <- estim$roc
      AUCc <- 1 - estim$auc
      AUC <- data.frame(AUC = round(AUCc, 4) , sd = NA, LCL = NA, UCL = NA)
    }

  if (B > 0){
    data <- data.frame(Y=Y, M=M, cen=censor)
    aucb <- NULL
    rocb <- matrix(NA, nrow = length(U), ncol = B)
    for (i in 1:B){
      bootsample <- sample(1:nrow(data), nrow(data), replace=TRUE)
      dat <- data[bootsample, ]
      Dt <- Csurv(Y = dat$Y, M = dat$M, censor = dat$cen, t = t, h = h, kernel = ktype1)$positive
      estim <- RocFun(U = U, D = Dt, M = dat$M, method = method, bw = bw, ktype = ktype)
      aucb[i] <- 1 - estim$auc
      rocb[, i] <- estim$roc
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
