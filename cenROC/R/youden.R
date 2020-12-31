#' Computes optimal cutoff point using the Youden index criteria
#'
#' @description  This function computes the optimal cutoff point using the Youden index criteria of both right and interval censored time-to-event data.
#'               The Youden index estimator can be either empirical (non-smoothed) or smoothed with/without boundary correction.
#' @param est The object returned either by \code{cenROC} or \code{IntROC}.
#' @param plot The logical parameter to see the ROC curve plot along with the Youden inex. The default is \code{TRUE}.
#' @details In medical decision-making, obtaining the optimal cutoff value is crucial to identify subject at high risk
#'  of experiencing the event of interest.  Therefore, it is necessary to select a marker value that classifies subjects into
#'  healthy and diseased groups. To this end, in the literature, several methods for selecting optimal cutoff point have been
#'  proposed. In this package, we only included the Youden index criteria.
#'
#' @references Youden, W.J. (1950). Index for rating diagnostic tests. \emph{Cancer} 3, 32–35.
#' @examples library(cenROC)
#'
#' data(mayo)
#' resu <- cenROC(Y=mayo$time, M=mayo$mayoscore5, censor=mayo$censor, t=365*6, plot="FALSE")
#' youden(resu,  plot="TRUE")
#' @export

youden <- function(est,  plot="FALSE"){
  M <- est$M
  D <- est$Dt
  if(is.na(est$bw)){
    ord <- order(M)
    M <- M[ord]
    D <- D[ord]
    #################################################################
    ########### Sensitivity and Specificity #########################
    sens <- spec <- NULL ;
    for (m in M) {
      sens <- c( sens, sum(D * as.numeric(M > m)) / sum(D) ) ;
      spec <- c( spec, sum((1 - D) * as.numeric(M <= m)) / sum(1 - D) ) ;
    }
    ################################################################
    ################# YODEN INDEX ##################################
    Jm <- sens + spec - 1;
    opt.Jm   <- max(Jm) ;
    opt.sens <- unique(sens[Jm == opt.Jm]);
    opt.spec <- unique(spec[Jm == opt.Jm]);
    opt.cut  <- unique(M[Jm == opt.Jm]);
    out <- data.frame(Youden.index = opt.Jm, cutopt = opt.cut, sens = opt.sens, spec = opt.spec)
    if(plot=="TRUE"){
      plot( c(1, 1 - spec, 0), c(1, sens, 0),  type="l", lwd=3, lty=1, col.lab="blue", col="blue", xlab="False positive rate", ylab="True positive rate")
      abline(0,1)
      segments(1-opt.spec, 1-opt.spec, 1-opt.spec, opt.sens, lwd=2, lty=1, col = "red")
      text(x=1 - opt.spec + 0.05, y=opt.sens - 0.1, "J", col = "red", font = 2, cex=1.2,srt=90)
      legend("bottomright", title=expression('J(cut'['opt']*')'), legend=paste(round(opt.Jm, 3), " (", round(opt.cut, 3), ")", sep=""),  cex=0.95)

    }
  }
  else
  {
    ord <- order(M)
    M <- M[ord]
    D <- D[ord]
    #################################################################
    ########### Sensitivity and Specificity #########################
    sens <- spec <- NULL ;
    for (m in M) {
      sens <- c( sens, sum(D * as.numeric(M > m)) / sum(D) ) ;
    }
    Jm <- est$ROC - est$U;
    opt.Jm   <- max(Jm) ;
    opt.sens <- est$ROC[Jm == opt.Jm];
    opt.spec <- 1 - est$U[Jm == opt.Jm];
    opt.cut <-  withCallingHandlers(approx(sens, M, xout = opt.sens)$y, warning=function(w){invokeRestart("muffleWarning")})
    out <- data.frame(Youden.index = opt.Jm, cutopt = opt.cut, sens = opt.sens, spec = opt.spec)
    if(plot=="TRUE"){
      plot( c(0, est$U, 1), c(0, (est$ROC), 1),  type="l", lwd=3, lty=1, col.lab="blue", col="blue", xlab="False positive rate", ylab="True positive rate")
      abline(0,1)
      segments(1-opt.spec, 1-opt.spec, 1-opt.spec, opt.sens, lwd=2, lty=1, col = "red")
      text(x=1-opt.spec + 0.05, y=opt.sens - 0.1, "J", col = "red", font = 2, cex=1.2,srt=90)
      legend("bottomright", title=expression('J(cut'['opt']*')'), legend=paste(round(opt.Jm, 3), " (", round(opt.cut, 3), ")", sep=""), cex=0.95)

    }
  }
  return(out)
}
