#' Mayo Marker Data
#'
#' Two marker values with event time and censoring status for the subjects in Mayo PBC data.
#' @docType data
#' @usage data(mayo)
#' @keywords datasets
#' @format A data frame with 312 observations and 4 variables: time (event time/censoring time), censor (censoring
#'        indicator), mayoscore4, mayoscore5. The two scores are derived from 4 and 5 covariates
#'        respectively.
#' @references Heagerty, P. J., and Zheng, Y. (2005). Survival model predictive accuracy and ROC curves. \emph{Biometrics}, 61(1), 92-105.
"mayo"


#' NASA Hypobaric Decompression Sickness Marker Data
#'
#' This data contains the marker values with the left and right limits of the observed time for the subjects in NASA Hypobaric Decompression Sickness Data.
#' 
#' @docType data
#' @usage data(hds)
#' @keywords datasets
#' @format This is a data frame with 238 observations and 3 variables: L (left limit of the observed time), R (right limit of the observed time) 
#'         and M (marker). The marker is a score derived by combining the covariates Age, Sex, TR360, and Noadyn.
#' @references Beyene, K. M. and El Ghouch A. (2020). Time-dependent ROC curves estimator for interval-censored survival data. 
"hds"