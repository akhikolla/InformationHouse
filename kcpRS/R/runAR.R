

#' Running Autocorrelations
#'
#' Extracts the running autocorrelations by sliding a window comprised of \code{wsize} time points, and in each window, the autocorrelation for each variable is computed.
#' Each time the window is slid, the oldest time point is discarded and the latest time point is added.
#'
#' @param data \emph{N} x \emph{v} dataframe where \emph{N} is the no. of time points and \emph{v} the no. of variables
#' @param wsize Window size
#'
#' @return Running autocorrelations time series
#' @importFrom roll roll_cor
#' @export
#'
#' @examples
#' phase1=cbind(rnorm(50,0,1),rnorm(50,0,1)) #phase1: AutoCorr=0
#' phase2=cbind(rnorm(50,0,1),rnorm(50,0,1))
#' phase2=filter(phase2,.50, method="recursive") #phase2: AutoCorr=.5
#' X=rbind(phase1,phase2)
#' RS=runAR(data=X,wsize=25)
#' ts.plot(RS, gpars=list(xlab="Window", ylab="Autocorrelation", col=1:2,lwd=2))


runAR <- function(data, wsize = 25) {
  data <- as.data.frame(data)
  labs <- colnames(data)
  data <- as.matrix(data) #the roll package needs a matrix

  N <- nrow(data)
  wstep <- 1 #set to 1 for the kcpRS package
  windows <- floor(((N - 1 - wsize) / wstep) + 1)

  X <- as.matrix(data[2:N, ])
  X_lagged <- as.matrix(data[1:(N - 1), ])

  RunAR_TS <- NULL

  rs_temp <- roll_cor(X_lagged, as.matrix(X), width = wsize)  #uses the roll package to compute the running statistics
  rs_temp <- rs_temp[, , seq(wsize, wsize + windows - 1, wstep)]

  if (ncol(X) > 1) {
    ind_mat = diag(rep(TRUE, ncol(data)))
    rs_temp = rs_temp[ind_mat]
  }

  RunAR_TS <- as.data.frame(matrix(rs_temp, nrow = windows, byrow = TRUE))

  colnames(RunAR_TS) <- labs
  return(RunAR_TS)
}
