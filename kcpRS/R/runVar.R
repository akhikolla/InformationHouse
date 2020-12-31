

#' Running Variances
#'
#' Extracts the running variances by sliding a window comprised of \code{wsize} time points, and in each window, the variance for each variable is computed.
#' Each time the window is slid, the oldest time point is discarded and the latest time point is added.
#'
#' @param data \emph{N} x \emph{v} dataframe where \emph{N} is the no. of time points and \emph{v} the no. of variables
#' @param wsize Window size
#'
#' @return Running variances time series
#' @importFrom roll roll_var
#' @export
#'
#' @examples
#' phase1=cbind(rnorm(50,0,1),rnorm(50,0,1)) #phase1: SD=1
#' phase2=cbind(rnorm(50,0,2),rnorm(50,0,2)) #phase2: SD=2
#' X=rbind(phase1,phase2)
#' RS=runVar(data=X,wsize=25)
#' ts.plot(RS, gpars=list(xlab="Window", ylab="Variances", col=1:2,lwd=2))
#'
runVar <- function(data, wsize = 25) {
  data <- as.data.frame(data)
  labs <- colnames(data)
  data <- as.matrix(data) #the roll package needs a matrix

  N <- nrow(data)
  wstep <- 1 #set to 1 for the kcpRS package
  windows <- floor(((N - wsize) / wstep) + 1)

  runVar_TS <- roll_var(data, width = wsize)  #uses the roll package to compute the running statistics
  runVar_TS <- runVar_TS[seq(wsize, wsize + windows - 1, wstep), ]

  runVar_TS <- as.data.frame(runVar_TS)
  colnames(runVar_TS) <- colnames(data)

  return(runVar_TS)
}
