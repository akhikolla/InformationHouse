

#' Running Means
#'
#' Extracts the running means by sliding a window comprised of \code{wsize} time points, and in each window, the mean for each variable is computed.
#' Each time the window is slid, the oldest time point is discarded and the latest time point is added.
#'
#' @param data \emph{N} x \emph{v} dataframe where \emph{N} is the no. of time points and \emph{v} the no. of variables
#' @param wsize Window size
#'
#' @return Running means time series
#' @importFrom roll roll_mean
#' @export
#'
#' @examples
#' phase1=cbind(rnorm(50,0,1),rnorm(50,0,1)) #phase1: Means=0
#' phase2=cbind(rnorm(50,1,1),rnorm(50,1,1)) #phase2: Means=1
#' X=rbind(phase1,phase2)
#' RS=runMean(data=X,wsize=25)
#' ts.plot(RS, gpars=list(xlab="Window", ylab="Means", col=1:2,lwd=2))
#'
runMean <- function(data, wsize = 25) {
  data <- as.data.frame(data)
  labs <- colnames(data)
  data <- as.matrix(data) #the roll package needs a matrix

  N <- nrow(data)
  wstep <- 1 #set to 1 for the kcpRS package
  windows <- floor(((N - wsize) / wstep) + 1)

  RunMean_TS <- roll_mean(data, width = wsize)  #uses an online algorithm to compute the running statistics
  RunMean_TS <- RunMean_TS[seq(wsize, wsize + windows - 1, wstep), ]

  RunMean_TS <- as.data.frame(RunMean_TS)
  colnames(RunMean_TS) <- labs

  return(RunMean_TS)
}
