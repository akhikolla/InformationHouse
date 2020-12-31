#' Running Correlations
#'
#' Extracts the running correlations by sliding a window comprised of \code{wsize} time points, and in each window,
#' the correlation of each pair of variables is computed.
#' Each time the window is slid, the oldest time point is discarded and the latest time point is added.
#'
#' @param data \emph{N} x \emph{v} dataframe where \emph{N} is the no. of time points and \emph{v} the no. of variables
#' @param wsize window size
#'
#' @return Running correlationa time series
#' @export
#' @importFrom roll roll_cor
#' @importFrom utils combn
#' @examples
#' data(MentalLoad)
#' RS<-runCorr(data=MentalLoad,wsize=25)
#' ts.plot(RS, gpars=list(xlab="Window", ylab="Correlations", col=1:3,lwd=2))

runCorr <- function(data, wsize = 25) {
  data = as.data.frame(data)
  v <- ncol(data)
  N <- nrow(data)
  labs <- colnames(data)
  data <- as.matrix(data) #the roll package needs a matrix
  wstep <- 1 #set to 1 for the kcpRS package

  if (v > 1) {
    windows <- floor(((N - wsize) / wstep) + 1)
    pairs <- combn(v, 2)
    npairs <- dim(pairs)[2]

    rs_temp <- roll_cor(data,
               width = wsize,
               center = TRUE,
               scale = TRUE)  #uses the roll to compute the running statistics

    rs_temp <- rs_temp[, , seq(wsize, wsize + windows - 1, wstep)]
    rs_temp <- rs_temp[lower.tri(rs_temp[, , 1])]
    rs_temp <- .5 * (log(1 + rs_temp) - log(1 - rs_temp)) #convert to fishers_z

    RunCorr_TS <- matrix(rs_temp, nrow = windows, byrow = TRUE)
    RunCorr_TS <- as.data.frame(RunCorr_TS)

    colnames(RunCorr_TS) <-
      paste0(labs[pairs[1, ]], "&", labs[pairs[2, ]])
  }

  if (v == 1) {
    RunCorr_TS <- NULL
    stop("Cannot compute correlations for a univariate time series.")
  }

  return(RunCorr_TS)
}
