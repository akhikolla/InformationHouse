#' KCP (Kernel Change Point) Detection
#'
#' Finds the most optimal change point(s) in the running statistic time series \code{RunStat} by
#' looking at their kernel-based pairwise similarities.
#'
#' @param RunStat  Dataframe of running statistics with rows corresponding to the windows and the
#' columns corresponding to the variable(s) on which these running statistics were computed.
#' @param Kmax Maximum number of change points
#' @param wsize Window size
#'
#' @return \item{kcpSoln}{A matrix comprised of the minimized variance criterion \emph{Rmin} and the optimal change point location(s) for each \emph{k} from 1 to \code{Kmax}}
#' @importFrom stats dist median
#' @importFrom Rcpp evalCpp
#' @useDynLib kcpRS
#' @references \cite{Arlot, S., Celisse, A., & Harchaoui, Z. (2012). Kernel change-point detection.} \url{http://arxiv.org/abs/1202.3878}
#'
#' \cite{Cabrieto, J., Tuerlinckx, F., Kuppens, P., Grassmann, M., & Ceulemans, E. (2017). Detecting correlation changes in multivariate time series:
#' A comparison of four non-parametric change point detection methods. Behavior Research Methods, 49, 988-1005.}
#' \url{https://doi.org/10.3758/s13428-016-0754-9}


kcpa <- function(RunStat,
                 Kmax = 10,
                 wsize = 25) {
  locate_cp <- function(H,ncp,wsize) {
      #to find change point locations given K no. of change points

      #H<- h
      #K<- k
      #ws<- wsize
      #wst<-wstep

      wstep <- 1 #set to 1 for the kcpRS package
      d <- ncp + 1
      cps <- ncol(H)
      cp <- cps

      while (d > 0) {
        cp = H[d, cp]
        cps <- c(cp, cps)
        d <- d - 1
      }

      cps = cps + 1 #cps in the window scale

      #convert to the real observation scale
      if (wsize %% 2 > 0) {
        cps = ceiling(wsize / 2) + (cps - 1) * wstep
      }  #conversion: "ceiling(ws/2)" is the midpoint of the first window (if window size is odd)
      #"(cps-1)*wst" is what we add to reach the change point window
      #essentially, we get the midpoint of the change point window

      if (wsize %% 2 == 0) {
        cps = ceiling((wsize / 2 + .5) + (cps - 1) * wstep)
      } #"ws/2 + .5" is the midpoint of the first window (if window size is even)
      # we round up to have a whole number as a change point: thus existing (next new observation is the change point)

      ncps <- length(cps) - 1
      cps <- cps[2:(ncps)]
      return(cps)
    }



  ###################################################
  #Derive the II and H matrices for the KCP solution
  ###################################################
  X_ = as.matrix(RunStat)
  N <- nrow(X_)
  D <- Kmax + 1            #no. of segments
  II_ <- matrix(Inf, D, N)  #variance at a time point
  H_ <- matrix(0, D, N)     #change point

  kcp_res <- getScatterMatrix(II_, X_, H_) #the c function that generates the kcp soln
  II <- kcp_res[[1]]      #minimized variance criterion (unstandardized)
  H <- kcp_res[[2]]       #candidate change point locations
  Rmin <- II[, N]/N       #standardized variance criterion

  ##################################################
  #find change point locations from k=1 to Kmax
  ##################################################

  optimal_CPs <- matrix(0, nrow = D, ncol = Kmax)

  for (ncp in 1:Kmax) {
    optimal_CPs[(ncp + 1), (1:ncp)] = locate_cp(H, ncp, wsize)
  } #gives the location of the change point given ncp=no. of change points

  kcpSoln <- cbind(0:Kmax, round(Rmin, 4), optimal_CPs)
  kcpSoln <- as.data.frame(kcpSoln, row.names = FALSE)

  colnames(kcpSoln) <- c("k", "Rmin", paste0("CP", 1:Kmax))

  return(kcpSoln)
}
