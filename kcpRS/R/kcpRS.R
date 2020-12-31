


#' KCP on the running statistics
#'
#' Given a user-specified function \code{RS_fun} to compute the running statistics (see \code{\link{runMean}}, \code{\link{runVar}}, \code{\link{runAR}} and \code{\link{runCorr}}), a KCP permutation test (see \code{\link{permTest}}) is first implemented to test whether
#' there is at least one significant change point, then through model selection most optimal number of change points is chosen.
#' @aliases kcpRS kcpRS.default plot.kcpRS print.kcpRS summary.kcpRS
#' @param data data \emph{N} x \emph{v} dataframe where \emph{N} is the number of time points and \emph{v} the number of variables
#' @param RS_fun Running statistics function: Should require \code{wsize} and \code{wstep} as input and return a dataframe of running statistics
#' as output. The rows of this dataframe should correspond to the windows and the columns should correspond to the variable(s) on which the running statistics were computed.
#' @param RS_name Name of the monitored running statistic.
#' @param wsize Window size
#' @param nperm Number of permutations used in the permutation test
#' @param Kmax Maximum number of change points desired
#' @param alpha Significance level of the permutation test
#' @param varTest If set to FALSE, only the variance DROP test is implemented, and if set to TRUE, both the variance test and the variance DROP tests are implemented.
#' @param ncpu number of cpu cores to use
#'
#' @return \item{RS_name}{Name indicated for the monitored running statistic}
#' \item{RS}{Dataframe of running statistics with rows corresponding to the time window and columns corresponding to
#' the (combination of) variable(s) on which the running statistics were computed}
#' \item{wsize}{Selected window size}
#' \item{varTest}{Selected choice of implementation for varTest}
#' \item{nperm}{Selected number of permutations}
#' \item{alpha}{Selected significance level of the permutation test}
#' \item{subTest_alpha}{Significance level of each subtest. If \code{varTest}=0, \code{subTest_alpha} is equal to \code{alpha} since only the variance drop test is implemented.
#' If \code{varTest}=1, \code{subTest_alpha}=\code{alpha}/2 since two subtests are carried out and Bonferonni correction is applied.}
#' \item{BestK}{Optimal number of change points}
#' \item{changePoints}{Change point location(s)}
#' \item{p_var_test}{P-value of the variance test}
#' \item{p_varDrop_test}{P-value of the variance drop test}
#' \item{CPs_given_K}{A matrix comprised of the minimized variance criterion \emph{Rmin} and the optimal change point location(s) for each \emph{k} from 1 to \code{Kmax}
#' }
#' @references \cite{Cabrieto, J., Tuerlinckx, F., Kuppens, P., Wilhelm, F., Liedlgruber, M., & Ceulemans, E. (2018). Capturing correlation changes by applying kernel change point
#' detection on the running correlations. Information Sciences, 447, 117-139.} \url{https://doi.org/10.1016/j.ins.2018.03.010}
#'
#' \cite{Cabrieto, J., Adolf, J., Tuerlinckx, F., Kuppens, P., & Ceulemans, E. (2018). Detecting long-lived autodependency changes in a multivariate system via change point detection
#' and regime switching models. Scientific Reports, 8, 15637, 1-15.} \url{https://doi.org/10.1038/s41598-018-33819-8}
#'
#' @export
#' @importFrom stats cov
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster  stopCluster detectCores
#'
#' @examples
#' phase1=cbind(rnorm(50,0,1),rnorm(50,0,1)) #phase1: Means=0
#' phase2=cbind(rnorm(50,1,1),rnorm(50,1,1)) #phase2: Means=1
#' X=rbind(phase1,phase2)
#' res=kcpRS(data=X,RS_fun=runMean,RS_name="Mean",wsize=25,
#' nperm=1000,Kmax=10,alpha=.05,varTest=FALSE,ncpu=1)
#'
#' summary(res)
#' plot(res)
#'
kcpRS <-
  function(data,
           RS_fun,
           RS_name,
           wsize = 25,
           nperm = 1000,
           Kmax = 10,
           alpha = .05,
           varTest = FALSE,
           ncpu = 1) {
    UseMethod("kcpRS")
  }
