
#' KCP on the Running Statistics Workflow
#'
#' Any of the four basic running statistics (i.e., running means, running variances, running autocorrelations and running correlations) or a combination thereof can be scanned for change points.
#'
#' The workflow proceeds in two steps: First, the mean change points are flagged using KCP on the running means. If there are significant change points,
#' the data is centered based on the yielded change points. Otherwise, the data remains untouched for the next analysis. Second, the remaining running
#' statistics are computed using the centered data in the first step. The user can specify which running statistics to scan change points for
#' (see \code{RS_funs} and examples).

#' Bonferonni correction for tracking multiple running statistics can be implemented using the \code{bcorr} option.
#' @aliases kcpRS_workflow  kcpRS_workflow.default plot.kcpRS_workflow print.kcpRS_workflow summary.kcpRS_workflow
#' @param data data \emph{N} x \emph{v} dataframe where \emph{N} is the number of time points and \emph{v} the number of variables
#' @param RS_funs a vector of names of the functions that correspond to the running statistics to be monitored.
#' Options available: "runMean"=running mean, "runVar"=running variance, "runAR"=running autocorrelation and "runCorr"=running correlation.
#' @param wsize Window size
#' @param nperm Number of permutations used in the permutation test
#' @param Kmax Maximum number of change points desired
#' @param alpha Significance level for the full KCP-RS workflow analysis if \code{bcorr}=1. Otherwise, this is the significance level for each running statistic.
#' @param varTest If set to TRUE, only the variance DROP test is implemented, and if set to FALSE, both the variance test and the variance DROP tests are implemented.
#' @param bcorr Set to TRUE if Bonferonni correction is desired for the workflow analysis and set to FALSE otherwise.
#' @param ncpu number of cpu cores to use
#'
#' @return \item{kcpMean}{\code{kcpRS} solution for the running means. See output of \code{\link{kcpRS}} for further details.}
#' \item{kcpVar}{\code{kcpRS} solution for the running variances. See output of \code{\link{kcpRS}} for further details.}
#' \item{kcpAR}{\code{kcpRS} solution for the running autocorrelations. See output of \code{\link{kcpRS}} for further details.}
#' \item{kcpCorr}{\code{kcpRS} solution for the running correlations. See output of \code{\link{kcpRS}} for further details.}
#'
#' @references \cite{Cabrieto, J., Adolf, J., Tuerlinckx, F., Kuppens, P., & Ceulemans, E. (in press). An objective, comprehensive and flexible statistical framework for
#'  detecting early warning signs of mental health problems. Psychotherapy and Psychosomatics.}\url{https://doi.org/10.1159/000494356}
#'
#' @importFrom parallel detectCores
#'
#' @export
#'
#' @examples
#' phase1=cbind(rnorm(50,0,1),rnorm(50,0,1)) #phase1: Means=0
#' phase2=cbind(rnorm(50,1,1),rnorm(50,1,1)) #phase2: Means=1
#' X=rbind(phase1,phase2)
#'
#' #scan all statistics
#' \donttest{
#' res=kcpRS_workflow(data=X,RS_funs=c("runMean","runVar","runAR","runCorr"),
#' wsize=25,nperm=1000,Kmax=10,alpha=.05, varTest=FALSE, bcorr=TRUE, ncpu=1)
#' summary(res)
#' plot(res)
#'
#'
#' #scan the mean and the correlation only
#' res=kcpRS_workflow(data=X,RS_funs=c("runMean","runCorr"),wsize=25,nperm=1000,Kmax=10,
#'     alpha=.05, varTest=FALSE, bcorr=TRUE, ncpu=1)
#' summary(res)
#' plot(res)
#' }

kcpRS_workflow<-function(data,RS_funs=c("runMean","runVar","runAR","runCorr"),wsize=25,nperm=1000,Kmax=10,alpha=.05, varTest=FALSE, bcorr=TRUE,ncpu=1){
    UseMethod("kcpRS_workflow")
}
