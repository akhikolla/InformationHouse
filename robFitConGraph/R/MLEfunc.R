#' Graph-constrained robust scatter estimation.
#'
#' The function computes a robust estimate of a scatter matrix subject to
#' zero-constraints in its inverse. The methodology is described in Vogel &
#' Tyler (2014).
#'
#' @param X A data matrix with \eqn{n} rows and \eqn{p} columns, representing
#'   \eqn{n} observations and \eqn{p} variables. Elements of \code{X} must be
#'   numeric and \eqn{n} must be at least \eqn{p+1}.
#'
#' @param amat A \eqn{p} times \eqn{p} matrix representing the adjacency matrix
#'   of a graphical model. \code{amat} must be symetric with numerical entries
#'   \code{0} or \code{1}. The entries on the diagonal are irrelevant, they may
#'   be anything.
#'
#' @param df the degrees of freedom of the t-distribution used (see Details
#'   below).
#'
#' @param tol tolerance for numerical convergence. Iteration stops if the
#'   maximal element-wise difference between two successive matrices is less
#'   than \code{tol}. Must be at least \code{10e-14}. Default is \code{10e-5}.
#'
#' @param plug.in logical. The function offers two types of estimates: the
#'   plug-in M-estimator and the direct M-estimator. If \code{plug.in} is
#'   \code{TRUE}, the plug-in estimate is computed. If \code{FALSE}, the direct
#'   M-estimator is computed. The plug-in estimator is faster, but has higher variance.
#'   Default is \code{TRUE}.
#'
#' @param direct logical. If \code{TRUE}, the direct estimate is computed,
#'   otherwise the plug-in estimate. Default is \code{FALSE}. In case of
#'   conflicting specifications of plug-in and direct, \code{plug.in} overrides
#'   \code{direct}.
#'
#' @return List with 5 elements:
#'   \item{\code{Shat}}{\code{p} x \code{p} scatter matrix estimate}
#'   \item{\code{mu}}{numerical \code{p}-vector (robust location estimate)}
#'   \item{\code{em.it}}{integer. Number of iterations of the t-MLE
#'   computation.}
#'   \item{\code{ips.it}}{integer. In the case of the plug-in estimate, this is
#'   the number of iterations of the Gaussian graphical model fitting procedure
#'   (Algorithm 17.1) in Hastie et al 2004). In the case of the direct estimate,
#'   the Gaussian graphical model fitting is executed \code{em.it} times and the
#'   average number of iterations is returned.}
#'   \item{\code{dev}}{numerical. Value of the deviance test statistic
#'   \eqn{D_n} as defined in Vogel & Tyler (2014, p. 866 bottom).
#'   Comparing the model fitted against the full model.}
#'
#'
#' @details The function \code{robFitConGraph} implements the methodology of
#'   Vogel & Tyler (2014). Two types of estimates based on maximum likelihood
#'   estimation for the t-distribution are proposed: the direct estimate and the
#'   plug-in estimate. The direct estimate is referred to as graphical
#'   M-estimator in Vogel & Tyler (2014).
#'
#'   The plug-in estimate is two algorithms performed sequentially: First an
#'   unconstrained t-maximum likelihood estimate of scatter is computed (the
#'   same as \code{\link[MASS]{cov.trob}} from \code{MASS}). This is then
#'   plugged into the Gaussian graphical model fitting routine (the same as
#'   \code{\link[ggm]{fitConGraph}} from \code{ggm}). Specifically the
#'   algorithm 17.1 from Hastie, Tibshirani, Friedman (2009) is used.
#'
#'   The direct estimate is the actual maximum-likelihood estimator within the
#'   elliptical graphical model based on the elliptical t-distribution. The
#'   algorithm is an iteratively-reweighted least-squares algorithm, where the
#'   Gaussian graphical model fitting procedure is nested into the t-estimation
#'   iteration. The direct estimate therefore takes longer to compute, but the
#'   estimator has a better statistical efficiency for small sample sizes. Both
#'   estimators are asymptotically equivalent. The estimates tend to be very
#'   close to each other for large sample sizes.
#'
#'   Although \code{robFitConGraph} combines the functionality of
#'   \code{\link[ggm]{fitConGraph}} and \code{\link[MASS]{cov.trob}} and
#'   contains both as special cases, it uses only the latter function. The
#'   algorithms are largely implemented in C++.
#'
#'   Input and output of \code{robFitConGraph} are similar to
#'   \code{\link[ggm]{fitConGraph}} from the package \code{ggm}. Some notable
#'   differences:
#'
#'   \itemize{ \item \code{\link[ggm]{fitConGraph}} takes as input the
#'   unconstrained covariance matrix, \code{robFitConGraph} takes the actual data.
#'   \item \code{\link[ggm]{fitConGraph}} returns the deviance (test statistic)
#'   and the degrees of freedom \eqn{r}. The degrees of freedom \eqn{r} are the
#'   number of sub-diagonal 1-entries in the adjacency matrix. The deviance is
#'   compared to a chi-square distribution with \eqn{r} degrees of freedom to
#'   assess the model fit. These degrees of freedom \eqn{r} are unrelated to the
#'   the parameter \code{df}, which refers to the degrees of freedom of the
#'   t-distribution. The function \code{robFitConGraph} does return the
#'   deviance, but no  degrees of freedom. The deviance must be divided by a
#'   constant (\eqn{\sigma_1} in Vogel & Tyler, 2014) before comparing it to the
#'   \eqn{\chi^2_r}-distribution.}
#'
#' @author Stuart Watt, Daniel Vogel
#'
#' @references Vogel, D., Tyler, D. E. (2014): Robust estimators
#'   for nondecomposable elliptical graphical models, \emph{Biometrika}, 101, 865-882\cr\cr
#'   Hastie, T., Tibshirani, R. and Friedman, J. (2004). \emph{The elements of
#'   statistical learning}. New York: Springer.
#'
#' @seealso \code{\link[ggm]{fitConGraph}} from package
#'   \code{ggm} for non-robust graph-constrained covariance estimation
#'   \cr\cr
#'   \code{\link[MASS]{cov.trob}} from package \code{MASS} for unconstrained
#'   \code{p} times \code{p} t-MLE scatter matrix
#'
#' @examples
#' # --- build a graphical model ---
#'
#' chordless.p.cycle <- function(rho,p){
#'   M <- diag(1,p)
#'   for (i in 1:(p-1)) M[i,i+1] <- M[i+1,i] <- -rho
#'   M[1,p] <- M[p,1] <- -rho
#'   return(M)
#' }
#' p <- 7                             # number of variables
#' rho <- 0.4                         # partial correlation
#' PCM <- chordless.p.cycle(rho,p)    # partial correlation matrix
#' SM <- cov2cor(solve(PCM))          # shape matrix (i.e covariance matrix up to scale)
#' model <- abs(sign(PCM))            # adjacency matrix of the chordless-7-cycle
#' # > model
#' #      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#' # [1,]    1    1    0    0    0    0    1
#' # [2,]    1    1    1    0    0    0    0
#' # [3,]    0    1    1    1    0    0    0
#' # [4,]    0    0    1    1    1    0    0
#' # [5,]    0    0    0    1    1    1    0
#' # [6,]    0    0    0    0    1    1    1
#' # [7,]    1    0    0    0    0    1    1
#'
#' # This is the cordless-7-cycle (p.872 Figure 1 (a) in Vogel & Tyler, 2014).
#' # All non-zero partial correlations are 0.4.
#' # The true covariance is (up to scale) 'SM'. This matrix is constructed such
#' # that it has zero entries in its inverse as specified by 'model'.
#'
#'
#' # --- generate data from the graphical model ---
#'
#' n <- 50            # number of observations
#' df.data <- 3       # degrees of freedom
#' library(mvtnorm)   # for rmvt function
#' set.seed(918273)   # for reproducability
#' X <- rmvt(n=n,sigma=SM,df=df.data)
#'
#' # X contains a data set of size n = 50 and dimension p = 7, sampled from the
#' # elliptical t-distribution with 3 degrees of freedom and shape matrix 'SM'
#'
#'
#' # --- compare estimates ---
#'
#' # We compute three scatter estimates:
#'
#' # 1) the direct graph-constrained t-MLE estimator:
#' S1 <- robFitConGraph(X, amat=model, df=df.data, plug.in=FALSE, direct=TRUE)$Shat
#' round(S1,d=2)
#'
#' # 2) the plug-in graph-constrained t-MLE estimator:
#' S2 <- robFitConGraph(X, amat=model, df=df.data, plug.in=TRUE, direct=FALSE)$Shat
#' round(S2,d=2)
#'
#' # 3) the saple covariance matrix:
#' round(cov(X),d=2)
#'
#' # S1 and S2 are very similar. In Vogel & Tyler, 2014, it is shown that they
#' # are asymptotically equivalent as n goes to infinity.
#' # The sample covariance matrix substantially differs from S1 and S2. Note that
#' # S1 and S2 just estimate a multiple of the true covariance matrix (similarly
#' # SM is just proportional to the true covariance matrix). Therefore, consider
#' # correlation estimates based on the various scatter estimators:
#'
#' # the true correlation matrix:
#' round(cov2cor(SM),d=2)
#'
#' # sample correlations:
#' round(cov2cor(cov(X)),d=2)
#'
#' # robust correlations based on the direct graph-constrained t-MLE:
#' round(cov2cor(S1),d=2)
#'
#' # robust correlations based on the plug.in graph-constrained t-MLE:
#' round(cov2cor(S2),d=2)
#'
#' # The correlation estimates based on S1 and S2 are close to the true
#' # correlations, whereas the sample correlations, again, differ strongly.
#' # Note: sample correlations are not asymptotically normal at the t3 distribution.
#'
#' @export
robFitConGraph <- function(X,amat,df,tol=10e-5,plug.in=TRUE,direct=FALSE){
  #--- CHECK ---

  #- X check -
  if(missing(X)) stop('data matrix X required')
  if(is.double(X) == FALSE) stop('data matrix X must be double data type')

  #- amat check -
  if(missing(amat)) stop('amat matrix needs to be specified')
  if(isSymmetric(amat) == FALSE) stop('amat is not symmetric')
  if(length(amat[col(amat)!=row(amat)][!(amat[col(amat)!=row(amat)] %in% c(0,1))]) != 0) stop('off-diagonal amat entries must be 0 or 1')

  #- dimension check -
  if((ncol(X) != ncol(amat))) stop('number of columns of X and amat must match')
  if(ncol(amat) < 2) stop('number of variables must be greater than or equal to 2')
  if(nrow(X) <= ncol(X)) stop('number of observations must be greater than number of variables')

  #- df check -
  if(missing(df)) stop('df needs to be specified')
  if(is.double(df) == FALSE) stop('data input df must be double data type')
  if(df <= 0) stop('data input df must be positive')

  #- tol check -
  if(is.double(tol) == FALSE) stop('data input tol must be double data type')
  if((tol < 10e-14)) stop('data input df must be >= 10e-14')

  #- algorithm specification check
  if((missing(plug.in) && missing(direct))) warning('Algorithm not specified. plug.in = TRUE chosen by default.')
  if((missing(plug.in) || missing(direct)) == FALSE && (plug.in == TRUE) && (direct == TRUE)) warning('Conflicting specifications for plug.in and direct. plug.in = TRUE chosen by default.')
  if((plug.in == FALSE) && (direct == FALSE) && (missing(plug.in) || missing(direct)) == FALSE) warning('Conflicting specifications for plug.in and direct. plug.in = FALSE chosen by default.')

  #--- Check Complete ---

  if(df == Inf){
    Cov <- myIPSC(model=amat,estimate=stats::cov(X),tol=tol,n=nrow(X))
    return(list(Shat=Cov$Shat,mu=Cov$center,em.it=0,ips.it=Cov$it,dev=Cov$dev))
  }else if(((missing(plug.in) == TRUE) && (missing(direct) == TRUE)) ||
      ((plug.in == TRUE) && (missing(direct) == TRUE)) ||
      ((direct == FALSE) && (missing(plug.in) == TRUE)) ||
      ((plug.in == TRUE) && (direct == FALSE) && (missing(plug.in) || missing(direct)) == FALSE) ||
      ((plug.in == TRUE) && (direct == TRUE) && (missing(plug.in) || missing(direct)) == FALSE)){
    # GCMLE
    Cov <- MASS::cov.trob(X,nu=df,tol=tol/10,maxit=1000)
    CovG <- myIPSC(model=amat,estimate=Cov$cov,tol=tol,n=nrow(X))
    return(list(Shat=CovG$Shat,mu=Cov$center,em.it=Cov$iter,ips.it=CovG$it,dev=CovG$dev))
  }else{
    # ECMLE
    d <- function(A,B){return(max(abs(A-B)))}
    p <- ncol(X)
    n <- nrow(X)
    S0 <- MASS::cov.trob(X,nu=df,tol=tol/10,maxit=1000)$cov
    S.alt <- S0 # diag(1,p)
    S.neu <- myIPSC(model=amat,estimate=stats::cov(X),tol=tol)$S
    mu.neu <- colMeans(X)
    count <- 0
    counts.ips <- list()
    while (d(S.neu, S.alt) > tol) {
      S.alt <- S.neu
      mu.alt <- mu.neu
      ss <- stats::mahalanobis(X,center=mu.alt,cov=S.alt)
      us <-  sapply(X=ss,FUN=u,p=p,d=df)
      mu.neu <- as.vector((us %*% X)/sum(us))
      S.neu.vector <- CapplyB(us=us, X=X, mu=mu.neu)/n
      Cov <- myIPSC(model=amat,estimate=matrix(S.neu.vector,ncol=p),tol=tol)
      S.neu <- Cov$S
      counts.ips[count] <- Cov$it
      count <- count+1
    }
    dev <- n*(log(det(S.neu)) - log(det(S0)))
    return(list(Shat=S.neu,mu=mu.neu,em.it=count,ips.it=mean(as.numeric(counts.ips)),dev=dev))
  }
}
