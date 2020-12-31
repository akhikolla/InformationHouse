##' @useDynLib robustBLME, .registration = TRUE
##' @import Rcpp
##' @import doParallel
##' @import foreach
##' @import iterators
##' @importFrom numDeriv jacobian
##' @importFrom utils txtProgressBar setTxtProgressBar
##' @importFrom graphics abline plot
##' @importFrom stats splinefun uniroot integrate dnorm median nlminb density
##' @importFrom mvtnorm rmvnorm
##' @importFrom parallel mclapply
##' @importFrom lme4 lmer VarCorr getME
##' @importFrom stats getCall
##' @title Fits robust Bayesian Linear Mixed-effects Models (LMM) to data via robust REML estimating functions.
##'
##' @description This function fits robust Bayesian LMMs to data via robust REML estimating functions. The latters are those proposed by Richardson & Welsh (1995), which are robustified versions of restricted maximum likelihood (REML) estimating equations. Posterior sampling is done with an ABC-MCMC algorithm, where the data are summarised through a rescaled version of the aforementioned estimating functions; see Ruli et al. (2017) for the properties and details of the method. The current package version (0.1.2) supports only models with a single random effects. Extensions to more general settings will be provided in the future versions of the package.
##'
##' @usage rblme(nabc, h.obj, chain.control = list(trace.init = NULL, thin.by = NULL),
##'         n.cores = 1)
##' @param nabc number of posterior samples.
##' @param h.obj list of objects as returned by the \code{\link[robustBLME]{tune.h}} function. Hence \code{tune.h} must be called first.
##' @param chain.control parameters that control the tracing and the thinning of the chain(s).
##' @param n.cores number of cores for parallel computation on non Windows machines. For \code{n.cores}>2, \code{n.cores} chains are run each on a different core with using the same parameters but with a different random seed.
##'
##' @return list or list of lists with elements \code{abc} and \code{effi}. In case of \code{n.cores}=1, \code{effi} is the actual acceptance rate of the ABC-MCMC algorithm whereas in \code{abc} are stored the posterior samples. The latters are stored as a \eqn{(q + c) \times}nabc matrix, where \eqn{q} is the number of fixed effects, i.e. the number of columns in the design matrix and \eqn{c = 2} is the number of variance components. Hence, the first \eqn{q} rows of the matrix \code{abc} give the posterior samples for the fixed effects and the last two rows give the posterior samples for the log-variances of the fixed effects and the residual term, respectively. If \code{n.cores} > 1, i.e. if simulations are performed in parallel, then a list of lists is returned, where each element of the list is a list with elements \code{abc} and \code{effi}, where \code{abc} and \code{effi} are as those aforementioned.
##'
##'
##' @references
##' Ruli E., Sartori N. & Ventura L. (2017)
##' Robust approximate Bayesian inference with an application to linear mixed models.
##' \url{https://arxiv.org/abs/1706.01752}
##'
##' Richardson A. M. & Welsh A. H. (1995) Robust restricted maximum likelihood in mixed linear models. \emph{Biometrics} \bold{51}, 1429-1439.
##'
##' @seealso \code{\link[robustBLME]{tune.h}}, \code{\link[robustBLME]{ergoStool}}.
##' @examples
##'
##' ## The following example is meant for function documentation.
##' ## For realistic use probably you'll need to take a larger sample and choose a
##' ## "better" bandwidth h.
##'
##' data(ergoStool)
##'
##' require(lme4)
##' fm1 <- lmer(effort~Type + (1| Subject), data = ergoStool)
##'
##' ## tune h to get 0.8% acceptance
##' hopt <- tune.h(effort~Type + (1|Subject), data = ergoStool, n.samp = 1e+4,
##'                acc.rate = 0.01, n.sim.HJ = 100, grid.h = seq(0.3, 0.7, len = 3),
##'                prior = list(beta.sd = 10, s2.scale = 5), n.cores = 1)
##'
##' ## draw posterior samples with hopt.
##' abc.tmp <- rblme(nabc = 1e+4, h.obj = hopt,
##'                  n.cores = 1)
##'
##' # process ABC samples
##' abc.sim <- t(abc.tmp$abc)
##' abc.sim[,c(5,6)] <- exp(abc.sim[,c(5,6)])
##'
##' # ABC posterior
##' colMeans(abc.sim)
##'
##' # REML estimates
##' summary(fm1)
##'
##' @export
rblme <- function(nabc,
                  h.obj,
                  chain.control = list(trace.init = NULL, thin.by = NULL),
                  n.cores = 1){

  if(missing(h.obj))
      stop("\'h.obj\' must be a valid object as returned by \'tune.h()\'")

  if(mode(chain.control$trace.int) != "numeric")
    chain.control$trace.int = floor(nabc/10)


  if(mode(chain.control$thin.by) != "numeric")
    chain.control$thin.by = floor(1/h.obj$h.opt)

  thin = which(rep(1:chain.control$thin.by, len = nabc) == chain.control$thin.by) - 1

  # parallel computing via forking under windows doesn't work!
  if(.Platform$OS.type != "unix" || n.cores < 2){
    out <- .Call('robustBLME_ABCkern_reml2',
                 PACKAGE = 'robustBLME',
                 nabc,
                 h.obj$h.opt,
                 h.obj$y,
                 h.obj$Xn,
                 h.obj$ZZt_b,
                 h.obj$ZZt_b_ii,
                 h.obj$ZZt_eps,
                 h.obj$ZZt_eps_ii,
                 h.obj$cHub,
                 h.obj$cHub2,
                 h.obj$K2n,
                 h.obj$param.hat,
                 h.obj$lch_prop_scale_mat,
                 h.obj$Psi.hat,
                 h.obj$lchol_J_Psi_hat,
                 h.obj$prior$beta.sd,
                 h.obj$prior$s2.scale,
                 thin,
                 trace_int = chain.control$trace.int,
                 h.obj$m_i,
                 h.obj$n_ind)

  } else {
    out <- mclapply(X = 1:n.cores,
                    FUN = function(i)
                      .Call('robustBLME_ABCkern_reml2',
                                PACKAGE = 'robustBLME',
                            nabc,
                            h.obj$h.opt,
                            h.obj$y,
                            h.obj$Xn,
                            h.obj$ZZt_b,
                            h.obj$ZZt_b_ii,
                            h.obj$ZZt_eps,
                            h.obj$ZZt_eps_ii,
                            h.obj$cHub,
                            h.obj$cHub2,
                            h.obj$K2n,
                            h.obj$param.hat,
                            h.obj$lch_prop_scale_mat,
                            h.obj$Psi.hat,
                            h.obj$lchol_J_Psi_hat,
                            h.obj$prior$beta.sd,
                            h.obj$prior$s2.scale,
                            thin,
                            trace_int = chain.control$trace.int,
                            h.obj$m_i,
                            h.obj$n_ind))
  }

  return(out)
}