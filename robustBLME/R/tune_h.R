#' @importFrom lme4 lmer getME VarCorr
#' @name tune.h
#' @title Tune ABC distance bandwidth
#' @description  Tunes the bandwidth \eqn{h} of the ABC distance to get the desired level of acceptance rate specified via \code{acc.rate}. Besides tuning \eqn{h}, the function also builds the relevant quantities needed for running \code{\link[robustBLME]{rblme}}. For generating such quantities an internal call to \code{\link[lme4]{lmer}} is performed.
#'
##' @usage tune.h(formula, data, ..., n.samp = 1e+5, n.sim.HJ = 500, acc.rate, grid.h, prior,
##'        cHub = 1.345, cHub2 = 2.07,
##'        init, n.cores = 1, use.h)
##' @param formula two-sided linear formula object describing the
##'   fixed-effects part of the model, with the response on the left of
##'   a \code{~} operator and the terms, separated by \code{+}
##'   operators on the right.  The \code{"|"} character separates an
##'   expression for a model matrix and a grouping factor.
##' @param data optional data frame containing the variables named
##'   in \code{formula}. By default the variables are taken from the
##'   environment of \code{\link[lme4]{lmer}} called internally.
##' @param ... other arguments to be passed to \code{lmer}. Currently none is used.
##' @param n.samp number of pilot posterior samples to be drawn with ABC for each value of \code{grid.h}.
##' @param n.sim.HJ number of simulations to be used for computing the sensitivity and variability matrices.
##' @param acc.rate desired acceptance rate of the ABC-MCMC algorithm.
##' @param grid.h grid of \eqn{h} values within which the "optimal" value is to be found.
##' @param prior named list of user-defined prior hyper-parameters. See "Details" below.
##' @param cHub tuning constant of the Huber function for the location parameter.
##' @param cHub2 tuning constant of the Huber proposal 2 function for the scale parameter.
##' @param init optional object to use for starting values. Currently ignored as initial values are taken from \code{lmer}.
##' @param n.cores number of cores for parallel computation on non Windows machines.
##' @param use.h bandwidth to be used for the ABC distance. If provided, no tuning for \eqn{h} is performed and \code{acc.rate} is ignored.
#' @return list.
#' @details Given a specification of the \code{formula} and \code{data} the function calls internally \code{rlmer} and extracts from the resulting object all the necessary quantities. Then proceeds by finding the solution of the REML II robust estimating equations (Richardson & Welsh 1995), with the REML estimate used as starting point. The sensitivity and the variability matrices are computed by simulation at the solution of the robust REML II estimating equation. Depending on whether \code{use.h} or \code{acc.rate} and \code{grid.h} are specified, the function has a different behaviour. If \code{acc.rate} and \code{grid.h} are provided, then an adaptive step is performed in order to get an "optimal" \eqn{h} which gives the desired acceptance rate \code{acc.rate}. IN particular, for each value of \code{grid.h}, the function draws \code{n.samp} posterior samples with the ABC-MCMC algorithm and saves the resulting acceptance rate. Lastly, a function is built via a smoothing spline with acceptance rates being the \eqn{x}s and \code{grid.h} being the \eqn{y}s. The "optimal" value of \eqn{h} is found, within \code{grid.h}, as the prediction the spline function at \code{acc.rate}. If you already have an \eqn{h} value in mind then specify it via \code{use.h} and leave \code{grid.h} and \code{acc.rate} unspecified. Note that, in this case the acceptance rate of the ABC-MCMC algorithm may not be the one you wish to obtain since it depends in some complicated way also from \eqn{use.h}. Currently, the prior for the \eqn{q} fixed effects is the product of \eqn{q} scalar normals with mean zero and user-specified variance \code{beta.sd} (see Examples) equal for all the parameters. For the variance components the prior is a halfCauchy with user-specified scale \code{s2.scale}. Both variance parameters are assumed to have equal prior scale.
#'
##' @references
##' Ruli E., Sartori N. & Ventura L. (2017)
##' Robust approximate Bayesian inference with an application to linear mixed models.
##' \url{https://arxiv.org/abs/1706.01752}
##'
##' Richardson A. M. & Welsh A. H. (1995) Robust restricted maximum likelihood in mixed linear models. \emph{Biometrics} \bold{51}, 1429-1439.
##'
##' @seealso \code{\link[robustBLME]{rblme}}, \code{\link[robustBLME]{ergoStool}}.

#' @export
tune.h <- function(formula, data, ..., n.samp = 1e+5, n.sim.HJ = 500, acc.rate, grid.h,
                   prior, cHub = 1.345, cHub2 = 2.07,
                   init, n.cores = 1, use.h){

  if((missing(acc.rate) && missing(use.h)) || (!missing(acc.rate) && !missing(use.h)))
    stop("Either \'acc.rate\' or \'use.h\' must be provided but not both.")

  if(!missing(acc.rate)){
    if((acc.rate>=1 || acc.rate <= 0)) {
      stop("\'acc.rate must be in (0,1)")
    }
    if(missing(grid.h) || !is.vector(grid.h, mode = "numeric")){
      stop("\'grid.h\' must be a numeric vector.")
    }
  }

  if(!missing(use.h) && (use.h<=0))
    stop("\'use.h\' must be strictly positive")

  if(missing(prior) || !is.list(prior) || length(prior) != 2)
    stop("prior hyper-parameters must be a named list of two elements")

  rcall = match.call()

  rcall[setdiff(names(formals(tune.h)), names(formals(lmer)))] = NULL
  rcall$doFit = NULL
  rcall$REML = TRUE
  rcall[[1]] = as.name("lmer")
  pf = parent.frame()

  # finding the REML estimates
  rinit = eval(rcall, pf)

  # getting the necessary quantites from lmer()
  obj = getME(rinit, name = c("Zt", "X", "y", "beta", "flist", "k"))

  if(obj$k != 1) stop("I'm sorry but models with >1 random effects are currently not supported in this version. We hope to relax this assumption in the future versions.")

  start.betas = obj$beta
  start.sig2 = as.data.frame(VarCorr(rinit))$vcov
  n.betas = length(start.betas)
  Xn = obj$X
  # p = ncol(Xn)
  n = nrow(Xn)
  n_ind = length(levels(obj$flist[[1]]))
  m_i = n/n_ind

  Xnt = t(Xn)

  Zt_b = as.matrix(obj$Zt)
  ZZt_b = crossprod(Zt_b)
  ZZt_b_ii = ZZt_b[1:m_i, 1:m_i]
  Z_eps = diag(1, nrow = n, ncol = n)
  ZZt_eps = crossprod(Z_eps)
  ZZt_eps_ii = ZZt_eps[1:m_i, 1:m_i]

  # n_sig2 = obj$k + 1

  y = t(matrix(obj$y, n_ind, m_i, byrow=T))

  ###### K2n matrix consistency correction ########
  K2n <- K2mat(cHub2, n)

  ####### solve the estimating equtions ###########
  sol <- rootREML2_GH(start = c(start.betas, start.sig2),
                      n_iter = 80,
                      sig2b_interval = c(0.0001, 50),
                      sig2eps_interval = c(0.0001, 50),
                      y = y,
                      Xn = Xn,
                      Xnt = Xnt,
                      ZZt_b = ZZt_b,
                      ZZt_b_ii = ZZt_b_ii,
                      ZZt_eps = ZZt_eps,
                      ZZt_eps_ii = ZZt_eps_ii,
                      c_hub = cHub,
                      c2_hub = cHub2,
                      K2n = K2n)
  hat.beta <- sol$hat_beta
  hat.sig2 <- sol$hat_sig2
  Psi.hat <- sol$Psi_hat

  ####### The H and J matrices ###########
  sim_HJ <- sim_SensVar_mat_lmm(betas = hat.beta,
                                sig2_b = hat.sig2[1],
                                sig2_eps = hat.sig2[2],
                                Xn = Xn,
                                Xnt = Xnt,
                                ZZt_b = ZZt_b,
                                ZZt_b_ii = ZZt_b_ii,
                                ZZt_eps = ZZt_eps,
                                ZZt_eps_ii = ZZt_eps_ii,
                                c_hub = cHub,
                                c2_hub = cHub2,
                                K2n = K2n,
                                n_sim = n.sim.HJ,
                                m_i = m_i,
                                n = n,
                                n_ind = n_ind)


  # Godembe information matrix
  Godambe_info <- solve(solve(sim_HJ$H)%*%sim_HJ$J%*%solve(t(sim_HJ$H)))

  # Jacobian matrix for the reparametrization sigma2 -> log(sigma2)
  Jacob <- matrix(0, n.betas+2, n.betas+2)
  diag(Jacob) <- c(rep(1,n.betas),hat.sig2)

  # asymptotic covariance matrix of the robust estimate
  # to be used as scale matrix for the proposal but also for Wald inference!
  Gvar.rep <- solve(Jacob%*%Godambe_info%*%Jacob)

  # covariance matrix of the REML II estimating equation
  # that we as scaling matrix for the ABC distance
  J_Psi.hat <- solve(sim_HJ$J)

  lch_prop_scale_mat <- t(chol(Gvar.rep))
  lchol_J_Psi_hat = chol(J_Psi.hat)

  n.h = length(grid.h)
  eff.h = rep(NA, n.h)

  thin.h <- which(rep(1:n.h, len = n.samp) == n.h) - 1

  if(!missing(acc.rate)){

    if(n.cores < 2){
      cat('\n#####################################################\n Optimising h to get',
          acc.rate*100, '% acceptance\n')

      progress_bar <- txtProgressBar(min = 0,
                                     max = n.h,
                                     style = 3)

      i = 1
      for(i in 1:n.h){
        eff.h[i] = .Call('robustBLME_ABCkern_reml2',
              PACKAGE = 'robustBLME',
              n.samp,
              grid.h[i],
              y,
              Xn,
              ZZt_b,
              ZZt_b_ii,
              ZZt_eps,
              ZZt_eps_ii,
              cHub,
              cHub2,
              K2n,
              c(hat.beta, log(hat.sig2)),
              lch_prop_scale_mat,
              Psi.hat,
              lchol_J_Psi_hat,
              prior$beta.sd,
              prior$s2.scale,
              thin.h,
              trace_int = n.samp,
              m_i,
              n_ind)$effi
              setTxtProgressBar(progress_bar, i)
      }


    } else {

  cat('\n#####################################################\n Optimising h to get',
      acc.rate*100, '% acceptance\n')

  registerDoParallel(cores = n.cores)
  eff.h = foreach(i = 1:n.h, .combine = c) %dopar% {
    .Call('robustBLME_ABCkern_reml2',
          PACKAGE = 'robustBLME',
          n.samp,
          grid.h[i],
          y,
          Xn,
          ZZt_b,
          ZZt_b_ii,
          ZZt_eps,
          ZZt_eps_ii,
          cHub,
          cHub2,
          K2n,
          c(hat.beta, log(hat.sig2)),
          lch_prop_scale_mat,
          Psi.hat,
          lchol_J_Psi_hat,
          prior$beta.sd,
          prior$s2.scale,
          thin.h,
          trace_int = n.samp,
          m_i,
          n_ind)$effi
    }
  stopImplicitCluster()
  }

    sp.effh <- splinefun(x = grid.h, y = eff.h)

    optimal.h <- uniroot(function(x) sp.effh(x) - acc.rate,
                       interval = range(grid.h))$root

  plot(sp.effh,
       xlim = range(grid.h), xlab = "h",
       ylab = "acceptance rate", lwd = 2,
       main = paste("Optimal h = ",
                    round(optimal.h,2), "; actual acc.rate = ",
                    round(sp.effh(optimal.h), 4)*100, "%"))
  abline(v = optimal.h, col = "gray", lwd = 2, lty = "dashed")
  abline(h = sp.effh(optimal.h), col = "gray", lwd = 2, lty = "dashed")

  cat('\n Optimal h = ', round(optimal.h,2), 'actual acceptance rate = ',
      round(sp.effh(optimal.h), 4)*100, "%")

  out <- list(h.opt = optimal.h,
              y = y,
              Xn = Xn,
              ZZt_b = ZZt_b,
              ZZt_b_ii = ZZt_b_ii,
              ZZt_eps = ZZt_eps,
              ZZt_eps_ii = ZZt_eps_ii,
              cHub = cHub,
              cHub2 = cHub2,
              K2n = K2n,
              param.hat = c(hat.beta, log(hat.sig2)),
              lch_prop_scale_mat = lch_prop_scale_mat,
              Psi.hat = Psi.hat,
              lchol_J_Psi_hat = lchol_J_Psi_hat,
              prior = prior,
              m_i = m_i,
              n_ind = n_ind)
  } else {
    cat('\n#####################################################\n Running preliminary ABC-MCMC with h = ', use.h, '\n')
    out <- list(h.opt = optimal.h,
                y = y,
                Xn = Xn,
                ZZt_b = ZZt_b,
                ZZt_b_ii = ZZt_b_ii,
                ZZt_eps = ZZt_eps,
                ZZt_eps_ii = ZZt_eps_ii,
                cHub = cHub,
                cHub2 = cHub2,
                K2n = K2n,
                param.hat = c(hat.beta, log(hat.sig2)),
                lch_prop_scale_mat = lch_prop_scale_mat,
                Psi.hat = Psi.hat,
                lchol_J_Psi_hat = lchol_J_Psi_hat,
                prior = prior,
                m_i = m_i,
                n_ind = n_ind)
  }
  class(out) <- "rblmeMod"

  return(out)
}