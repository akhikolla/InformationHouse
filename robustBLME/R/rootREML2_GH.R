# function used internally by "rootREML2_GH"

.beta_given_rest<- function(betas,
                           rest,
                           y,
                           Xn,
                           ZZt_b,
                           ZZt_b_ii,
                           ZZt_eps,
                           ZZt_eps_ii,
                           c_hub,
                           n_betas,
                           m_i,
                           n,
                           n_ind) {

  oo = Psi_reml2_betas(betas = betas,
                       sig2_b = rest[1],
                       sig2_eps = rest[2],
                       y = y,
                       Xn = Xn,
                       ZZt_b = ZZt_b,
                       ZZt_b_ii = ZZt_b_ii,
                       ZZt_eps = ZZt_eps,
                       ZZt_eps_ii = ZZt_eps_ii,
                       c_hub = c_hub,
                       n_betas = n_betas,
                       m_i = m_i,
                       n = n,
                       n_ind = n_ind)
  crossprod(oo)
}

# function used internally by "rootREML2_GH"
.sig2b_given_rest <- function(sig2_b,
                             rest,
                             y,
                             Xn,
                             Xnt,
                             ZZt_b,
                             ZZt_b_ii,
                             ZZt_eps,
                             ZZt_eps_ii,
                             c2_hub,
                             K2n,
                             n_betas,
                             m_i,
                             n,
                             n_ind) {

  Psi_reml2_sig2_b(betas = rest[1:n_betas],
                   sig2_b = sig2_b,
                   sig2_eps = rest[n_betas+1],
                   y = y,
                   Xn = Xn,
                   Xnt = Xnt,
                   ZZt_b = ZZt_b,
                   ZZt_b_ii = ZZt_b_ii,
                   ZZt_eps = ZZt_eps,
                   ZZt_eps_ii = ZZt_eps_ii,
                   c2_hub = c2_hub,
                   K2n = K2n,
                   n_betas = n_betas,
                   m_i = m_i,
                   n = n,
                   n_ind = n_ind)
}


# function used internally by "rootREML2_GH"
.sig2eps_given_rest <- function(sig2_eps,
                               rest,
                               y,
                               Xn,
                               Xnt,
                               ZZt_b,
                               ZZt_b_ii,
                               ZZt_eps,
                               ZZt_eps_ii,
                               c2_hub,
                               K2n,
                               n_betas,
                               m_i,
                               n,
                               n_ind) {

  Psi_reml2_sig2_eps(betas = rest[1:n_betas],
                     sig2_b = rest[n_betas+1],
                     sig2_eps = sig2_eps,
                     y = y,
                     Xn = Xn,
                     Xnt = Xnt,
                     ZZt_b = ZZt_b,
                     ZZt_b_ii = ZZt_b_ii,
                     ZZt_eps = ZZt_eps,
                     ZZt_eps_ii = ZZt_eps_ii,
                     c2_hub = c2_hub,
                     K2n = K2n,
                     n_betas = n_betas,
                     m_i = m_i,
                     n = n,
                     n_ind = n_ind)
}


##' @title Solves robust REML II by Gauss-Seidel iterations.
##' @description Implements the Gauss-Seidel algorithm for finding the root of robust REML II estimating equation.
##' @usage rootREML2_GH(start, n_iter, sig2b_interval, sig2eps_interval, y,
##'              Xn, Xnt, ZZt_b, ZZt_b_ii, ZZt_eps, ZZt_eps_ii, c_hub,
##'              c2_hub, K2n, trace = TRUE)
##'
##' @param start Starting value
##' @param n_iter Number of allowed iterations.
##' @param sig2b_interval Interval in which to look for the solution of variance of the random effects.
##' @param sig2eps_interval Interval in which to look for the solution of variance of the error term.
##' @param y The response variable.
##' @param Xn The design matrix for all statistical units.
##' @param Xnt The transpose of Xn.
##' @param ZZt_b The ZZt matrix for the random effects.
##' @param ZZt_b_ii The ii block of ZZt_b corresponding to the i-th unit.
##' @param ZZt_eps The ZZt matrix for the error term.
##' @param ZZt_eps_ii The ii block of ZZt_eps corresponding to the i-th unit.
##' @param c_hub The tuning constant of the Huber function for the location.
##' @param c2_hub The tuning constant of the Huber function for the scale.
##' @param K2n The K2n matrix needed for consistency correction.
##' @param trace Should the output be prit on the screen?
##' @details The algorithm works by solving the blocks by part. In particular, given a starting value the algorithm starts by solving for the fixed effects, then iterates among the variance components with the fixed effects held fixed at the updated value. The interation proceeds until the solution stabilises withing a pre-specified threshold level.
##' @noRd
rootREML2_GH <- function(start,
                         n_iter,
                         sig2b_interval,
                         sig2eps_interval,
                         y,
                         Xn,
                         Xnt,
                         ZZt_b,
                         ZZt_b_ii,
                         ZZt_eps,
                         ZZt_eps_ii,
                         c_hub,
                         c2_hub,
                         K2n,
                         trace = TRUE) {

  n_betas = ncol(Xn)
  npar = n_betas + 2
  if(npar != length(start))
    stop("Too many/much values in 'start'!")
  m_i = nrow(y)
  n_ind = ncol(y)
  n = nrow(Xn)
  hat_beta <- matrix(NA, n_iter, n_betas)
  hat_sig2 <- matrix(NA, n_iter, 2)
  hat_sig2[1, ] <- start[c(npar - 1, npar)]
  hat_beta[1, ] <- start[1:n_betas]
  Psi_hat <- matrix(NA, n_iter, npar)
  hat_fun <- rep(NA, n_iter)

  cat(
    '##################################################\n  Starting Gauss-Seidel iteration'
  )
  if (trace)
    cat("\nPsi: ")
  for (i in 2:n_iter) {
    # updating of beta
    opt.beta <- nlminb(
      c(hat_beta[i - 1, ]),
      objective = .beta_given_rest,
      rest = hat_sig2[i - 1, ],
      y = y,
      Xn = Xn,
      ZZt_b = ZZt_b,
      ZZt_b_ii = ZZt_b_ii,
      ZZt_eps = ZZt_eps,
      ZZt_eps_ii = ZZt_eps_ii,
      c_hub = c_hub,
      n_betas = n_betas,
      m_i = m_i,
      n = n,
      n_ind = n_ind
    )

    #updating of sigma2_b
    opt.sig2b <- uniroot(
      .sig2b_given_rest,
      interval = sig2b_interval,
      rest = c(opt.beta$par, hat_sig2[i - 1, 2]),
      y = y,
      Xn = Xn,
      Xnt = Xnt,
      ZZt_b = ZZt_b,
      ZZt_b_ii = ZZt_b_ii,
      ZZt_eps = ZZt_eps,
      ZZt_eps_ii = ZZt_eps_ii,
      c2_hub = c2_hub,
      K2n = K2n,
      n_betas = n_betas,
      m_i = m_i,
      n = n,
      n_ind = n_ind
    )

    #updating of sigma2_eps
    opt.sig2eps <- uniroot(
      .sig2eps_given_rest,
      interval = sig2eps_interval,
      rest = c(opt.beta$par, opt.sig2b$root),
      y = y,
      Xn = Xn,
      Xnt = Xnt,
      ZZt_b = ZZt_b,
      ZZt_b_ii = ZZt_b_ii,
      ZZt_eps = ZZt_eps,
      ZZt_eps_ii = ZZt_eps_ii,
      c2_hub = c2_hub,
      K2n = K2n,
      n_betas = n_betas,
      m_i = m_i,
      n = n,
      n_ind = n_ind
    )

    #update all and go the next iteration
    hat_beta[i, ] <- opt.beta$par
    hat_sig2[i, ] <- c(opt.sig2b$root, opt.sig2eps$root)
    Psi_hat[i, ] <- Psi_reml2(
      betas = hat_beta[i, ],
      sig2_b = hat_sig2[i, 1],
      sig2_eps = hat_sig2[i, 2],
      y = y,
      Xn = Xn,
      Xnt = Xnt,
      ZZt_b = ZZt_b,
      ZZt_b_ii = ZZt_b_ii,
      ZZt_eps = ZZt_eps,
      ZZt_eps_ii = ZZt_eps_ii,
      c_hub = c_hub,
      c2_hub = c2_hub,
      K2n = K2n,
      n_betas = n_betas,
      m_i = m_i,
      n = n,
      n_ind = n_ind
    )

    hat_fun[i] <- sqrt(crossprod(Psi_hat[i, ]))

    if (trace) {
      cat(paste(" ", hat_fun[i]))
    }
    if (i > 3) {
      cond1 <- hat_fun[i]
      cond3 <- hat_fun[i - 2]

      if (all.equal(cond1, cond3) == TRUE)
        break
    }
  }
  cat(
    '\n Finished after:',
    i,
    "iterations",
    '\n##################################################\n'
  )
  return(
    list(
      hat_beta = hat_beta[i, ],
      hat_sig2 = hat_sig2[i, ],
      hat_Psi = hat_fun[2:i],
      Psi_hat = Psi_hat[i, ]
    )
  )
}