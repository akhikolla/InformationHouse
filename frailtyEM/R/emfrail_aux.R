#' dist_to_pars
#'
#' @param dist One of gamma, stable, pvf
#' @param logfrailtypar The log of theta
#' @param pvfm The pvfm
#'
#' @keywords internal
#' @return A list with 3 elements: alpha, beta (the parameters of the Laplace transform) and dist_id.
dist_to_pars <- function(dist, logfrailtypar, pvfm) {

    if (dist == "gamma") {
        alpha <- bbeta <- exp(logfrailtypar)
        dist_id <- 0L
    }

    # if (dist == "stable") {
    #     theta <- exp(logfrailtypar) + 1  # so theta >1
    #     bbeta <- 1 - 1/theta  # so bbeta in (0,1), that's what's important
    #     alpha <- theta / (theta - 1)  # alpha = 1/beta for scaling
    #     dist_id <- 1L
    # }

    if (dist == "stable") {
      # theta <- exp(logfrailtypar) + 1 # so theta >1
      # bbeta <- 1 - 1/theta
      alpha <- 1
      #bbeta <- 1 - exp(logfrailtypar) / (exp(logfrailtypar) + 1)
      bbeta <- exp(logfrailtypar) / (exp(logfrailtypar) + 1)
      dist_id <- 1L
    }

    if (dist == "pvf") {
        alpha <- abs((pvfm + 1)/pvfm * exp(logfrailtypar))
        bbeta <- (pvfm + 1) * exp(logfrailtypar)
        dist_id <- 2L
    }

    list(alpha = alpha, bbeta = bbeta, dist = dist_id)
}


#' Laplace transform calculation
#'
#' @param x A vector of positive values where to calculate the Laplace transform
#' @param distribution An \code{emfrail_dist} object. See \code{?emfrail_dist}.
#'
#' @return A vector of the same length as \code{x} with the Laplace transform of \code{x}
#'
#' @details This is a simple function which calculates the Laplace transform for the gamma, positive stable and PVF distribution.
#' It is intended to be used to calculate marginal quantities from an \code{emfrail} object.
#' Note that the \code{left_truncation} argument is ignored here;
#' the marginal survival or hazard are given for the Laplace transform of a baseline subject entered at time 0.
#' @keywords internal
laplace_transform <- function(x, distribution) {
  # if(missing(.distribution) & missing())
  if(!inherits(distribution, "emfrail_dist"))
    stop("distribution argument misspecified; see ?emfrail_dist()")

  getpars <- dist_to_pars(distribution$dist, log(distribution$frailtypar), distribution$pvfm)

  if(getpars$dist == 0L) {
    L <- with(getpars, (bbeta / (bbeta + x))^alpha)
  }

  if(getpars$dist == 1L) {
    L <- with(getpars, exp(-1 * x^bbeta))
  }

  if(getpars$dist == 2L) {
    L <- with(getpars, exp(-alpha * sign(distribution$pvfm) * (1 - (bbeta / (bbeta + x))^distribution$pvfm )))
  }

  L

}


#' Profile log-likelihood calculation
#'
#' @param data Same as in \code{emfrail}
#' @param formula Same as in \code{emfrail}
#' @param distribution Same as in \code{emfrail}
#' @param values A vector of values on where to calculate the profile likelihood. See details.
#'
#' @return The profile log-likelihood at the specific value of the frailty parameter
#' @export
#'
#' @details This function can be used to calculate the profile log-likelihood for different values of \eqn{\theta}.
#' The scale is that of \code{theta} as defined in \code{emfrail_dist()}.
#' For the gamma and pvf frailty, that is the inverse of the frailty variance.
#'
#' @note This function is just a simple wrapper for \code{emfrail()} with the \code{control} argument
#' a call from \code{emfrail_control} with the option \code{opt_fit = FALSE}. More flexibility can be obtained
#' by calling \code{emfrail} with this option, especially
#' for setting other \code{emfrail_control} parameters.
#'
#' @examples
#'
#' fr_var <- seq(from = 0.01, to = 1.4, length.out = 20)
#' pll_gamma <- emfrail_pll(formula = Surv(time, status) ~  rx + sex + cluster(litter),
#'  data =  rats,
#'  values = 1/fr_var )
#'  plot(fr_var, pll_gamma,
#'      type = "l",
#'      xlab = "Frailty variance",
#'      ylab = "Profile log-likelihood")
#'
#' # check with coxph;
#' # attention: theta is the the inverse frailty variance in emfrail,
#' # but theta is the frailty variance in coxph.
#'
#' pll_cph <- sapply(fr_var, function(th)
#'   coxph(data =  rats, formula = Surv(time, status) ~ rx + sex + frailty(litter, theta = th),
#'         method = "breslow")$history[[1]][[3]])
#'
#' lines(fr_var, pll_cph, col = 2)
#'
#' # Same for inverse gaussian
#' pll_if <- emfrail_pll(Surv(time, status) ~  rx + sex + cluster(litter),
#'                       rats,
#'                       distribution = emfrail_dist(dist = "pvf"),
#'                       values = 1/fr_var )
#'
#' # Same for pvf with a positive pvfm parameter
#' pll_pvf <- emfrail_pll(Surv(time, status) ~  rx + sex + cluster(litter),
#'                        rats,
#'                        distribution = emfrail_dist(dist = "pvf", pvfm = 1.5),
#'                        values = 1/fr_var )
#'
#' miny <- min(c(pll_gamma, pll_cph, pll_if, pll_pvf))
#' maxy <- max(c(pll_gamma, pll_cph, pll_if, pll_pvf))
#'
#' plot(fr_var, pll_gamma,
#'      type = "l",
#'      xlab = "Frailty variance",
#'      ylab = "Profile log-likelihood",
#'      ylim = c(miny, maxy))
#' points(fr_var, pll_cph, col = 2)
#' lines(fr_var, pll_if, col = 3)
#' lines(fr_var, pll_pvf, col = 4)
#'
#' legend(legend = c("gamma (emfrail)", "gamma (coxph)", "inverse gaussian", "pvf, m=1.5"),
#'        col = 1:4,
#'        lty = 1,
#'        x = 0,
#'        y = (maxy + miny)/2)

emfrail_pll <- function(formula, data,
                        distribution = emfrail_dist(),
                        values) {
  sapply(values, function(fp) {
    -emfrail(formula = formula,
             data = data,
             distribution =  emfrail_dist(dist = distribution$dist,
                                                  theta = fp,
                                                  pvfm = distribution$pvfm,
                                                  left_truncation = distribution$left_truncation),
             control = emfrail_control(opt_fit = FALSE))
  })



}



# this one gives the baseline cumulative hazard at all the time points;
#
# getchz <- function(Y, newrisk, explp) {
#     death <- (Y[, ncol(Y)] == 1) # this is a TRUE FALSE for the status column
#     dtime <- Y[, ncol(Y) - 1] # this is the tstop
#
#     time <- sort(unique(dtime)) # unique tstops
#
#     nevent <- as.vector(rowsum(1 * death, dtime))
#
#     nrisk <- rev(cumsum(rev(rowsum(explp, dtime)))) # This gives the sum
#     delta <- min(diff(time))/2
#     etime <- c(sort(unique(Y[, 1])), max(Y[, 1]) + delta)  #unique entry times
#
#     indx <- approx(etime, 1:length(etime), time, method = "constant", rule = 2, f = 1)$y
#
#     esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))  #not yet entered
#     nrisk <- nrisk - c(esum, 0)[indx]
#
#     haz <- nevent/nrisk
#     cumhaz <- cumsum(haz)
#
#     chz2 <- cumhaz * newrisk
#
#     tev = time[haz > 0]
#     haz_ret = haz * newrisk
#
#     haz_tev = haz_ret[haz_ret > 0]
#
#
#     list(time = time, cumhaz = chz2, haz = haz_ret, tev = tev, haz_tev = haz_tev)
#
#
# }
#
#
