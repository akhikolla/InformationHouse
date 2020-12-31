#' Compute middle p-value
#'
#' Computes P(Y > cases) + P(Y = cases)/2 when Y ~
#' Poisson(ex) or Y ~ Binomial(n = pop, p = ex/pop).  This
#' is middle p-value computed by Tango and Takahashi (2012).
#'
#' @inheritParams scan.test
#' @return A vector of middle p-values
#' @export
#' @author Joshua French
#' @references Tango, T. and Takahashi, K. (2012), A
#'   flexible spatial scan statistic with a restricted
#'   likelihood ratio for detecting disease clusters.
#'   Statist. Med., 31: 4207-4218. <doi:10.1002/sim.5478>
#'
#' @examples
#' data(nydf)
#' cases = floor(nydf$cases)
#' pop = nydf$pop
#' ex = pop * sum(cases)/sum(pop)
#' # zones for poisson model
#' pp = rflex.midp(cases, ex)
#' # zones for binomial model
#' bp = rflex.midp(cases, ex, type = "binomial", pop = pop)
rflex.midp = function(cases, ex, type = "poisson", pop = NULL) {
  arg_check_rflex_midp(cases, ex, type, pop)

  if (type == "poisson") {
    p = stats::ppois(cases, ex, lower.tail = FALSE) +
      stats::dpois(cases, ex) / 2
  } else if (type == "binomial") {
    p = stats::pbinom(cases, size = pop, prob = ex / pop,
                      lower.tail = FALSE) +
      stats::dbinom(cases, size = pop, prob = ex / pop) / 2
  }
  return(p)
}

arg_check_rflex_midp = function(cases, ex, type, pop) {
  if (length(cases) != length(ex)) {
    stop("length(cases) != length(ex)")
  }
  if (!is.numeric(cases) | !is.numeric(ex)) {
    stop("cases and ex must be numeric vectors")
  }
  if (length(type) != 1 |
      !is.element(type, c("poisson", "binomial"))) {
    stop("type must be poisson or binomial")
  }
  if (type == "binomial" & is.null(pop)) {
    stop("pop must be provided when type = 'binomial'")
  }
  if (!is.null(pop)) {
    if (length(cases) != length(pop)) {
      stop("length(case) != length(pop)")
    }
  }
}
