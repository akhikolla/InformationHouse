#' ca_test_fit
#'
#' @param mcox An object as returned by  \code{agreg.fit}
#' @param X The \code{X} matrix of covariates as obtained from \code{model.matrix} without the cluster or intercept
#' @param atrisk A list with some fields from which the at-risk indicator can be deduced, as calculated in \code{emfrail()}
#' @param exp_g_x A vector of exponentiated linear predictor for each row of the data
#' @param cumhaz An estimate of the cumulative hazard at each time point in the data.
#'
#' @return A list with the test statistic, variance, and p-value
#'
#' @details This is an implementation of Commenges & Andersen (1995) test for heterogeneity, with a few adjustments to make it work with recurrent events data (?).
#' It could be made much faster and more efficient, but since it's not such an essential part, I'll let someone else do it.
#' @keywords internal
#'
#' @references Commenges, D. and Andersen, P.K., 1995. Score test of homogeneity for survival data. Lifetime Data Analysis, 1(2), pp.145-156.
#'
if(getRversion() >= "2.15.1")  utils::globalVariables(".")
ca_test_fit <- function(mcox, X, atrisk, exp_g_x, cumhaz) {

  # numerator

  S0_t <-
    lapply(seq_along(atrisk$time), function(x)
      exp_g_x[atrisk$indx2 < x & x <= atrisk$time_to_stop ]) %>%
    lapply(sum) %>%
    do.call(c, .)
  #
  #
  # S0 <- (nrisk/newrisk)[atrisk$time_to_stop]
  #
  #
  #

  # the basic ingredients
  # pij_t is a list of length (tev) where at each event time point we have the sum of
  # elp of people at risk at that time point
  # all divided by S0_t
  pij_t <-
    lapply(seq_along(atrisk$time), function(x)
      exp_g_x * as.numeric(atrisk$indx2 < x & x <= atrisk$time_to_stop )) %>%
    lapply(as.numeric) %>%
    mapply(function(a,b) a/ b, ., S0_t, SIMPLIFY = FALSE)

  pi_t <- lapply(pij_t, function(x) rowsum(x, atrisk$order_id)) %>%
    lapply(as.numeric)

  T_stat <- sum(rowsum(mcox$residuals , atrisk$order_id)  ^2) -
    sum(atrisk$death) +
  sum(pi_t %>%
    lapply(function(x) sum(x^2)) %>%
    do.call(c, .) *
    atrisk$nevent )


  # this is the martingale path of each line

  mt_ij <- mapply(function(pos_left, pos_right) {

    if(pos_left == 0) {
      cumhaz[pos_right:length(cumhaz)] <- cumhaz[pos_right]
    return(-cumhaz)}


    cumhaz[(pos_left+1):length(cumhaz)] <- cumhaz[(pos_left+1):length(cumhaz)] - cumhaz[pos_left]
    cumhaz[1:(pos_left+1)] <- 0
    cumhaz[pos_right:length(cumhaz)] <- cumhaz[pos_right]

    return(-cumhaz)
    # stop it after pos_right

  }, atrisk$indx2, atrisk$time_to_stop, SIMPLIFY = FALSE) %>%
    mapply(function(a,b) a*b, ., exp_g_x, SIMPLIFY = FALSE) %>%
    mapply(function(a, pos, delta) {
      a[pos:length(a)] <- delta + a[pos]
      a
    }, ., atrisk$time_to_stop, atrisk$death, SIMPLIFY  = FALSE)


  m0t_ij <-
    lapply(mt_ij, function(x) c(0, x[-length(x)]))

  mij_t <- m0t_ij %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    as.list()

  mi_t <- lapply(mij_t, function(x) rowsum(x, atrisk$order_id)) %>%
    lapply(as.numeric)


  mp_t <- mapply(function(a,b) a * b, mi_t, pi_t, SIMPLIFY = FALSE) %>%
    lapply(sum)

  pp_t <- pi_t %>%
    lapply(function(x) x^2)  %>%
    lapply(sum)


  qi_t <- mapply(function(a,b,c,d) 2 * (a - b - c + d), mi_t, mp_t, pi_t, pp_t, SIMPLIFY = FALSE)

  # Main part of V
  V1 <- qi_t %>%
    lapply(function(x) x^2) %>%
    mapply(function(a,b,c) a * b * c, ., pi_t, atrisk$nevent, SIMPLIFY = FALSE) %>%
    do.call(sum, .)

  # second part of V
  if(!is.null(mcox$var)) {
  theta2i_t <- pij_t %>%
    lapply(function(x) x * X) %>%
    lapply(as.data.frame) %>%
    lapply(function(x) rowsum.data.frame(x, atrisk$order_id))

  theta_h <- atrisk$nevent %>%
    mapply(function(a,b,c) a * b * c, ., qi_t, theta2i_t, SIMPLIFY = FALSE) %>%
    do.call(rbind, .) %>%
    apply(., 2, sum)

  V2 <- t(theta_h) %*% mcox$var %*% theta_h
  } else V2 <- 0



  V <- V1 - V2

  c(tstat = T_stat, var = V, pval = pchisq(T_stat^2 / V, 1, lower.tail = FALSE))
}

