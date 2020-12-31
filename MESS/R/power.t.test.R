#' Power calculations for one and two sample t tests with unequal sample size
#'
#' Compute power of test, or determine parameters to obtain target
#' power for equal and unequal sample sizes.
#'
#' @param n Number of observations (in the smallest group if two groups)
#' @param delta True difference in means
#' @param sd Standard deviation
#' @param sig.level Significance level (Type I error probability)
#' @param power Power of test (1 minus Type II error probability)
#' @param ratio The ratio n2/n1 between the larger group and the smaller group. Should be a value equal to or greater than 1 since n2 is the larger group. Defaults to 1 (equal group sizes). If ratio is set to NULL (i.e., find the ratio) then the ratio might be smaller than 1 depending on the desired power and ratio of the sd's.
#' @param sd.ratio The ratio sd2/sd1 between the standard deviations in the larger group and the smaller group. Defaults to 1 (equal standard deviations in the two groups)
#' @param type Type of t test
#' @param alternative One- or two-sided test
#' @param df.method Method for calculating the degrees of default. Possibilities are welch (the default) or classical.
#' @param strict Use strict interpretation in two-sided case. Defaults to TRUE unlike the standard power.t.test function.
#'
#' @return  Object of class \code{power.htest}, a list of the arguments (including the computed one)
#' augmented with \code{method} and \code{note} elements.
#' @details Exactly one of the parameters \code{n}, \code{delta}, \code{power}, \code{sd}, \code{sig.level}, \code{ratio} \code{sd.ratio}
#' must be passed as NULL,
#' and that parameter is determined from the others. Notice that the last two have non-NULL defaults
#' so NULL must be explicitly passed if you want to compute them.
#'
#' The default \code{strict = TRUE} ensures that the power will include the probability
#' of rejection in the opposite direction of the true effect, in the
#' two-sided case. Without this the power will be half the
#' significance level if the true difference is zero.
#' @note \code{uniroot} is used to solve power equation for unknowns, so you may
#' see errors from it, notably about inability to bracket the root
#' when invalid arguments are given.
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{\link{power.t.test}}, \code{\link{power_prop_test}}, \code{\link{power.prop.test}}
#' @keywords htest
#' @examples
#' # Sampling with a ratio of 1:4
#' power_t_test(delta=300, sd=450, power=.8, ratio=4)
#'
#' # Equal group sizes but different sd's
#' # The sd in the first group is twice the sd in the second group
#' power_t_test(delta=300, sd=450, power=.8, sd.ratio=.5)
#'
#' # Fixed group one size to 50 individuals, but looking for the number of individuals in the
#' # second group. Different sd's with twice the sd in the larger group
#' power_t_test(n=50, delta=300, sd=450, power=.8, ratio=NULL, sd.ratio=2)
#' @export
power_t_test <-
  function (n = NULL, delta = NULL, sd = 1, sig.level = 0.05, power = NULL,
            ratio = 1, sd.ratio = 1,
            type = c("two.sample", "one.sample", "paired"),
            alternative = c("two.sided", "one.sided"),
            df.method = c("welch", "classical"),
            strict = TRUE)
{
  type <- match.arg(type)
  if (type == "two.sample") {

    if (sum(sapply(list(n, delta, sd, power, sig.level, ratio, sd.ratio), is.null)) != 1)
      stop("exactly one of n, delta, sd, power, sig.level, ratio and sd.ratio must be NULL")

    if (!is.null(ratio) && ratio < 1)
      stop("ratio between group sizes cannot be less than 1")
    if (!is.null(sd.ratio) && sd.ratio < 0)
      stop("sd.ratio between group sd's cannot be less than 1")
  }
  else {
      ratio <- 1
      sd.ratio <- 1
      if (sum(sapply(list(n, delta, sd, power, sig.level), is.null)) != 1)
          stop("exactly one of n, delta, sd, power, and sig.level must be NULL")
  }

  alternative <- match.arg(alternative)
  df.method <- match.arg(df.method)
  tsample <- switch(type, one.sample = 1, two.sample = 2, paired = 1)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (tside == 2 && !is.null(delta))
    delta <- abs(delta)
  p.body <- quote({
    nu <- switch(tsample, n-1, switch(df.method, welch=(sd^2/n + (sd*sd.ratio)^2/(n*ratio))^2/((sd^2/n)^2/(n-1) + ((sd*sd.ratio)^2/(ratio*n))^2/(n*ratio-1)),
classical=(1+ratio)*n-2))
    pt(qt(sig.level/tside, nu, lower = FALSE), nu, ncp = switch(tsample, sqrt(n/tsample), sqrt(n/(1+sd.ratio^2/ratio))) * delta/sd, lower = FALSE)
  })
  if (strict & tside == 2)
    p.body <- quote({
      nu <- switch(tsample, n-1, switch(df.method, welch=(sd^2/n + (sd*sd.ratio)^2/(n*ratio))^2/((sd^2/n)^2/(n-1) + ((sd*sd.ratio)^2/(ratio*n))^2/(n*ratio-1)),
classical=(1+ratio)*n-2))
      qu <- qt(sig.level/tside, nu, lower = FALSE)
      ncp <- switch(tsample, sqrt(n/tsample), sqrt(n/(1+sd.ratio^2/ratio))) * delta/sd
      pt(qu, nu, ncp = ncp, lower = FALSE) +
        pt(-qu, nu, ncp = ncp, lower = TRUE)
    })
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(sd))
    sd <- uniroot(function(sd) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body) - power, sd * c(1e-07, 1e+07))$root
  else if (is.null(sig.level))
    sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else if (is.null(ratio))
    ratio <- uniroot(function(ratio) eval(p.body) - power, c(2/n, 1e+07))$root
  else if (is.null(sd.ratio))
    sd.ratio <- uniroot(function(sd.ratio) eval(p.body) - power, c(1e-07, 1e+07))$root
  else stop("internal error")
  NOTE <- switch(type,
                 paired = "n is number of *pairs*, sd is std.dev. of *differences* within pairs",
                 two.sample = ifelse(ratio==1, "n is number in *each* group", "n is vector of number in each group"),
                 NULL)
#  n <- switch(type, paired=n, two.sample=c(n, ifelse(ratio==1, NULL, n*ratio)), one.sample=n)
#  sd <- switch(type, paired=sd, two.sample=c(sd, ifelse(ratio==1, NULL, sd*sd.ratio)), one.sample=sd)

  if (type=="two.sample" & (ratio!=1 | sd.ratio !=1)) {
      n <- c(n, n*ratio)
      sd <- c(sd, sd*sd.ratio)
  }


  METHOD <- paste(switch(type, one.sample = "One-sample t test power calculation",
                               two.sample = ifelse(ratio==1, "Two-sample t test power calculation", "Two-sample t test power calculation with unequal sample sizes"),
                         paired = "Paired t test power calculation"))
  if (type=="two.sample" & sd.ratio != 1) {
      METHOD <- paste0(METHOD, ifelse(ratio==1, " with", " and"), " unequal variances")
  }
  structure(list(n = n, delta = delta, sd = sd, sig.level = sig.level,
                 power = power, alternative = alternative, note = NOTE,
                 method = METHOD), class = "power.htest")
}
