#' Power Calculations for Two-Sample Test for Proportions with unequal sample size
#'
#' Compute power of test, or determine parameters to obtain target
#' power for equal and unequal sample sizes.
#'
#' @param n Number of observations (in group 1)
#' @param p1 Probability in one group
#' @param p2 Probability in other group
#' @param sig.level Significance level (Type I error probability)
#' @param power Power of test (1 minus Type II error probability)
#' @param ratio The ratio n2/n1 between the larger group and the smaller group. Should be a value equal to or greater than 1 since n2 is the larger group. Defaults to 1 (equal group sizes)                                                                                            
#' @param alternative String. Can be one- or two-sided test. Can be abbreviated.
#' @param tol Numerical tolerance used in root finding, the default providing (at least) four significant digits
#' @return Object of class \code{power.htest}, a list of the arguments (including the computed one) augmented with \code{method} and \code{note} elements.                                                                                    
#' @details Exactly one of the parameters \code{n}, \code{delta}, \code{power}, \code{sd}, \code{sig.level}, \code{ratio} \code{sd.ratio}
#' must be passed as NULL, and that parameter is determined from the others. Notice that the last two have non-NULL defaults
#' so NULL must be explicitly passed if you want to compute them.                        
#' @note \code{uniroot} is used to solve power equation for unknowns, so you may
#' see errors from it, notably about inability to bracket the root
#' when invalid arguments are given.
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{\link{power.prop.test}}, \code{\link{power_t_test}}, \code{\link{power.t.test}}
#' @keywords htest
#' @examples
#' power_prop_test(n=NULL, p1=.65, p2=.85, power=.8, ratio=2)
#' @export
power_prop_test <- function (n = NULL, p1 = NULL, p2 = NULL, sig.level = 0.05, power = NULL, ratio=1,
    alternative = c("two.sided", "one.sided"), 
    tol = .Machine$double.eps^0.25) 
{
    if (sum(sapply(list(n, p1, p2, power, sig.level, ratio), is.null)) != 1) 
        stop("exactly one of 'n', 'p1', 'p2', 'power', 'sig.level', and 'ratio' must be NULL")
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
        stop("'sig.level' must be numeric in [0, 1]")
    if (!is.null(ratio) && ratio <= 0)
        stop("ratio between group sizes must be positive")
    
    alternative <- match.arg(alternative)
    tside <- switch(alternative, one.sided = 1, two.sided = 2)
    p.body <- quote({
            qu <- qnorm(sig.level/tside, lower.tail = FALSE)
            d <- abs(p1 - p2)
            q1 <- 1 - p1
            q2 <- 1 - p2
            pbar <- (p1 + ratio*p2)/(1+ratio)
            qbar <- 1 - pbar

            (  qnorm(sig.level/tside)* sqrt(pbar*qbar*(1+ratio))  + qnorm(1-power)* sqrt(ratio*p1*q1 + p2*q2) )^2  / (ratio*d^2)

        })
    if (is.null(n)) 
        n <- eval(p.body)
    else if (is.null(power)) 
        power <- uniroot(function(power) eval(p.body) - n, c(0.00001, .99999), 
            tol = tol, extendInt = "upX")$root    
    else if (is.null(p1)) 
        p1 <- uniroot(function(p1) eval(p.body) - n, c(0, 
            p2), tol = tol, extendInt = "yes")$root
    else if (is.null(p2)) 
        p2 <- uniroot(function(p2) eval(p.body) - n, c(p1, 
            1), tol = tol, extendInt = "yes")$root
    else if (is.null(ratio))
        ratio <- uniroot(function(ratio) eval(p.body) - n, c(2/n, 1e+07))$root
    else if (is.null(sig.level)) 
        sig.level <- uniroot(function(sig.level) eval(p.body) - 
            n, c(1e-10, 1 - 1e-10), tol = tol, extendInt = "upX")$root
    else stop("internal error", domain = NA)

    if (ratio!=1) {
      n <- c(n, n*ratio)
    }

    NOTE <- ifelse(ratio==1, "n is number in *each* group", "n is vector of number in each group")

    METHOD <- ifelse(ratio==1,
                     "Two-sample comparison of proportions power calculation",
                     "Two-sample comparison of proportions power calculation with unequal sample sizes")
    structure(list(n = n, p1 = p1, p2 = p2, sig.level = sig.level, 
        power = power, alternative = alternative, note = NOTE, 
        method = METHOD), class = "power.htest")
}
