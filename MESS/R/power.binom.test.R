#' Power Calculations for Exact Test of a simple null hypothesis in a Bernoulli
#' experiment
#'
#' Compute power of test, or determine parameters to obtain target power.
#'
#' The procedure uses uniroot to find the root of a discontinuous function so
#' some errors may pop up due to the given setup that causes the root-finding
#' procedure to fail. Also, since exact binomial tests are used we have
#' discontinuities in the function that we use to find the root of but despite
#' this the function is usually quite stable.
#'
#' @param n Number of observations
#' @param p0 Probability under the null
#' @param pa Probability under the alternative
#' @param sig.level Significance level (Type I error probability)
#' @param power Power of test (1 minus Type II error probability)
#' @param alternative One- or two-sided test
#' @return Object of class \code{power.htest}, a list of the arguments
#' (including the computed one) augmented with method and note elements.
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{\link{binom.test}}
#' @keywords htest
#' @examples
#'
#' power_binom_test(n = 50, p0 = .50, pa = .75)      ## => power = 0.971
#' power_binom_test(p0 = .50, pa = .75, power = .90) ## =>     n = 41
#' power_binom_test(n = 50, p0 = .25, power = .90, alternative="less")  ## => pa = 0.0954
#'
#' @export
power_binom_test <- function(n = NULL, p0 = NULL, pa = NULL, sig.level = 0.05, power = NULL,
                               alternative = c("two.sided", "less", "greater")) {

    if (sum(sapply(list(n, p0, pa, power, sig.level), is.null)) != 1)
        stop("exactly one of 'n', 'p0', 'pa', 'power', and 'sig.level' must be NULL")
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 >
        sig.level | sig.level > 1))
        stop("'sig.level' must be numeric in [0, 1]")

    alternative <- match.arg(alternative)

    # Check boundary conditions!
    # Think they should be okay

    pfun <- function(n, p0, pa, sig.level, alternative) {

        n <- ceiling(n)

    power <- switch(alternative,
        less = { # Slightly counterintuitive to get the border with
            pbinom(qbinom(1-sig.level, size=n, prob=p0, lower.tail=FALSE)-1, size=n, prob=pa)
        },
        greater = {
            pbinom(qbinom(1-sig.level, size=n, prob=p0), size=n, prob=pa, lower.tail=FALSE)
        },
        two.sided = {
            # Stat by finding a large set of values that do contain the proper quantile/probability where
            # These are intentionally a bit too large to reduce computations but still make sure all necessary values are in the set
            lx <- qbinom(sig.level, size=n, prob=p0)
            ux <- qbinom(sig.level, size=n, prob=p0, lower.tail=FALSE)

            x <- c(seq(0,lx), seq(ux,n))
            d <- dbinom(x, size=n, prob=p0)
            # Order to rank the probabilities accodring to size (small to high)
            ordd <- order(d)
            # Calculate cumulative probabilities
            cs <- cumsum(sort(d))
            # Find position in x where the significance level is too much
            xval <- which.min(cs<sig.level)-1
            #cat("asdasd")
            #print(ordd)
            #print(c(xval, ordd[xval]))
            # This position corresponds to this probability
            ssh <- d[ordd[xval]]
            #cat("Order of the cumsum match")
            #print(ordd[xval])
            #cat("SSH to match")
            #print(ssh)

            relErr <- 1 + 1e-07

            m <- n*p0

            if (xval==0)
                return(0)

            if (x[ordd[xval]] < m) {
                i <- seq.int(from = ux, to = n)
                y <- sum(dbinom(i, n, p0) <= ssh * relErr)
                pbinom(x[ordd[xval]], size=n, prob=pa) + pbinom(n - y, size=n, prob=pa, lower.tail = FALSE)
            } else {
                i <- seq.int(from = 0, to = lx)
                y <- sum(dbinom(i, n, p0) <= ssh * relErr)
                pbinom(y-1, size=n, prob=pa) + pbinom(x[ordd[xval]]-1, n, pa, lower.tail = FALSE)
            }
        }
        )
    power
    }

#    p.body <- Vectorize(pfun)
#    ppp <- body(p.body)

    qqq <- quote({ do.call("mapply", c(FUN=pfun, list(n, p0, pa, sig.level, alternative),  SIMPLIFY = TRUE,
       USE.NAMES = TRUE)) })

    if (is.null(power))
        power <- eval(qqq)
    else if (is.null(n)) {
        ans <- uniroot(function(n) eval(qqq) - power, c(2, 1e+06))
        # Make small correction for shooting low to ensure that the power is at least as desired
        n <- ans$root + (ans$f.root<0)
    }
    else if (is.null(p0)) {
        p0 <- switch(alternative,
                     less = {
                         uniroot(function(p0) eval(qqq) - power, c(1e-07, p0-1e-07))$root
                     },
                     greater = {
                         uniroot(function(p0) eval(qqq) - power, c(p0+1e-07, 1-1e-07))$root
                     },
                     two.sided = {
                         stop("Not making a lot of sense without some assumptions of symmetry or something")
                         uniroot(function(p0) eval(qqq) - power, c(1e-07, 1-1e-07))$root
                     }
                     )
    }
    else if (is.null(pa)) {
        pa <- switch(alternative,
                     less = {
                         uniroot(function(pa) eval(qqq) - power, c(1e-07, p0-1e-07))$root
                     },
                     greater = {
                         uniroot(function(pa) eval(qqq) - power, c(p0+1e-07, 1-1e-07))$root
                     },
                     two.sided = {
                         stop("Not making a lot of sense without some assumptions of symmetry or something")
                         uniroot(function(pa) eval(qqq) - power, c(1e-07, 1-1e-07))$root
                     }
                     )
    }
    else if (is.null(sig.level)) {
        sig.level <- switch(alternative,
                     less = {
                         uniroot(function(sig.level) eval(qqq) - power, c(1e-07, 1-1e-07))$root
                     },
                     greater = {
                         uniroot(function(sig.level) eval(qqq) - power, c(1e-07, 1-1e-07))$root
                     },
                     two.sided = {
                         stop("Not making a lot of sense without some assumptions of symmetry or something")
                         uniroot(function(sig.level) eval(qqq) - power, c(1e-07, 1-1e-07))$root
                     }
                     )
    }
    else stop("internal error", domain = NA)

    NOTE <- NULL
    METHOD <- "One-sample exact binomial power calculation"
    structure(list(n = n, p0 = p0, pa = pa, sig.level = sig.level,
        power = power, alternative = alternative, note = NOTE,
        method = METHOD), class = "power.htest")
}


