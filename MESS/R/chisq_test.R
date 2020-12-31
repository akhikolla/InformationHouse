#' Two-sided table test with fixed margins
#'
#' @description Monte Carlo test in a two-way contingency table with the total number of observations fixed, row margin fixed, or both margins fixed. 
#' @param x A matrix representing the contingency table.
#' @param margin A string that determines which margin is fixed: Either "N" for the total number of observations (the default), "rows" for fixed row sums, and "both" for simultaneously fixed row and column sums.
#' @param B The number of simulations used to compute the p-value.
#' @details Simulation is done by random sampling from the set of all tables with given marginal(s), and works only if the relevant marginal(s) are strictly positive. Continuity correction is never used, and the statistic is quoted without it.
#' @return A list of class "htest" giving the simulation results.
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @keywords manip
#' @examples
#'
#' m <- matrix(c(12, 4, 8, 6), 2)
#' chisq.test(m)
#' chisq.test(m, correct=FALSE)
#' monte_carlo_chisq_test(m)
#' 
#' fisher.test(m)
#' monte_carlo_chisq_test(m, margin="both")
#'
#' m2 <- matrix(c(9, 3, 3, 7), 2)
#' monte_carlo_chisq_test(m, margin="N")
#' monte_carlo_chisq_test(m, margin="both") 
#' 
#' @export
monte_carlo_chisq_test <- function(x,
                                   margin=c("N", "rows", "both"),
                                   B = 100000L) {
    
    DNAME <- deparse(substitute(x))
    margin <- match.arg(margin)
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    if (is.matrix(x)) {
        if (min(dim(x)) == 1L)
            stop("Only works for 2-way tables at the moment")
    }
    if (!is.matrix(x)) {
        stop("x must be convertible with a matrix")
    }

    ## Sanity checks
    if (any(x < 0) || anyNA(x)) 
        stop("all entries of 'x' must be nonnegative and finite")
    if ((n <- sum(x)) == 0) 
        stop("at least one entry of 'x' must be positive")
 
    METHOD <- "Monte Carlo Pearson's Chi-squared test"

    rsum <- rowSums(x)
    csum <- colSums(x)
    expected <- outer(rsum, csum)/sum(rsum)            


    marg <- 0
    fixed = "N (total number of observations)"
    if (margin=="both") {
        ## Do something special for Fisher's test situation
        fixed <- "row AND column margins"

        ## Check that we do not have any zeros in the marginals
        if (any(rsum<=0))
            stop("All row marginals must be non-negative")
        if (any(csum<=0))
            stop("All column marginals must be non-negative")

        ## Use the simulate p value argument from chisq.test
        res <- chisq.test(x, correct=FALSE, B=B, simulate.p.value=TRUE)
    } else {
        
        if (margin=="rows") {
            if (any(rsum<=0))
                stop("All row marginals must be non-negative")
            marg <- 1
            fixed = "row margins"
            res <- .chisq_test_cpp(x, margin=marg, B=B);            
        } else if (margin=="N") {
            res <- .chisq_test_cpp(x, margin=marg, B=B);            
        }        
    }

    names(B) <- "B"
    
    METHOD <- paste(METHOD, " \n\t (based on ", B, " replicates with fixed ", fixed, ")", sep="")
    
    structure(list(statistic = res$statistic,
                   parameter = B, 
                   p.value = res$p.value,
                   method = METHOD,
                   data.name = DNAME,
                   observed = x, 
                   expected = expected),
              class = "htest")
    
    
}
