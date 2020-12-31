#' Goodman-Kruskal's gamma statistic for a two-dimensional table
#'
#' Compute Goodman-Kruskal's gamma statistic for a two-dimensional table of
#' ordered categories
#'
#'
#' @param x A matrix or table representing the two-dimensional ordered
#' contingency table of observations
#' @param conf.level Level of confidence interval
#' @return A list with class \code{htest} containing the following components:
#'
#' \item{statistic }{the value the test statistic for testing no association}
#' \item{p.value }{the p-value for the test} \item{estimate }{the value the
#' gamma estimate} \item{conf.int }{the confidence interval for the gamma
#' estimate} \item{method }{a character string indicating the type of test
#' performed} \item{data.name }{a character string indicating the name of the
#' data input} \item{observed }{the observed counts} \item{s0 }{the SE used
#' when computing the test statistics} \item{s1 }{the SE used when computing
#' the confidence interval}
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{\link{chisq.test}}
#' @references Goodman, Leo A. and Kruskal, William H. (1954). "Measures of
#' Association for Cross Classifications". Journal of the American Statistical
#' Association 49 (268): 732-764.
#' @keywords htest
#' @examples
#'
#' # Data from the Glostrup study comparing smoking to overall health in males
#' smoke <- matrix(c(16, 15, 13, 10, 1, 73, 75, 59, 81, 29, 6, 6, 7, 17, 3, 1, 0, 1, 3, 1), ncol=4)
#' colnames(smoke) <- c("VGood", "Good", "Fair", "Bad") # General health status
#' rownames(smoke) <- c("Never", "No more", "1-14", "15-24", "25+")  # Smoke amount
#' gkgamma(smoke)
#' chisq.test(smoke)
#'
#' @export gkgamma
gkgamma <- function(x, conf.level = 0.95) {

    if (! (is.matrix(as.matrix(x)))) {
        stop("x must be a table or matrix")
    }

    DNAME <- deparse(substitute(x))

    x <- as.matrix(x)

    rows <- nrow(x)
    cols <- ncol(x)
    n <- sum(x)

    con <- x
    dis <- x

    rseq <- 1:rows
    cseq <- 1:cols

    for (i in 1:rows) {
        for (j in 1:cols) {
            con[i,j] <- sum(x[rseq<i,cseq<j]) + sum(x[rseq>i,cseq>j])
            dis[i,j] <- sum(x[rseq>i,cseq<j]) + sum(x[rseq<i,cseq>j])
        }
    }

    CC <- sum(x*con)
    DC <- sum(x*dis)

    ESTIMATE <- (CC-DC)/(CC+DC)

    se1 <- sqrt(sum(x*(DC*con - CC*dis)^2))*4/((CC+DC)^2) # sqrt(sum(x*(DC*con - CC*dis)^2)*16/(CC+DC)^4)
    se0 <- 2 / (CC + DC) * sqrt(sum(x*(con-dis)^2) - (CC-DC)^2/n)

    STATISTIC <- ESTIMATE / se0

    PVAL <- 2*pnorm(-abs(STATISTIC))
    CINT <- ESTIMATE + c(1,-1)*qnorm((1-conf.level)/2)*se1
    if (!is.null(CINT))
        attr(CINT, "conf.level") <- conf.level

    METHOD <- "Goodman-Kruskal's gamma for ordinal categorical data"
    names(STATISTIC) <- "Z"
    names(ESTIMATE) <- "Goodman-Kruskal's gamma"
                                        #  names(CINT) <- "Confidence "
                                        #  names(PARAMETER) <- "df"
                                        #  if (any(E < 5) && is.finite(PARAMETER))
                                        #    warning("Chi-squared approximation may be incorrect")
    structure(list(statistic = STATISTIC,
                                        #                 parameter = PARAMETER,
                   p.value = PVAL,
                   estimate = ESTIMATE,
                   conf.int = CINT,
                   method = METHOD,
                   data.name = DNAME,
                   observed = x,
                   se0 = se0,
                   se1 = se1),
              class = "htest")

}
