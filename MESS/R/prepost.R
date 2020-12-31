#' Pretest-posttest RCT for quantitative observations with possible missing values
#'
#' In a typical pretest-posttest RCT, subjects are randomized to two treatments, and response is measured at baseline,
#' prior to intervention with the randomized treatment (pretest), and at prespecified follow-up time (posttest).
#' Interest focuses on the effect of treatments on the change between mean baseline and follow-up response.
#' Missing posttest response for some subjects is routine, and disregarding missing cases can lead to invalid inference.
#'
#' @param baseline A vector of quantitative baseline measurements
#' @param post A vector of quantitative post-test measurements with same length as baseline. May contain missing values
#' @param treatment A vector of 0s and 1s corresponding to treatment indicator. 1 = treated, Same length as baseline
#' @param conf.level confidence level of the interval
#' @param delta A numeric between 0 and 1 OR the string "estimate" (the default). The proportion of observation treated.
#' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk}
#' @seealso \code{\link{chisq.test}}
#' @references Marie Davidian, Anastasios A. Tsiatis and Selene Leon (2005).
#' "Semiparametric Estimation of Treatment Effect in a Pretest-Posttest Study
#' with Missing Data". Statistical Science 20, 261-301.
#' @keywords htest
#'
#' @examples
#' # From Altman
#' expo = c(rep(1,9),rep(0,7))
#' bp1w = c(137,120,141,137,140,144,134,123,142,139,134,136,151,147,137,149)
#' bp_base = c(147,129,158,164,134,155,151,141,153,133,129,152,161,154,141,156)
#' diff = bp1w-bp_base
#' prepost.test(bp_base, bp1w, expo)
#'
#' @export
prepost.test <- function(baseline, post, treatment, conf.level = 0.95, delta="estimate") {

    ## Define Z to bypass R CMD check grievances with with below
    Z <- NULL

    ## Check factor
    if ("factor" %in% class(treatment)) {
        treat <- !(treatment==levels(treatment)[1])
    }
    else {
        treat <- !(treatment==0)
    }

    if (length(unique(treatment))!=2)
        stop("Can only handle exactly two treatments")

    if (!(class(treatment) %in% c("integer", "numeric", "logical"))) {
        stop("treatment must be numeric (or logical)")
    }

    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
                                     conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")


    ## Check missing

    ## Check delta

    ## Handle missing baseline/treatment

    ## Convert to notation
    DF <- data.frame(Z  = treatment,    # 0/1, 1 treated
                     R  = !is.na(post), # 0/1 1 if post observed
                     Y1 = baseline,
                     Y2 = post)

    ## Remove missing baselines
    DF <- DF[!is.na(DF$Y1),]

    pi1 <- glm(R ~ Y1, family="binomial", data=DF, subset=(Z==1))
    pi0 <- glm(R ~ Y1, family="binomial", data=DF, subset=(Z==0))

    pihat1 <- predict(pi1, newdata=DF, type="response")
    pihat0 <- predict(pi0, newdata=DF, type="response")

    reg1 <- lm(Y2 ~ Y1, data=DF, subset=(Z==1))
    reg0 <- lm(Y2 ~ Y1, data=DF, subset=(Z==0))

    ehat1 <- predict(reg1, newdata=DF)
    ehat0 <- predict(reg0, newdata=DF)

    N <- nrow(DF)
    N1 <- sum(DF$Z)
    N0 <- N-N1

    deltahat <- N1/N

    ## Use fixed delta if provided
    if (is.numeric(delta) && delta>0 && delta <1) {
        deltahat <- delta
    }

    DF$Y2[is.na(DF$Y2)] <- 0
    mu21 <- with(DF, sum(R*Z*Y2/pihat1 - (Z-deltahat)*ehat1 -(R-pihat1)*Z*ehat1/pihat1 )/N1)
    mu20 <- ifelse(N0==0, 0, with(DF, sum(R*(1-Z)*Y2/pihat0 + (Z-deltahat)*ehat0 - (R-pihat0)*(1-Z)*ehat1/pihat0)/N0))


    mu21 <- with(DF, sum(R*Z*Y2/pihat1 - (Z-deltahat)*ehat1 -(R-pihat1)*Z*ehat1/pihat1 )/ (deltahat*N))
    mu20 <- ifelse(N0==0, 0, with(DF, sum(R*(1-Z)*Y2/pihat0 + (Z-deltahat)*ehat0 - (R-pihat0)*(1-Z)*ehat1/pihat0)/((1-deltahat)*N)))


    betahat <- mu21 - mu20

    iFun1 <- with(DF, R*Z*(Y2-mu21)/(deltahat*pihat1) -
                      (Z-deltahat)*(ehat1 - mu21)/(deltahat) -
                      (R-pihat1)*Z*(ehat1 - mu21)/(deltahat*pihat1))
    iFun2 <- with(DF, -(R*(1-Z)*(Y2-mu20))/((1-deltahat)*pihat0) -
                      ((Z-deltahat)*(ehat0-mu20))/(1-deltahat) -
                      (R-pihat0)*(1-Z)*(ehat1-mu20)/((1-deltahat)*pihat0))

    iFun <- iFun1 + iFun2

    sebetahat <- sqrt(sum((iFun/N)^2))


    method <- "Semiparametric Estimation of Treatment Effect in a Pretest-Posttest Study with Missing Data"

    tstat <- betahat/sebetahat
    alpha <- 1 - conf.level
    cint <- betahat + sebetahat*qnorm(1 - alpha/2) * c(-1, 1)

    names(tstat) <- "z"
    names(betahat) <- c("estimated treatment effect")
    attr(cint, "conf.level") <- conf.level

#    names(mu) <- if (paired || !is.null(y))
#                     "difference in means"
#                 else "mean"


#    rval <- list(statistic = tstat, parameter = df, p.value = 1-pchisq(tstat^2, df=1),
#                 conf.int = cint, estimate = estimate, null.value = mu,
#                 alternative = alternative, method = method, data.name = dname)

    rval <- list(statistic = tstat, se=sebetahat, p.value = 1-pchisq(tstat^2, df=1),
                 conf.int = cint, estimate=betahat, method = method)


    class(rval) <- "htest"
    return(rval)

}

