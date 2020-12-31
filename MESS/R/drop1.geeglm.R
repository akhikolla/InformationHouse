#' Drop All Possible Single Terms to a geeglm Model Using Wald or Score Test
#'
#' Compute all the single terms in the scope argument that can dropped from the
#' model, and compute a table of the corresponding Wald test statistics.
#'
#'
#' @param object a fitted object of class geese.
#' @param scope a formula giving the terms to be considered for adding or dropping.
#' @param test  the type of test to include.
#' @param method Indicates which method is used for computing the standard
#' error. \code{robust} is the default and corresponds to the modified sandwich
#' estimator. \code{naive} is the classical naive cariance estimate.
#' \code{sandwich} is an alias for \code{robust}.
#' @param \dots other arguments. Not currently used
#' @return An object of class "anova" summarizing the differences in fit
#' between the models.
#' @author Claus Ekstrom \email{claus@@ekstroem.dk}
#' @seealso \code{\link{drop1}}, \code{geeglm}, \code{geese}
#' @keywords htest
#' @examples
#'
#' library(geepack)
#' data(ohio)
#' fit <- geeglm(resp ~ age + smoke + age:smoke, id=id, data=ohio,
#'              family=binomial, corstr="exch", scale.fix=TRUE)
#' drop1(fit)
#'
#' @import geepack
#' @export
drop1.geeglm <- function(object, scope, test = c("Wald", "none", "score", "sasscore"),
                         method=c("robust", "naive", "sandwich"), ...) {
  test <- match.arg(test)
  method <- match.arg(method)

  if (! any(c("geeglm", "geem") %in% class(object)) ) {
    stop("Presently drop1.geeglm only works for geeglm or geem objects")
  }

  x <- model.matrix(object)
#  y <- object$y
#  id <- object$id
#  n <- nrow(x)
  asgn <- attr(x, "assign")

  tl <- attr(object$terms, "term.labels")
  if (missing(scope)) {
    scope <- drop.scope(object)
  } else {
    if (!is.character(scope))
      scope <- attr(terms(update.formula(object, scope)),
                    "term.labels")
    if (!all(match(scope, tl, 0L) > 0L))
      stop("scope is not a subset of term labels")
  }

  ndrop <- match(scope, tl)
  ns <- length(scope)

  score <- numeric(ns)
  dfs <- numeric(ns)
  pvals <- numeric(ns)

  for (i in 1L:ns) {
    ii <- seq_along(asgn)[asgn == ndrop[i]]
    jj <- setdiff(seq(ncol(x)), ii)

    if (test %in% c("score", "sasscore")) {
      ans <- update(object, paste(". ~  . -",  scope[i], collapse=""))
      newbeta <- coef(object)
      newbeta[ii] <- 0
      newbeta[jj] <- coef(ans)
      sasscoretest <- switch(test, score = FALSE, sasscore=TRUE)
      score[i] <- scorefct(object, newbeta, testidx=ii, sas=sasscoretest)
    } else {
      param <- coef(object)[ii]
      varmat <- switch(method,
                       robust = object$geese$vbeta[ii,ii],
                       sandwich = object$geese$vbeta[ii,ii],
                       naive = object$geese$vbeta.naiv[ii,ii]
      )
      score[i] <- as.numeric(param %*% solve(varmat) %*% param)
    }
    dfs[i] <- length(ii)
    pvals[i] <- pchisq(score[i], df=dfs[i], lower.tail=FALSE)
  }

  testname <- switch(test, Wald="Wald", score = "Score", sasscore = "Score(SAS)")

  # Start buiding the table
  aod <- data.frame(DF=dfs, row.names=scope)
  aod[, testname] <- score
  if (test!="none")
    aod[,c("Pr(>Chi)")] <- list(pvals)

  head <- c("Single term deletions", "\nModel:", deparse(formula(object)))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}


#' Drop All Possible Single Terms to a geem Model Using Wald or Score Test
#'
#' Compute all the single terms in the scope argument that can dropped from the
#' model, and compute a table of the corresponding Wald test statistics.
#'
#'
#' @param object a fitted object of class geese.
#' @param scope a formula giving the terms to be considered for adding or dropping.
#' @param test  the type of test to include.
#' @param method Indicates which method is used for computing the standard
#' error. \code{robust} is the default and corresponds to the modified sandwich
#' estimator. \code{naive} is the classical naive cariance estimate.
#' \code{sandwich} is an alias for \code{robust}.
#' @param \dots other arguments. Not currently used
#' @return An object of class "anova" summarizing the differences in fit
#' between the models.
#' @author Claus Ekstrom \email{claus@@ekstroem.dk}
#' @seealso \code{\link{drop1}}, \code{geem}
#' @keywords htest
#' @examples
#'
#' library(geeM)
#' library(geepack)
#' data(ohio)
#' \dontrun{
#' fit <- geem(resp ~ age + smoke + age:smoke, id=id, data=ohio,
#'             family="binomial", corstr="exch", scale.fix=TRUE)
#' drop1(fit)
#' }
#'
#' @import geeM
#' @export
drop1.geem <- function(object, scope, test = c("Wald", "none", "score", "sasscore"),
                       method=c("robust", "naive", "sandwich"), ...) {
    test <- match.arg(test)
    method <- match.arg(method)

    if (! any(c("geem") %in% class(object)) ) {
        stop("Presently drop1.geem only works for geeglm or geem objects")
    }

    x <- model.matrix(object)
                                        #  y <- object$y
                                        #  id <- object$id
                                        #  n <- nrow(x)
    asgn <- attr(x, "assign")

    tl <- attr(object$terms, "term.labels")
    if (missing(scope)) {
        scope <- drop.scope(object)
    } else {
        if (!is.character(scope))
            scope <- attr(terms(update.formula(object, scope)),
                          "term.labels")
        if (!all(match(scope, tl, 0L) > 0L))
            stop("scope is not a subset of term labels")
    }

    ndrop <- match(scope, tl)
    ns <- length(scope)

    score <- numeric(ns)
    dfs <- numeric(ns)
    pvals <- numeric(ns)

    for (i in 1L:ns) {
        ii <- seq_along(asgn)[asgn == ndrop[i]]
        jj <- setdiff(seq(ncol(x)), ii)

        if (test %in% c("score", "sasscore")) {
            ans <- update(object, paste(". ~  . -",  scope[i], collapse=""))
            newbeta <- coef(object)
            newbeta[ii] <- 0
            newbeta[jj] <- coef(ans)
            sasscoretest <- switch(test, score = FALSE, sasscore=TRUE)
            score[i] <- scorefct(object, newbeta, testidx=ii, sas=sasscoretest)
        } else {
            param <- coef(object)[ii]
            varmat <- switch(method,
                             robust = object$var[ii,ii],
                             sandwich = object$var[ii,ii],
                             naive = object$naiv.var[ii,ii]
                             )
            score[i] <- as.numeric(param %*% solve(varmat) %*% param)
        }
        dfs[i] <- length(ii)
        pvals[i] <- pchisq(score[i], df=dfs[i], lower.tail=FALSE)
    }

    testname <- switch(test, Wald="Wald", score = "Score", sasscore = "Score(SAS)")

                                        # Start buiding the table
    aod <- data.frame(DF=dfs, row.names=scope)
    aod[, testname] <- score
    if (test!="none")
        aod[,c("Pr(>Chi)")] <- list(pvals)

    head <- c("Single term deletions", "\nModel:", deparse(formula(object)))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}


drop1.geese <- function(object, scope, test = c("Wald", "none"),
                        method=c("robust", "naive", "sandwich"), ...) {
  test <- match.arg(test)
  method <- match.arg(method)

  if (! ("geese" %in% class(object)) ) {
    stop("Requires a geese object")
  }

  glmcall <- object$call
  glmcall$id <- glmcall$jack <- glmcall$control <- glmcall$corstr <- glmcall$waves <- glmcall$zcor <- glmcall$std.err <- glmcall$scale.fix <- glmcall$scale.value<- glmcall$z <- glmcall$family <- NULL
  glmcall[[1]] <- as.name("model.frame")
  mf <- eval(glmcall, parent.frame())
  x <- model.matrix(formula(object), data=mf)
#  n <- nrow(x)
  asgn <- attr(x, "assign")

  tl <- attr(terms(formula(object)), "term.labels")

  if (missing(scope)) {
    scope <- drop.scope(object)
  } else {
    if (!is.character(scope))
      scope <- attr(terms(update.formula(formula(object), scope)),
                    "term.labels")
    if (!all(match(scope, tl, 0L) > 0L))
      stop("scope is not a subset of term labels")
  }

  ndrop <- match(scope, tl)
  ns <- length(scope)

  score <- numeric(ns)
  dfs <- numeric(ns)
  pvals <- numeric(ns)

  for (i in 1L:ns) {
    ii <- seq_along(asgn)[asgn == ndrop[i]]
    param <- coef(object)[ii]
    varmat <- switch(method,
                     robust = object$geese$vbeta[ii,ii],
                     sandwich = object$geese$vbeta[ii,ii],
                     naive = object$geese$vbeta.naiv[ii,ii]
    )
    score[i] <- as.numeric(param %*% solve(varmat) %*% param)
    dfs[i] <- length(ii)
    pvals[i] <- pchisq(score[i], df=dfs[i], lower.tail=FALSE)
  }

  aod <- data.frame(DF=dfs, Wald=score, row.names=scope)
  if (test=="Wald")
    aod[,c("Pr(>Chi)")] <- list(pvals)

  head <- c("Single term deletions", "\nModel:", deparse(formula(object)))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

