#' Factor Analysis for data (high or low dimensional).
#' @rdname fad
#' @description Perform fast matrix-free maximum-likelihood factor analysis on a
#'   covariance matrix or data matrix, works if number of variables is more than
#'   number of observations.
#' @param x A formula or a numeric matrix or an object that can be coerced to a
#'   numeric matrix.
#' @param factors The number of factors to be fitted.
#' @param data An optional data frame (or similar: see \code{\link[stats]{model.frame}}),
#'   used only if \code{x} is a formula.  By default the variables are taken
#'   from \code{environment(formula)}.
#' @param covmat A covariance matrix, or a covariance list as returned by
#'   \code{\link[stats]{cov.wt}}. Of course, correlation matrices are covariance
#'   matrices.
#' @param n.obs The number of observations, used if \code{covmat} is a
#'   covariance matrix.
#' @param subset A specification of the cases to be used, if \code{x} is used as
#'   a matrix or formula.
#' @param na.action The \code{na.action} to be used if \code{x} is used as a
#'   formula.
#' @param start \code{NULL} or a matrix of starting values, each column giving
#'   an initial set of uniquenesses.
#' @param scores Type of scores to produce, if any.  The default is none,
#'   \code{"regression"} gives Thompson's scores, \code{"Bartlett"} given
#'   Bartlett's weighted least-squares scores. Partial matching allows these
#'   names to be abbreviated. Also note that some of the scores-types are not
#'   applicable when \code{p > n}.
#' @param rotation character. \code{"none"} or the name of a function to be used
#'   to rotate the factors: it will be called with first argument the loadings
#'   matrix, and should return a list with component \code{loadings} giving the
#'   rotated loadings, or just the rotated loadings. The options included in the package are:
#'   \code{varimax}, \code{promax}, \code{quartimax}, \code{equamax}.
#' @param control A list of control values:
#' \describe{
#' \item{nstart}{The number of starting values to be tried if \code{start = NULL}. Default 1.}
#' \item{trace}{logical. Output tracing information? Default \code{FALSE}.}
#' \item{opt}{A list of control values to be passed to \code{\link[stats]{optim}}'s
#' \code{control} argument.}
#' \item{rotate}{a list of additional arguments for the rotation function.}
#'}
#' @param lower The lower bound for uniquenesses during optimization. Should be > 0. Default 0.005.
#' @param \dots Components of \code{control} can also be supplied as named arguments to \code{fad}.
#'
#' @return An object of class \code{"fad"} with components
#'\item{loadings}{A matrix of loadings, one column for each factor.  The
#' factors are ordered in decreasing order of sums of squares of
#' loadings, and given the sign that will make the sum of the loadings
#' positive.  This is of class \code{"loadings"}}
#' \item{uniquenesses}{The uniquenesses computed.}
#' \item{criteria}{The results of the optimization: the value of the
#' criterion (a linear function of the negative log-likelihood) and information
#' on the iterations used.}
#' \item{factors}{The argument \code{factors}.}
#' \item{dof}{The number of degrees of freedom of the factor analysis model.}
#' \item{method}{The method: always \code{"mle"}.}
#' \item{rotmat}{The rotation matrix if relevant.}
#' \item{scores}{If requested, a matrix of scores.  \code{napredict} is
#' applied to handle the treatment of values omitted by the \code{na.action}.}
#' \item{n.obs}{The number of observations if available, or \code{NA}.}
#' \item{call}{The matched call.}
#' \item{na.action}{If relevant.}
#' \item{loglik, BIC}{The maximum log-likelihood and the Bayesian Information Criteria.}
#'
#'
#' @seealso \code{\link[stats]{factanal}}
#'  
#'@examples
#'set.seed(1234)
#'
#'## Simulate a 200 x 3 loadings matrix ~i.i.d N(0,1)
#'L <- matrix(rnorm(200*3),200,3) 
#'
#'## Simulate the uniquenesses i.i.d U(0.2,0.9)
#'D <- runif(200,0.2,0.9) 
#'
#'## Generate a data matrix of size 50 x 200 with rows
#'## ~i.i.d. N(0,LL'+diag(D))
#'X <- tcrossprod(matrix(rnorm(50*3),50,3),L) + matrix(rnorm(50*200),50,200) %*% diag(sqrt(D))
#'
#'## Fit a factor model with 3 factors:
#'fit = fad(X,3)
#'
#'
#'## Print the loadings:
#' print(fit$loadings)
#'
#' @export
fad <- function (x, factors, data = NULL, covmat = NULL, n.obs = NA,
                  subset, na.action, start = NULL,
                  scores = c("none", "regression", "Bartlett"),
                  rotation = "varimax",
                  control = NULL, lower=0.005,...)
{

  cl <- match.call()
  na.act <- NULL
  if (is.list(covmat)) {
    if (any(is.na(match(c("cov", "n.obs"), names(covmat)))))
      stop("'covmat' is not a valid covariance list")
    isds <- 1/sqrt(diag(covmat$cov))
    cv <- cov2cor(covmat$cov)
    p <- ncol(cv)
    n.obs <- covmat$n.obs
    have.x <- FALSE
  }
  else if (is.matrix(covmat) || is.Matrix(covmat)) {
    isds <- 1/sqrt(diag(covmat))
    cv <- cov2cor(covmat)
    p <- ncol(cv)
    have.x <- FALSE
  }
  else if (is.null(covmat)) {
    if(missing(x)) stop("neither 'x' nor 'covmat' supplied")
    have.x <- TRUE
    if(inherits(x, "formula")) {
      ## this is not a `standard' model-fitting function,
      ## so no need to consider contrasts or levels
      mt <- terms(x, data = data)
      if(attr(mt, "response") > 0)
        stop("response not allowed in formula")
      attr(mt, "intercept") <- 0
      mf <- match.call(expand.dots = FALSE)
      names(mf)[names(mf) == "x"] <- "formula"
      mf$factors <- mf$covmat <- mf$scores <- mf$start <-
        mf$rotation <- mf$control <- mf$... <- NULL
      ## need stats:: for non-standard evaluation
      mf[[1L]] <- quote(stats::model.frame)
      mf <- eval.parent(mf)
      na.act <- attr(mf, "na.action")
      if (.check_vars_numeric(mf))
        stop("factor analysis applies only to numerical variables")
      z <- model.matrix(mt, mf)
      n.obs <- nrow(z)
      p <- ncol(z)
      means <- colMeans(z)
      if(is.Matrix(z))
      {
        isds <- 1/colSD(z,means)
      } else      isds <- 1/apply(z,2L,sd)*sqrt(n.obs/{n.obs-1})



    } else {
      z <- x
      p <- ncol(z)
      n.obs <- nrow(z)

      if(!is.numeric(z) && !is.Matrix(z))
        stop("factor analysis applies only to numerical variables")
      if(!missing(subset)) z <- z[subset, , drop = FALSE]

      means <- colMeans(z)
      if(is.Matrix(z))
      {
        isds <- 1/colSD(z,means)
      } else      isds <- 1/apply(z,2L,sd)*sqrt(n.obs/{n.obs-1})
    }
  }
  else stop("'covmat' is of unknown type")

  scores <- match.arg(scores)
  if(scores != "none" && !have.x)
    stop("requested scores without an 'x' matrix")



  if(p < 3) stop("factor analysis requires at least three variables")
  dof <- 0.5 * ((p - factors)^2 - p - factors)
  if(dof < 0)
    stop(sprintf(ngettext(factors,
                          "%d factor is too many for %d variables",
                          "%d factors are too many for %d variables"),
                 factors, p), domain = NA)
  
  if(factors == 0){
   if(have.x){
     p <- ncol(z)
     n.obs <- nrow(z)
     mu <- colMeans(z)

     if(is.Matrix(z))
      {
        SD <- colSD(z,mu)/sqrt(n.obs)
      } else{
        SD <- apply(sweep(z,2,mu), 2,function(v) sqrt(mean(v^2)))/sqrt(n.obs)
      }
      
    } else{
      SD <- 1/isds
    }
                    
    loglik = - 0.5*n.obs*p
    BIC = -2*loglik + factors*p*log(n.obs)
                    
    ans <- list(loadings = matrix(0,p,0), uniquenesses = rep(1,p),
                gerr = 0,
                sd = SD,
                loglik = loglik,
                BIC = BIC,
                factors = q, method = "mle")
    class(ans) <- "fad"
    return(ans);
  }


  cn <- list(nstart = 1, trace = FALSE, lower = 0.005)
  cn[names(control)] <- control
  more <- list(...)[c("nstart", "trace", "lower", "opt", "rotate")]
  if(length(more)) cn[names(more)] <- more

  if(have.x) isdsn <- isds/sqrt(n.obs);

  if(is.null(start)) {


    if(have.x)
    {
      start <- cbind(get.start.pc(z,isdsn,means,factors,lower = lower))
    } else
    {
      eigres <- eigs_sym(cv,factors)
      start <- 1 - rowSums( .postmdiag(eigres$vectors,sqrt(eigres$values))^2);
      start[start <= 0] = lower;
      start = cbind(start)
    }
    if((ns <- cn$nstart) > 1)
      start <- cbind(start, matrix(runif(ns-1), p, ns-1, byrow=TRUE))
  }

  start <- as.matrix(start)
  if(nrow(start) != p)
    stop(sprintf(ngettext(p,
                          "'start' must have %d row",
                          "'start' must have %d rows"),
                 p), domain = NA)
  nc <- ncol(start)
  if(nc < 1) stop("no starting values supplied")
  best <- Inf
  for (i in 1L:nc) {
    if(have.x)
    {
      nfit <- fad.fit.X(X = z,q = factors,iSD = isdsn,mu = means,
                         start=start[,i],lower = max(cn$lower,0),control = cn$opt)
    } else
    {
      nfit <- fad.fit.cor(R = cv,q = factors,start=start[,i],
                           lower = max(cn$lower,0),control = cn$opt)
    }

    if(cn$trace)
      cat("start", i, "value:", format(nfit$criteria[1L]),
          "uniqs:", format(as.vector(round(nfit$uniquenesses, 4))), "\n")
    if(nfit$criteria[1L] < best) {
      if(!nfit$converged && nfit$gerr > p*1e-5) warning("Algorithm may not have converged. Try another starting value.")
      fit <- nfit
      best <- fit$criteria[1L]
    }
  }
  if(best == Inf)
    stop(ngettext(nc,
                  "unable to optimize from this starting value",
                  "unable to optimize from these starting values"),
         domain = NA)
  load <- fit$loadings
  if(rotation != "none") {
    rot <- do.call(rotation, c(list(load), cn$rotate))
    load <- if (is.list(rot)) {
      load <- rot$loadings
      fit$rotmat <-
        if(inherits(rot, "GPArotation")) t(solve(rot$Th))
      else rot$rotmat
      rot$loadings
    } else rot
  }
  fit$loadings <- sortLoadings(load)
  class(fit$loadings) <- "loadings"
  fit$na.action <- na.act # not used currently
  if(have.x && scores != "none") {
    Lambda <- fit$loadings
    switch(scores,
           regression = {
             if(p >= n.obs) {
               sprintf("Cannot compute regression score when p >= n");
               sc = NULL
             } else {
               zz <- scale(z, TRUE, TRUE)
               sc <- zz %*% solve(cv, Lambda)
               if(!is.null(Phi <- attr(Lambda, "covariance")))
                 sc <- sc %*% Phi
             }
           },
           Bartlett = {
             d <- 1/fit$uniquenesses
             tmp <- Lambda * d
             tmp <- isds*{tmp %*% solve(crossprod(tmp, Lambda))}
             sc <- z %*% tmp - sum(means*tmp)
           })
    rownames(sc) <- rownames(z)
    colnames(sc) <- colnames(Lambda)
    if(!is.null(na.act)) sc <- napredict(na.act, sc)
    fit$scores <- sc
  }
  # if(!is.na(n.obs) && dof > 0) {
  #   fit$STATISTIC <- (n.obs - 1 - (2 * p + 5)/6 -
  #                       (2 * factors)/3) * fit$criteria["objective"]
  #   fit$PVAL <- pchisq(fit$STATISTIC, dof, lower.tail = FALSE)
  # }


  Gamma = diag(factors) + crossprod(fit$loadings,{1/fit$uniquenesses}*fit$loadings);
  R = chol(Gamma)
  logdet = sum(log(fit$uniquenesses)) + 2*sum(log(diag(R)))
  logdet0 = 2*sum(log(isds));
  fit$loglik = {-logdet + logdet0}*0.5*n.obs - n.obs*p/2;
  if(is.na(n.obs)) warning("Number of obs. not supplied. BIC cannot be computed")
  fit$BIC = -2*fit$loglik + factors*p*log(n.obs);
  fit$sd = 1/isds;
  fit$n.obs <- n.obs
  fit$call <- cl
  fit
}


sortLoadings <- function(Lambda)
{
  cn <- colnames(Lambda)
  Phi <- attr(Lambda, "covariance")
  ssq <- -colSums(Lambda^2) #apply(Lambda, 2L, function(x) -sum(x^2))
  Lambda <- Lambda[, order(ssq), drop = FALSE]
  colnames(Lambda) <- cn
  neg <- colSums(Lambda) < 0
  Lambda[, neg] <- -Lambda[, neg]
  if(!is.null(Phi)) {
    unit <- ifelse(neg, -1, 1)
    attr(Lambda, "covariance") <-
      unit %*% Phi[order(ssq), order(ssq)] %*% unit
  }
  Lambda
}

.check_vars_numeric <- function(mf) {
  mt <- attr(mf, "terms")
  mterms <- attr(mt, "factors")
  mterms <- rownames(mterms)[apply(mterms, 1L, function(x) any(x > 0L))]
  any(sapply(mterms, function(x) is.factor(mf[, x]) || !is.numeric(mf[, x])))
}


