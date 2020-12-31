
# Package: survival
# Version: 2.44-1.1
# Authors@R: c(person(c("Terry", "M"), "Therneau",
#                     email="therneau.terry@mayo.edu",
#                     role=c("aut", "cre")),
#              person("Thomas", "Lumley", role=c("ctb", "trl"),
#                     comment="original S->R port and R maintainer until 2009"))
# License: LGPL (>= 2)

# coxph function from survival package version 2.44-1.1
# followed by non-exported functions of the same package 
# (potentially used by coxph)

coxph<-function (formula, data, weights, subset, na.action, init, control, 
          ties = c("efron", "breslow", "exact"), singular.ok = TRUE, 
          robust = FALSE, model = FALSE, x = FALSE, y = TRUE, tt, method = ties, 
          ...) 
{
  ties <- match.arg(ties)
  Call <- match.call()
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(coxph.control))
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L)
    if (any(indx == 0L)) 
      stop(gettextf("Argument %s not matched", names(extraArgs)[indx == 
                                                                  0L]), domain = NA)
  }
  if (missing(control)) 
    control <- coxph.control(...)
  indx <- match(c("formula", "data", "weights", "subset", "na.action"), 
                names(Call), nomatch = 0)
  if (indx[1] == 0) 
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  special <- c("strata", "cluster", "tt")
  temp$formula <- if (missing(data)) 
    terms(formula, special)
  else terms(formula, special, data = data)
  if (!is.null(attr(temp$formula, "specials")$tt)) {
    coxenv <- new.env(parent = environment(formula))
    assign("tt", function(x) x, env = coxenv)
    environment(temp$formula) <- coxenv
  }
  mf <- eval(temp, parent.frame())
  if (nrow(mf) == 0) 
    stop("No (non-missing) observations")
  Terms <- terms(mf)
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")
  type <- attr(Y, "type")
  if (type != "right" && type != "counting") 
    stop(paste("Cox model doesn't support \"", type, "\" survival data", 
               sep = ""))
  data.n <- nrow(Y)
  if (control$timefix) 
    Y <- aeqSurv(Y)
  if (length(attr(Terms, "variables")) > 2) {
    ytemp <- terms.inner(formula[1:2])
    xtemp <- terms.inner(formula[-2])
    if (any(!is.na(match(xtemp, ytemp)))) 
      warning("a variable appears on both the left and right sides of the formula")
  }
  strats <- attr(Terms, "specials")$strata
  if (length(strats)) {
    stemp <- untangle.specials(Terms, "strata", 1)
    if (length(stemp$vars) == 1) 
      strata.keep <- mf[[stemp$vars]]
    else strata.keep <- strata(mf[, stemp$vars], shortlabel = TRUE)
    strats <- as.numeric(strata.keep)
  }
  timetrans <- attr(Terms, "specials")$tt
  if (missing(tt)) 
    tt <- NULL
  if (length(timetrans)) {
    timetrans <- untangle.specials(Terms, "tt")
    ntrans <- length(timetrans$terms)
    if (is.null(tt)) {
      tt <- function(x, time, riskset, weights) {
        obrien <- function(x) {
          r <- rank(x)
          (r - 0.5)/(0.5 + length(r) - r)
        }
        unlist(tapply(x, riskset, obrien))
      }
    }
    if (is.function(tt)) 
      tt <- list(tt)
    if (is.list(tt)) {
      if (any(!sapply(tt, is.function))) 
        stop("The tt argument must contain function or list of functions")
      if (length(tt) != ntrans) {
        if (length(tt) == 1) {
          temp <- vector("list", ntrans)
          for (i in 1:ntrans) temp[[i]] <- tt[[1]]
          tt <- temp
        }
        else stop("Wrong length for tt argument")
      }
    }
    else stop("The tt argument must contain a function or list of functions")
    if (ncol(Y) == 2) {
      if (length(strats) == 0) {
        sorted <- order(-Y[, 1], Y[, 2])
        newstrat <- rep.int(0L, nrow(Y))
        newstrat[1] <- 1L
      }
      else {
        sorted <- order(strats, -Y[, 1], Y[, 2])
        newstrat <- as.integer(c(1, 1 * (diff(strats[sorted]) != 
                                           0)))
      }
      if (storage.mode(Y) != "double") 
        storage.mode(Y) <- "double"
      counts <- .Call(Ccoxcount1, Y[sorted, ], as.integer(newstrat))
      tindex <- sorted[counts$index]
    }
    else {
      if (length(strats) == 0) {
        sort.end <- order(-Y[, 2], Y[, 3])
        sort.start <- order(-Y[, 1])
        newstrat <- c(1L, rep(0, nrow(Y) - 1))
      }
      else {
        sort.end <- order(strats, -Y[, 2], Y[, 3])
        sort.start <- order(strats, -Y[, 1])
        newstrat <- c(1L, as.integer(diff(strats[sort.end]) != 
                                       0))
      }
      if (storage.mode(Y) != "double") 
        storage.mode(Y) <- "double"
      counts <- .Call(Ccoxcount2, Y, as.integer(sort.start - 
                                                  1L), as.integer(sort.end - 1L), as.integer(newstrat))
      tindex <- counts$index
    }
    Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
    type <- "right"
    mf <- mf[tindex, ]
    strats <- rep(1:length(counts$nrisk), counts$nrisk)
    weights <- model.weights(mf)
    if (!is.null(weights) && any(!is.finite(weights))) 
      stop("weights must be finite")
    tcall <- attr(Terms, "variables")[timetrans$terms + 2]
    pvars <- attr(Terms, "predvars")
    pmethod <- sub("makepredictcall.", "", as.vector(methods("makepredictcall")))
    for (i in 1:ntrans) {
      newtt <- (tt[[i]])(mf[[timetrans$var[i]]], Y[, 1], 
                         strats, weights)
      mf[[timetrans$var[i]]] <- newtt
      nclass <- class(newtt)
      if (any(nclass %in% pmethod)) {
        dummy <- as.call(list(as.name(class(newtt)[1]), 
                              tcall[[i]][[2]]))
        ptemp <- makepredictcall(newtt, dummy)
        pvars[[timetrans$terms[i] + 2]] <- ptemp
      }
    }
    attr(Terms, "predvars") <- pvars
  }
  cluster <- attr(Terms, "specials")$cluster
  if (length(cluster)) {
    robust <- TRUE
    tempc <- untangle.specials(Terms, "cluster", 1:10)
    ord <- attr(Terms, "order")[tempc$terms]
    if (any(ord > 1)) 
      stop("Cluster can not be used in an interaction")
    cluster <- strata(mf[, tempc$vars], shortlabel = TRUE)
    dropterms <- tempc$terms
    xlevels <- .getXlevels(Terms[-tempc$terms], mf)
  }
  else {
    if (missing(robust)) 
      robust <- FALSE
    xlevels <- .getXlevels(Terms, mf)
    dropterms <- NULL
  }
  contrast.arg <- NULL
  attr(Terms, "intercept") <- 1
  stemp <- untangle.specials(Terms, "strata", 1)
  hasinteractions <- FALSE
  if (length(stemp$vars) > 0) {
    for (i in stemp$vars) {
      if (any(attr(Terms, "order")[attr(Terms, "factors")[i, 
      ] > 0] > 1)) 
        hasinteractions <- TRUE
    }
    if (!hasinteractions) 
      dropterms <- c(dropterms, stemp$terms)
  }
  if (length(dropterms)) {
    Terms2 <- Terms[-dropterms]
    X <- model.matrix(Terms2, mf, constrasts = contrast.arg)
    temp <- attr(X, "assign")
    shift <- sort(dropterms)
    temp <- temp + 1 * (shift[1] <= temp)
    if (length(shift) == 2) 
      temp + 1 * (shift[2] <= temp)
    attr(X, "assign") <- temp
  }
  else X <- model.matrix(Terms, mf, contrasts = contrast.arg)
  if (!all(is.finite(X))) 
    stop("data contains an infinite predictor")
  Xatt <- attributes(X)
  if (hasinteractions) 
    adrop <- c(0, untangle.specials(Terms, "strata")$terms)
  else adrop <- 0
  xdrop <- Xatt$assign %in% adrop
  X <- X[, !xdrop, drop = FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  attr(X, "contrasts") <- Xatt$contrasts
  offset <- model.offset(mf)
  if (is.null(offset) | all(offset == 0)) 
    offset <- rep(0, nrow(mf))
  else if (any(!is.finite(offset) | !is.finite(exp(offset)))) 
    stop("offsets must lead to a finite risk score")
  weights <- model.weights(mf)
  if (!is.null(weights) && any(!is.finite(weights))) 
    stop("weights must be finite")
  assign <- attrassign(X, Terms)
  contr.save <- attr(X, "contrasts")
  if (missing(init)) 
    init <- NULL
  else {
    if (length(init) != ncol(X)) 
      stop("wrong length for init argument")
    temp <- X %*% init - sum(colMeans(X) * init) + offset
    if (any(exp(temp) > .Machine$double.xmax) || all(exp(temp) == 
                                                     0)) 
      stop("initial values lead to overflow or underflow of the exp function")
  }
  if (sum(Y[, ncol(Y)]) == 0) {
    ncoef <- ncol(X)
    ctemp <- rep(NA, ncoef)
    names(ctemp) <- colnames(X)
    concordance = c(concordant = 0, discordant = 0, tied.x = 0, 
                    tied.y = 0, tied.xy = 0, concordance = NA, std = NA, 
                    timefix = FALSE)
    rval <- list(coefficients = ctemp, var = matrix(0, ncoef, 
                                                    ncoef), loglik = c(0, 0), score = 0, iter = 0, linear.predictors = offset, 
                 residuals = rep(0, data.n), means = colMeans(X), 
                 method = method, n = data.n, nevent = 0, terms = Terms, 
                 assign = assign, concordance = concordance, y = Y, 
                 call = Call)
    class(rval) <- "coxph"
    return(rval)
  }
  pterms <- sapply(mf, inherits, "coxph.penalty")
  if (any(pterms)) {
    pattr <- lapply(mf[pterms], attributes)
    pname <- names(pterms)[pterms]
    ord <- attr(Terms, "order")[match(pname, attr(Terms, 
                                                  "term.labels"))]
    if (any(ord > 1)) 
      stop("Penalty terms cannot be in an interaction")
    pcols <- assign[match(pname, names(assign))]
    fit <- coxpenal.fit(X, Y, strats, offset, init = init, 
                        control, weights = weights, method = method, row.names(mf), 
                        pcols, pattr, assign)
  }
  else {
    if (method == "breslow" || method == "efron") {
      if (type == "right") 
        fitter <- get("coxph.fit")
      else fitter <- get("agreg.fit")
    }
    else if (method == "exact") {
      if (type == "right") 
        fitter <- get("coxexact.fit")
      else fitter <- get("agexact.fit")
    }
    else stop(paste("Unknown method", method))
    fit <- fitter(X, Y, strats, offset, init, control, weights = weights, 
                  method = method, row.names(mf))
  }
  if (is.character(fit)) {
    fit <- list(fail = fit)
    class(fit) <- "coxph"
  }
  else {
    if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
      vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
      msg <- paste("X matrix deemed to be singular; variable", 
                   paste(vars, collapse = " "))
      if (!singular.ok) 
        stop(msg)
    }
    fit$n <- data.n
    fit$nevent <- sum(Y[, ncol(Y)])
    fit$terms <- Terms
    fit$assign <- assign
    class(fit) <- fit$class
    fit$class <- NULL
    if (robust && !is.null(fit$coefficients) && !all(is.na(fit$coefficients))) {
      fit$naive.var <- fit$var
      fit2 <- c(fit, list(x = X, y = Y, weights = weights))
      if (length(strats)) 
        fit2$strata <- strats
      if (length(cluster)) {
        temp <- residuals.coxph(fit2, type = "dfbeta", 
                                collapse = cluster, weighted = TRUE)
        if (is.null(init)) 
          fit2$linear.predictors <- 0 * fit$linear.predictors
        else fit2$linear.predictors <- c(X %*% init)
        temp0 <- residuals.coxph(fit2, type = "score", 
                                 collapse = cluster, weighted = TRUE)
      }
      else {
        temp <- residuals.coxph(fit2, type = "dfbeta", 
                                weighted = TRUE)
        fit2$linear.predictors <- 0 * fit$linear.predictors
        temp0 <- residuals.coxph(fit2, type = "score", 
                                 weighted = TRUE)
      }
      fit$var <- t(temp) %*% temp
      u <- apply(as.matrix(temp0), 2, sum)
      fit$rscore <- coxph.wtest(t(temp0) %*% temp0, u, 
                                control$toler.chol)$test
    }
    if (length(fit$coefficients) && is.null(fit$wald.test)) {
      nabeta <- !is.na(fit$coefficients)
      if (is.null(init)) 
        temp <- fit$coefficients[nabeta]
      else temp <- (fit$coefficients - init[1:length(fit$coefficients)])[nabeta]
      fit$wald.test <- coxph.wtest(fit$var[nabeta, nabeta], 
                                   temp, control$toler.chol)$test
    }
    if (length(cluster)) 
      temp <- concordancefit(Y, fit$linear.predictors, 
                             strats, weights, cluster = cluster, reverse = TRUE, 
                             timefix = FALSE)
    else temp <- concordancefit(Y, fit$linear.predictors, 
                                strats, weights, reverse = TRUE, timefix = FALSE)
    if (is.matrix(temp$count)) 
      fit$concordance <- c(colSums(temp$count), concordance = temp$concordance, 
                           std = sqrt(temp$var))
    else fit$concordance <- c(temp$count, concordance = temp$concordance, 
                              std = sqrt(temp$var))
    na.action <- attr(mf, "na.action")
    if (length(na.action)) 
      fit$na.action <- na.action
    if (model) {
      if (length(timetrans)) {
        mf[[".surv."]] <- Y
        mf[[".strata."]] <- strats
        stop("'model=TRUE' not supported for models with tt terms")
      }
      fit$model <- mf
    }
    if (x) {
      fit$x <- X
      if (length(strats)) {
        if (length(timetrans)) 
          fit$strata <- strats
        else fit$strata <- strata.keep
      }
    }
    if (y) 
      fit$y <- Y
  }
  if (!is.null(weights) && any(weights != 1)) 
    fit$weights <- weights
  names(fit$means) <- names(fit$coefficients)
  fit$formula <- formula(Terms)
  if (length(xlevels) > 0) 
    fit$xlevels <- xlevels
  fit$contrasts <- contr.save
  if (any(offset != 0)) 
    fit$offset <- offset
  fit$call <- Call
  fit
}


terms.inner<-function (x) 
{
  if (inherits(x, "formula")) {
    if (length(x) == 3) 
      c(terms.inner(x[[2]]), terms.inner(x[[3]]))
    else terms.inner(x[[2]])
  }
  else if (inherits(x, "call") && (x[[1]] != as.name("$") && 
                                   x[[1]] != as.name("["))) {
    if (x[[1]] == "+" || x[[1]] == "*" || x[[1]] == "-") {
      if (length(x) == 3) 
        c(terms.inner(x[[2]]), terms.inner(x[[3]]))
      else terms.inner(x[[2]])
    }
    else if (x[[1]] == as.name("Surv")) 
      unlist(lapply(x[-1], terms.inner))
    else terms.inner(x[[2]])
  }
  else (deparse(x))
}

coxpenal.df<-function (hmat, hinv, fdiag, assign.list, ptype, nvar, pen1, 
                       pen2, sparse) 
{
  if (ptype == 1 & nvar == 0) {
    hdiag <- 1/fdiag
    list(fvar2 = (hdiag - pen1) * fdiag^2, df = sum((hdiag - 
                                                       pen1) * fdiag), fvar = fdiag, trH = sum(fdiag))
  }
  else if (ptype == 2) {
    hmat.full <- t(hmat) %*% (ifelse(fdiag == 0, 0, 1/fdiag) * 
                                hmat)
    hinv.full <- hinv %*% (fdiag * t(hinv))
    if (length(pen2) == length(hmat.full)) 
      imat <- hmat.full - pen2
    else imat <- hmat.full - diag(pen2)
    var <- hinv.full %*% imat %*% hinv.full
    if (length(assign.list) == 1) 
      list(var2 = var, df = sum(imat * hinv.full), trH = sum(diag(hinv.full)), 
           var = hinv.full)
    else {
      df <- trH <- NULL
      d2 <- diag(hinv.full)
      for (i in assign.list) {
        temp <- coxph.wtest(hinv.full[i, i], var[i, i])$solve
        if (is.matrix(temp)) 
          df <- c(df, sum(diag(temp)))
        else df <- c(df, sum(temp))
        trH <- c(trH, sum(d2[i]))
      }
      list(var2 = var, df = df, trH = trH, var = hinv.full)
    }
  }
  else {
    nf <- length(fdiag) - nvar
    nr1 <- 1:nf
    nr2 <- (nf + 1):(nf + nvar)
    d1 <- fdiag[nr1]
    d2 <- fdiag[nr2]
    temp <- t(hinv[nr1, ])
    temp2 <- t(hinv[nr2, , drop = FALSE])
    A.diag <- d1 + c(rep(1, nvar) %*% (temp^2 * d2))
    B <- hinv[nr1, ] %*% (d2 * temp2)
    C <- hinv[nr2, ] %*% (d2 * temp2)
    var2 <- C - t(B) %*% (pen1 * B)
    if (ptype == 3) {
      hmat.22 <- t(hmat) %*% (ifelse(fdiag == 0, 0, 1/fdiag) * 
                                hmat)
      temp <- C - coxph.wtest(hmat.22, diag(nvar))$solve
      if (nvar == 1) {
        var2 <- var2 - C * pen2 * C
        temp2 <- c(temp * pen2)
      }
      else if (length(pen2) == nvar) {
        var2 <- var2 - C %*% (pen2 * C)
        temp2 <- sum(diag(temp) * pen2)
      }
      else {
        var2 <- var2 - C %*% matrix(pen2, nvar) %*% C
        temp2 <- sum(diag(temp * pen2))
      }
    }
    else temp2 <- 0
    df <- trH <- NULL
    cdiag <- diag(C)
    for (i in 1:length(assign.list)) {
      if (sparse == i) {
        df <- c(df, nf - (sum(A.diag * pen1) + temp2))
        trH <- c(trH, sum(A.diag))
      }
      else {
        j <- assign.list[[i]]
        temp <- coxph.wtest(C[j, j], var2[j, j])$solve
        if (is.matrix(temp)) 
          df <- c(df, sum(diag(temp)))
        else df <- c(df, sum(temp))
        trH <- c(trH, sum(cdiag[j]))
      }
    }
    list(var = C, df = df, trH = trH, fvar = A.diag, var2 = var2)
  }
}

coxpenal.fit<-function (x, y, strata, offset, init, control, weights, method, 
                        rownames, pcols, pattr, assign) 
{
  eps <- control$eps
  n <- nrow(y)
  if (is.matrix(x)) 
    nvar <- ncol(x)
  else if (length(x) == 0) 
    stop("Must have an X variable")
  else nvar <- 1
  if (missing(offset) || is.null(offset)) 
    offset <- rep(0, n)
  if (missing(weights) || is.null(weights)) 
    weights <- rep(1, n)
  else {
    if (any(weights <= 0)) 
      stop("Invalid weights, must be >0")
  }
  if (ncol(y) == 3) {
    if (length(strata) == 0) {
      sorted <- cbind(order(-y[, 2], y[, 3]), order(-y[, 
                                                       1])) - 1L
      newstrat <- as.integer(n)
    }
    else {
      sorted <- cbind(order(strata, -y[, 2], y[, 3]), order(strata, 
                                                            -y[, 1])) - 1L
      newstrat <- as.integer(cumsum(table(strata)))
    }
    status <- y[, 3]
    andersen <- TRUE
  }
  else {
    if (length(strata) == 0) {
      sorted <- order(-y[, 1], y[, 2]) - 1L
      newstrat <- as.integer(n)
    }
    else {
      sorted <- order(strata, -y[, 1], y[, 2]) - 1L
      newstrat <- as.integer(cumsum(table(strata)))
    }
    status <- y[, 2]
    andersen <- FALSE
  }
  n.eff <- sum(y[, ncol(y)])
  npenal <- length(pattr)
  if (npenal == 0 || length(pcols) != npenal) 
    stop("Invalid pcols or pattr arg")
  sparse <- sapply(pattr, function(x) !is.null(x$sparse) && 
                     x$sparse)
  if (sum(sparse) > 1) 
    stop("Only one sparse penalty term allowed")
  pterms <- rep(0, length(assign))
  names(pterms) <- names(assign)
  pindex <- rep(0, npenal)
  for (i in 1:npenal) {
    temp <- unlist(lapply(assign, function(x, y) (length(x) == 
                                                    length(y) && all(x == y)), pcols[[i]]))
    if (sparse[i]) 
      pterms[temp] <- 2
    else pterms[temp] <- 1
    pindex[i] <- (seq(along.with = temp))[temp]
  }
  if ((sum(pterms == 2) != sum(sparse)) || (sum(pterms > 0) != 
                                            npenal)) 
    stop("pcols and assign arguments disagree")
  if (any(pindex != sort(pindex))) {
    temp <- order(pindex)
    pindex <- pindex[temp]
    pcols <- pcols[temp]
    pattr <- pattr[temp]
  }
  ptype <- any(sparse) + 2 * (any(!sparse))
  f.expr1 <- function(coef) NULL
  f.expr2 <- function(coef) NULL
  if (any(sparse)) {
    sparse.attr <- (pattr[sparse])[[1]]
    fcol <- unlist(pcols[sparse])
    if (length(fcol) > 1) 
      stop("Sparse term must be single column")
    xx <- x[, -fcol, drop = FALSE]
    for (i in 1:length(assign)) {
      j <- assign[[i]]
      if (j[1] > fcol) 
        assign[[i]] <- j - 1
    }
    for (i in 1:npenal) {
      j <- pcols[[i]]
      if (j[1] > fcol) 
        pcols[[i]] <- j - 1
    }
    frailx <- x[, fcol]
    frailx <- match(frailx, sort(unique(frailx)))
    nfrail <- max(frailx)
    nvar <- nvar - 1
    pfun1 <- sparse.attr$pfun
    f.expr1 <- function(coef) {
      coxlist1$coef <- coef
      if (is.null(extra1)) 
        temp <- pfun1(coef, theta1, n.eff)
      else temp <- pfun1(coef, theta1, n.eff, extra1)
      if (!is.null(temp$recenter)) 
        coxlist1$coef <- coxlist1$coef - as.double(temp$recenter)
      if (!temp$flag) {
        coxlist1$first <- -as.double(temp$first)
        coxlist1$second <- as.double(temp$second)
      }
      coxlist1$penalty <- -as.double(temp$penalty)
      coxlist1$flag <- as.logical(temp$flag)
      if (any(sapply(coxlist1, length) != c(rep(nfrail, 
                                                3), 1, 1))) 
        stop("Incorrect length in coxlist1")
      coxlist1
    }
    if (!is.null(getOption("survdebug"))) 
      debug(f.expr1)
    coxlist1 <- list(coef = double(nfrail), first = double(nfrail), 
                     second = double(nfrail), penalty = 0, flag = FALSE)
  }
  else {
    xx <- x
    frailx <- 0
    nfrail <- 0
  }
  if (sum(!sparse) > 0) {
    full.imat <- !all(unlist(lapply(pattr, function(x) x$diag)))
    ipenal <- (1:length(pattr))[!sparse]
    f.expr2 <- function(coef) {
      coxlist2$coef <- coef
      pentot <- 0
      for (i in ipenal) {
        pen.col <- pcols[[i]]
        coef <- coxlist2$coef[pen.col]
        if (is.null(extralist[[i]])) 
          temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]], 
                                      n.eff)
        else temp <- ((pattr[[i]])$pfun)(coef, thetalist[[i]], 
                                         n.eff, extralist[[i]])
        if (!is.null(temp$recenter)) 
          coxlist2$coef[pen.col] <- coxlist2$coef[pen.col] - 
          temp$recenter
        if (temp$flag) 
          coxlist2$flag[pen.col] <- TRUE
        else {
          coxlist2$flag[pen.col] <- FALSE
          coxlist2$first[pen.col] <- -temp$first
          if (full.imat) {
            tmat <- matrix(coxlist2$second, nvar, nvar)
            tmat[pen.col, pen.col] <- temp$second
            coxlist2$second <- c(tmat)
          }
          else coxlist2$second[pen.col] <- temp$second
        }
        pentot <- pentot - temp$penalty
      }
      coxlist2$penalty <- as.double(pentot)
      if (any(sapply(coxlist2, length) != length2)) 
        stop("Length error in coxlist2")
      coxlist2
    }
    if (!is.null(getOption("survdebug"))) 
      debug(f.expr2)
    if (full.imat) {
      coxlist2 <- list(coef = double(nvar), first = double(nvar), 
                       second = double(nvar * nvar), penalty = 0, flag = rep(FALSE, 
                                                                             nvar))
      length2 <- c(nvar, nvar, nvar * nvar, 1, nvar)
    }
    else {
      coxlist2 <- list(coef = double(nvar), first = double(nvar), 
                       second = double(nvar), penalty = 0, flag = rep(FALSE, 
                                                                      nvar))
      length2 <- c(nvar, nvar, nvar, 1, nvar)
    }
  }
  else full.imat <- FALSE
  if (nfrail > 0) 
    finit <- rep(0, nfrail)
  else finit <- 0
  if (!missing(init) && !is.null(init)) {
    if (length(init) != nvar) {
      if (length(init) == (nvar + nfrail)) {
        finit <- init[-(1:nvar)]
        init <- init[1:nvar]
      }
      else stop("Wrong length for inital values")
    }
  }
  else init <- double(nvar)
  cfun <- lapply(pattr, function(x) x$cfun)
  parmlist <- lapply(pattr, function(x, eps) c(x$cparm, eps2 = eps), 
                     sqrt(eps))
  extralist <- lapply(pattr, function(x) x$pparm)
  iterlist <- vector("list", length(cfun))
  thetalist <- vector("list", length(cfun))
  printfun <- lapply(pattr, function(x) x$printfun)
  for (i in 1:length(cfun)) {
    temp <- (cfun[[i]])(parmlist[[i]], iter = 0)
    if (sparse[i]) {
      theta1 <- temp$theta
      extra1 <- extralist[[i]]
    }
    thetalist[[i]] <- temp$theta
    iterlist[[i]] <- temp
  }
  temp1 <- c("x", "coef", "plik", "loglik", "status", "neff", 
             "df", "trH")
  temp2 <- c("frailx", "coxfit$fcoef", "loglik1", "coxfit$loglik", 
             "status", "n.eff")
  temp3 <- c("xx[,pen.col]", "coxfit$coef[pen.col]", "loglik1", 
             "coxfit$loglik", "status", "n.eff")
  calls <- vector("expression", length(cfun))
  cargs <- lapply(pattr, function(x) x$cargs)
  for (i in 1:length(cfun)) {
    tempchar <- paste("(cfun[[", i, "]])(parmlist[[", i, 
                      "]], iter,", "iterlist[[", i, "]]")
    temp2b <- c(temp2, paste("pdf[", i, "]"), paste("trH[", 
                                                    i, "]"))
    temp3b <- c(temp3, paste("pdf[", i, "]"), paste("trH[", 
                                                    i, "]"))
    if (length(cargs[[i]]) == 0) 
      calls[i] <- parse(text = paste(tempchar, ")"))
    else {
      temp <- match(cargs[[i]], temp1)
      if (any(is.na(temp))) 
        stop(paste((cargs[[i]])[is.na(temp)], "not matched"))
      if (sparse[i]) 
        temp4 <- paste(temp2b[temp], collapse = ",")
      else temp4 <- paste(temp3b[temp], collapse = ",")
      calls[i] <- parse(text = paste(paste(tempchar, temp4, 
                                           sep = ","), ")"))
    }
  }
  need.df <- any(!is.na(match(c("df", "trH"), unlist(cargs))))
  varnames <- dimnames(xx)[[2]]
  for (i in 1:length(cfun)) {
    if (!is.null(pattr[[i]]$varname)) 
      varnames[pcols[[i]]] <- pattr[[i]]$varname
  }
  rho <- environment()
  storage.mode(y) <- storage.mode(weights) <- "double"
  storage.mode(xx) <- storage.mode(offset) <- "double"
  if (andersen) 
    coxfit <- .C(Cagfit5a, as.integer(n), as.integer(nvar), 
                 y, xx, offset, weights, newstrat, sorted, means = double(nvar), 
                 coef = as.double(init), u = double(nvar), loglik = double(1), 
                 as.integer(method == "efron"), as.integer(ptype), 
                 as.integer(full.imat), as.integer(nfrail), as.integer(frailx), 
                 f.expr1, f.expr2, rho)
  else coxfit <- .C(Ccoxfit5a, as.integer(n), as.integer(nvar), 
                    y, xx, offset, weights, newstrat, sorted, means = double(nvar), 
                    coef = as.double(init), u = double(nvar), loglik = double(1), 
                    as.integer(method == "efron"), as.integer(ptype), as.integer(full.imat), 
                    as.integer(nfrail), as.integer(frailx), f.expr1, f.expr2, 
                    rho)
  loglik0 <- coxfit$loglik
  means <- coxfit$means
  iter2 <- 0
  iterfail <- NULL
  thetasave <- unlist(thetalist)
  for (outer in 1:control$outer.max) {
    if (andersen) 
      coxfit <- .C(Cagfit5b, iter = as.integer(control$iter.max), 
                   as.integer(n), as.integer(nvar), as.integer(newstrat), 
                   coef = as.double(init), u = double(nvar + nfrail), 
                   hmat = double(nvar * (nvar + nfrail)), hinv = double(nvar * 
                                                                          (nvar + nfrail)), loglik = double(1), flag = integer(1), 
                   as.double(control$eps), as.double(control$toler.chol), 
                   as.integer(method == "efron"), as.integer(nfrail), 
                   fcoef = as.double(finit), fdiag = double(nfrail + 
                                                              nvar), f.expr1, f.expr2, rho)
    else coxfit <- .C(Ccoxfit5b, iter = as.integer(control$iter.max), 
                      as.integer(n), as.integer(nvar), as.integer(newstrat), 
                      coef = as.double(init), u = double(nvar + nfrail), 
                      hmat = double(nvar * (nvar + nfrail)), hinv = double(nvar * 
                                                                             (nvar + nfrail)), loglik = double(1), flag = integer(1), 
                      as.double(control$eps), as.double(control$toler.chol), 
                      as.integer(method == "efron"), as.integer(nfrail), 
                      fcoef = as.double(finit), fdiag = double(nfrail + 
                                                                 nvar), f.expr1, f.expr2, rho)
    iter <- outer
    iter2 <- iter2 + coxfit$iter
    if (coxfit$iter >= control$iter.max) 
      iterfail <- c(iterfail, iter)
    temp <- rep(FALSE, nvar + nfrail)
    if (nfrail > 0) 
      temp[1:nfrail] <- coxlist1$flag
    if (ptype > 1) 
      temp[nfrail + 1:nvar] <- coxlist2$flag
    fdiag <- ifelse(temp, 0, coxfit$fdiag)
    if (need.df) {
      if (nfrail > 0) 
        temp1 <- coxlist1$second
      else temp1 <- 0
      if (ptype > 1) 
        temp2 <- coxlist2$second
      else temp2 <- 0
      dftemp <- coxpenal.df(matrix(coxfit$hmat, ncol = nvar), 
                            matrix(coxfit$hinv, ncol = nvar), fdiag, assign, 
                            ptype, nvar, temp1, temp2, pindex[sparse])
      df <- dftemp$df
      var <- dftemp$var
      var2 <- dftemp$var2
      pdf <- df[pterms > 0]
      trH <- dftemp$trH[pterms > 0]
    }
    if (nfrail > 0) 
      penalty <- -coxlist1$penalty
    else penalty <- 0
    if (ptype > 1) 
      penalty <- penalty - coxlist2$penalty
    loglik1 <- coxfit$loglik + penalty
    if (iter == 1) 
      penalty0 <- penalty
    done <- TRUE
    for (i in 1:length(cfun)) {
      pen.col <- pcols[[i]]
      temp <- eval(calls[i])
      if (sparse[i]) 
        theta1 <- temp$theta
      thetalist[[i]] <- temp$theta
      iterlist[[i]] <- temp
      done <- done & temp$done
    }
    if (done) 
      break
    if (iter == 1) {
      init <- coefsave <- coxfit$coef
      finit <- fsave <- coxfit$fcoef
      thetasave <- cbind(thetasave, unlist(thetalist))
    }
    else {
      temp <- as.vector(unlist(thetalist))
      coefsave <- cbind(coefsave, coxfit$coef)
      fsave <- cbind(fsave, coxfit$fcoef)
      howclose <- apply((thetasave - temp)^2, 2, sum)
      which <- min((1:iter)[howclose == min(howclose)])
      if (nvar > 0) 
        init <- coefsave[, which]
      if (nfrail > 0) 
        finit <- fsave[, which]
      thetasave <- cbind(thetasave, temp)
    }
  }
  if (nfrail > 0) {
    lp <- offset + coxfit$fcoef[frailx]
    if (nvar > 0) 
      lp <- lp + x[, -fcol, drop = FALSE] %*% coxfit$coef - 
        sum(means * coxfit$coef)
  }
  else lp <- offset + as.vector(x %*% coxfit$coef) - sum(means * 
                                                           coxfit$coef)
  if (andersen) {
    .C(Cagfit5c, as.integer(nvar))
    resid <- .Call(Cagmart3, y, exp(lp), weights, newstrat, 
                   sorted, as.integer(method == "efron"))
  }
  else {
    expect <- .C(Ccoxfit5c, as.integer(n), as.integer(nvar), 
                 as.integer(newstrat), as.integer(method == "efron"), 
                 expect = double(n))$expect
    resid <- status - expect
  }
  names(resid) <- rownames
  if (!need.df) {
    if (nfrail > 0) 
      temp1 <- coxlist1$second
    else temp1 <- 0
    if (ptype > 1) 
      temp2 <- coxlist2$second
    else temp2 <- 0
    dftemp <- coxpenal.df(matrix(coxfit$hmat, ncol = nvar), 
                          matrix(coxfit$hinv, ncol = nvar), fdiag, assign, 
                          ptype, nvar, temp1, temp2, pindex[sparse])
    df <- dftemp$df
    trH <- dftemp$trH
    var <- dftemp$var
    var2 <- dftemp$var2
  }
  if (control$iter.max > 1 && length(iterfail) > 0) 
    warning(paste("Inner loop failed to coverge for iterations", 
                  paste(iterfail, collapse = " ")))
  which.sing <- (fdiag[nfrail + 1:nvar] == 0)
  coef <- coxfit$coef
  names(coef) <- varnames
  coef[which.sing] <- NA
  names(iterlist) <- names(pterms[pterms > 0])
  if (nfrail > 0) {
    if (nvar > 0) {
      list(coefficients = coef, var = var, var2 = var2, 
           loglik = c(loglik0, loglik1), iter = c(iter, 
                                                  iter2), linear.predictors = as.vector(lp), 
           residuals = resid, means = means, method = method, 
           class = c("coxph.penal", "coxph"), frail = coxfit$fcoef, 
           fvar = dftemp$fvar, df = df, df2 = dftemp$df2, 
           penalty = c(penalty0, penalty), pterms = pterms, 
           assign2 = assign, history = iterlist, coxlist1 = coxlist1, 
           printfun = printfun)
    }
    else {
      list(loglik = c(loglik0, loglik1), iter = c(iter, 
                                                  iter2), linear.predictors = as.vector(lp), residuals = resid, 
           means = means, method = method, class = c("coxph.penal", 
                                                     "coxph"), frail = coxfit$fcoef, fvar = dftemp$fvar, 
           df = df, df2 = dftemp$df2, penalty = c(penalty0, 
                                                  penalty), pterms = pterms, assign2 = assign, 
           history = iterlist, printfun = printfun)
    }
  }
  else {
    list(coefficients = coef, var = var, var2 = var2, loglik = c(loglik0, 
                                                                 loglik1), iter = c(iter, iter2), linear.predictors = lp, 
         residuals = resid, means = means, method = method, 
         class = c("coxph.penal", "coxph"), df = df, df2 = dftemp$df2, 
         penalty = c(penalty0, penalty), pterms = pterms, 
         assign2 = assign, history = iterlist, coxlist2 = coxlist2, 
         printfun = printfun)
  }
}

parsecovar1<-function (flist, statedatanames) 
{
  if (any(sapply(flist, function(x) !inherits(x, "formula")))) 
    stop("flist must be a list of formulas")
  if (any(sapply(flist, length) != 3)) 
    stop("all formulas must have a left and right side")
  lhs <- lapply(flist, function(x) x[-3])
  rhs <- lapply(flist, function(x) x[-2])
  rh2 <- lapply(rhs, function(form) {
    parts <- strsplit(deparse(form, width.cutoff = 300, control = NULL), 
                      "/", fixed = TRUE)[[1]]
    if (length(parts) == 1) {
      ival <- NULL
      common <- FALSE
      fixed <- FALSE
      clear <- FALSE
    }
    else {
      optterms <- terms(formula(paste("~", parts[2])))
      ff <- rownames(attr(optterms, "factors"))
      index <- match(ff, c("common", "fixed", "init", "clear"))
      if (any(is.na(index))) 
        stop("option not recognized in a covariates formula: ", 
             paste(ff[is.na(index)], collapse = ", "))
      common <- any(index == 1)
      fixed <- any(index == 2)
      clear <- any(index == 3)
      if (any(index == 3)) {
        optatt <- attributes(optterms)
        j <- optatt$variables[1 + which(index == 3)]
        j[[1]] <- as.name("list")
        ival <- unlist(eval(j, parent.frame()))
      }
      else ival <- NULL
    }
    form <- formula(paste("~ -1 +", parts[1]))
    list(common = common, fixed = fixed, clear = clear, ival = ival, 
         formula = form)
  })
  pcut <- function(form) {
    if (length(form) == 3) {
      if (form[[1]] == "+") 
        c(pcut(form[[2]]), pcut(form[[3]]))
      else if (form[[1]] == "~") 
        pcut(form[[2]])
      else list(form)
    }
    else list(form)
  }
  lcut <- lapply(lhs, function(x) pcut(x[[2]]))
  env1 <- new.env(parent = parent.frame(2))
  env2 <- new.env(parent = env1)
  if (missing(statedatanames)) {
    assign("state", function(...) list(stateid = "state", 
                                       values = c(...)), env1)
    assign("state", list(stateid = "state"))
  }
  else {
    for (i in statedatanames) {
      assign(i, eval(list(stateid = i)), env2)
      tfun <- eval(parse(text = paste0("function(...) list(stateid='", 
                                       i, "', values=c(...)")))
      assign(i, tfun, env1)
    }
  }
  lterm <- lapply(lcut, function(x) {
    lapply(x, function(z) {
      if (length(z) == 1) {
        temp <- eval(z, envir = env2)
        if (is.list(temp) && names(temp)[[1]] == "stateid") 
          temp
        else temp
      }
      else if (length(z) == 3 && z[[1]] == ":") 
        list(left = eval(z[[2]], envir = env2), right = eval(z[[3]], 
                                                             envir = env2))
      else stop("invalid term: ", deparse(z))
    })
  })
  list(rhs = rh2, lhs = lterm)
}

parsecovar2<-function (covar1, statedata, dformula, Terms, transitions, states) 
{
  if (is.null(statedata)) 
    statedata <- data.frame(state = states)
  else {
    if (is.null(statedata$state)) 
      stop("the statedata data set must contain a variable 'state'")
    indx1 <- match(states, statedata$state, nomatch = 0)
    if (any(indx1 == 0)) 
      stop("statedata does not contain all the possible states: ", 
           states[indx1 == 0])
    statedata <- statedata[indx1, ]
  }
  allterm <- attr(Terms, "term.labels")
  nterm <- length(allterm)
  nstate <- length(states)
  tmap <- array(0L, dim = c(nterm + 1, nstate, nstate))
  dterms <- match(attr(terms.formula(dformula), "term.labels"), 
                  allterm)
  dterms <- c(1L, 1L + dterms)
  k <- seq(along = dterms)
  for (i in 1:nstate) {
    for (j in 1:nstate) {
      tmap[dterms, i, j] <- k
      k <- k + length(k)
    }
  }
  ncoef <- max(tmap)
  inits <- NULL
  if (!is.null(covar1)) {
    for (i in 1:length(covar1$rhs)) {
      rhs <- covar1$rhs[[i]]
      lhs <- covar1$lhs[[i]]
      rterm <- terms.formula(rhs$formula)
      rindex <- 1L + match(attr(rterm, "term.labels"), 
                           allterm, nomatch = 0)
      if (any(rindex == 1L)) 
        stop("dterm mismatch bug 2")
      if (attr(rterm, "intercept") == 1) 
        rindex <- c(1L, rindex)
      state1 <- state2 <- NULL
      for (x in lhs) {
        if (is.null(x$left)) 
          stop("term found without a :", x)
        if (!is.list(x$left) && length(x$left) == 1 & 
            x$left == 1) 
          temp1 <- 1:nrow(statedata)
        else if (is.list(x$left) && names(x$left)[1] == 
                 "stateid") {
          if (is.null(x$left$value)) 
            stop("state variable with no list of values: ", 
                 x$left$stateid)
          else temp1 <- which(statedata[[x$left$stateid]] %in% 
                                x$left$value)
        }
        else temp1 <- which(statedata$state %in% x$left)
        if (!is.list(x$right) && length(x$right) == 1 && 
            x$right == 1) 
          temp2 <- 1:nrow(statedata)
        else if (is.list(x$right) && names(x$right)[1] == 
                 "stateid") {
          if (is.null(x$right$value)) 
            stop("state variable with no list of values: ", 
                 x$right$stateid)
          else temp2 <- which(statedata[[x$right$stateid]] %in% 
                                x$right$value)
        }
        else temp2 <- which(statedata$state %in% x$right)
        state1 <- c(state1, rep(temp1, length(temp2)))
        state2 <- c(state2, rep(temp2, each = length(temp1)))
      }
      if (rhs$clear) 
        tmap[-1, state1, state2] <- 0
      if (length(rhs$ival)) 
        inits <- c(inits, list(term = rindex, state1 = state1, 
                               state2 = state2, init = rhs$ival))
      if (rhs$common) 
        j <- ncoef + seq_len(length(rindex))
      else j <- ncoef + seq_len(length(rindex) * length(state1))
      tmap[rindex, state1, state2] <- j
      ncoef <- max(j)
    }
  }
  t2 <- transitions[, -1, drop = FALSE]
  indx1 <- match(rownames(t2), states)
  indx2 <- match(colnames(t2), states)
  tmap2 <- matrix(0L, nrow = 1 + nterm, ncol = sum(t2 > 0))
  trow <- row(t2)[t2 > 0]
  tcol <- col(t2)[t2 > 0]
  for (i in 1:length(trow)) tmap2[, i] <- tmap[, indx1[trow[i]], 
                                               indx2[tcol[i]]]
  dimnames(tmap2) <- list(c("Intercept", allterm), paste(indx1[trow], 
                                                         indx2[tcol], sep = ":"))
  list(tmap = tmap2, inits = inits, mapid = cbind(indx1[trow], 
                                                  indx2[tcol]))
}

parsecovar3<-function (tmap, Xcol, Xassign) 
{
  hasintercept <- (Xassign[1] == 0)
  cmap <- matrix(0L, length(Xcol) + !hasintercept, ncol(tmap))
  cmap[1, ] <- match(tmap[1, ], sort(c(0, unique(tmap[1, ])))) - 
    1L
  xcount <- table(factor(Xassign, levels = 1:max(Xassign)))
  mult <- 1 + max(xcount)
  j <- 1
  for (i in 2:nrow(tmap)) {
    k <- seq_len(xcount[i - 1])
    cmap[j + k, ] <- ifelse(tmap[i, ] == 0, 0, tmap[i, ] * 
                              mult + rep(k, ncol(tmap)))
    j <- j + max(k)
  }
  cmap[-1, ] <- match(cmap[-1, ], sort(unique(c(0L, cmap[-1, 
  ])))) - 1L
  colnames(cmap) <- colnames(tmap)
  if (hasintercept) 
    rownames(cmap) <- Xcol
  else rownames(cmap) <- c("(Intercept)", Xcol)
  cmap
}

residuals.coxph<-function (object, type = c("martingale", "deviance", "score", 
                           "schoenfeld", "dfbeta", "dfbetas", "scaledsch", "partial"), 
          collapse = FALSE, weighted = FALSE, ...) 
{
  type <- match.arg(type)
  otype <- type
  if (type == "dfbeta" || type == "dfbetas") {
    otype <- type
    type <- "score"
    if (missing(weighted)) 
      weighted <- TRUE
  }
  if (type == "scaledsch") 
    type <- "schoenfeld"
  n <- length(object$residuals)
  rr <- object$residuals
  y <- object$y
  x <- object[["x"]]
  vv <- drop(object$naive.var)
  if (is.null(vv)) 
    vv <- drop(object$var)
  weights <- object$weights
  if (is.null(weights)) 
    weights <- rep(1, n)
  strat <- object$strata
  method <- object$method
  if (method == "exact" && (type == "score" || type == "schoenfeld")) 
    stop(paste(otype, "residuals are not available for the exact method"))
  if (type == "martingale" || type == "partial") 
    rr <- object$residuals
  else {
    Terms <- object$terms
    if (!inherits(Terms, "terms")) 
      stop("invalid terms component of object")
    strats <- attr(Terms, "specials")$strata
    if (is.null(y) || (is.null(x) && type != "deviance")) {
      temp <- coxph.getdata(object, y = TRUE, x = TRUE, 
                            stratax = TRUE)
      y <- temp$y
      x <- temp$x
      if (length(strats)) 
        strat <- temp$strata
    }
    ny <- ncol(y)
    status <- y[, ny, drop = TRUE]
    if (type != "deviance") {
      nstrat <- as.numeric(strat)
      nvar <- ncol(x)
      if (is.null(strat)) {
        ord <- order(y[, ny - 1], -status)
        newstrat <- rep(0, n)
      }
      else {
        ord <- order(nstrat, y[, ny - 1], -status)
        newstrat <- c(diff(as.numeric(nstrat[ord])) != 
                        0, 1)
      }
      newstrat[n] <- 1
      x <- x[ord, ]
      y <- y[ord, ]
      score <- exp(object$linear.predictors)[ord]
    }
  }
  if (type == "schoenfeld") {
    if (ny == 2) {
      mintime <- min(y[, 1])
      if (mintime < 0) 
        y <- cbind(2 * mintime - 1, y)
      else y <- cbind(-1, y)
    }
    temp <- .C(Ccoxscho, n = as.integer(n), as.integer(nvar), 
               as.double(y), resid = as.double(x), as.double(score * 
                                                               weights[ord]), as.integer(newstrat), as.integer(method == 
                                                                                                                 "efron"), double(3 * nvar))
    deaths <- y[, 3] == 1
    if (nvar == 1) 
      rr <- temp$resid[deaths]
    else rr <- matrix(temp$resid[deaths], ncol = nvar)
    if (weighted) 
      rr <- rr * weights[deaths]
    if (length(strats)) 
      attr(rr, "strata") <- table((strat[ord])[deaths])
    time <- c(y[deaths, 2])
    if (is.matrix(rr)) 
      dimnames(rr) <- list(time, names(object$coefficients))
    else names(rr) <- time
    if (otype == "scaledsch") {
      ndead <- sum(deaths)
      coef <- ifelse(is.na(object$coefficients), 0, object$coefficients)
      rr <- drop(rr %*% vv * ndead + rep(coef, each = nrow(rr)))
    }
    return(rr)
  }
  if (type == "score") {
    if (ny == 2) {
      resid <- .C(Ccoxscore, as.integer(n), as.integer(nvar), 
                  as.double(y), x = as.double(x), as.integer(newstrat), 
                  as.double(score), as.double(weights[ord]), as.integer(method == 
                                                                          "efron"), resid = double(n * nvar), double(2 * 
                                                                                                                       nvar))$resid
    }
    else {
      resid <- .C(Cagscore, as.integer(n), as.integer(nvar), 
                  as.double(y), as.double(x), as.integer(newstrat), 
                  as.double(score), as.double(weights[ord]), as.integer(method == 
                                                                          "efron"), resid = double(n * nvar), double(nvar * 
                                                                                                                       6))$resid
    }
    if (nvar > 1) {
      rr <- matrix(0, n, nvar)
      rr[ord, ] <- matrix(resid, ncol = nvar)
      dimnames(rr) <- list(names(object$residuals), names(object$coefficients))
    }
    else rr[ord] <- resid
    if (otype == "dfbeta") {
      if (is.matrix(rr)) 
        rr <- rr %*% vv
      else rr <- rr * vv
    }
    else if (otype == "dfbetas") {
      if (is.matrix(rr)) 
        rr <- (rr %*% vv) %*% diag(sqrt(1/diag(vv)))
      else rr <- rr * sqrt(vv)
    }
  }
  if (weighted) 
    rr <- rr * weights
  if (!is.null(object$na.action)) {
    rr <- naresid(object$na.action, rr)
    if (is.matrix(rr)) 
      n <- nrow(rr)
    else n <- length(rr)
    if (type == "deviance") 
      status <- naresid(object$na.action, status)
  }
  if (type == "partial") {
    rr <- rr + predict(object, type = "terms")
  }
  if (!missing(collapse)) {
    if (length(collapse) != n) 
      stop("Wrong length for 'collapse'")
    rr <- drop(rowsum(rr, collapse))
    if (type == "deviance") 
      status <- drop(rowsum(status, collapse))
  }
  if (type == "deviance") 
    sign(rr) * sqrt(-2 * (rr + ifelse(status == 0, 0, status * 
                                        log(status - rr))))
  else rr
}

coxph.getdata<-function (fit, y = TRUE, x = TRUE, stratax = TRUE, offset = FALSE) 
{
  ty <- fit[["y"]]
  tx <- fit[["x"]]
  strat <- fit$strata
  Terms <- fit$terms
  if (is.null(attr(Terms, "offset"))) 
    offset <- FALSE
  if (offset) 
    x <- TRUE
  if (!inherits(Terms, "terms")) 
    stop("invalid terms component of fit")
  strats <- attr(Terms, "specials")$strata
  if (length(strats) == 0) 
    stratax <- FALSE
  if ((y && is.null(ty)) || (x && is.null(tx)) || (stratax && 
                                                   is.null(strat)) || offset) {
    m <- stats::model.frame(fit)
    if (y && is.null(ty)) 
      ty <- model.extract(m, "response")
    if (offset) 
      toff <- model.extract(m, "offset")
    if ((x || stratax) && is.null(tx)) {
      if (stratax) {
        temp <- untangle.specials(Terms, "strata", 1)
        strat <- strata(m[temp$vars], shortlabel = TRUE)
      }
      if (x) 
        tx <- model.matrix(fit, data = m)
    }
  }
  else if (offset) 
    toff <- fit$linear.predictors - (c(tx %*% fit$coef) - 
                                       sum(fit$means * fit$coef))
  temp <- list()
  if (y) 
    temp$y <- ty
  if (x) 
    temp$x <- tx
  if (stratax) 
    temp$strata <- strat
  if (offset) 
    temp$offset <- toff
  temp
}

# Create a strata variable, possibly from many objects
#
strata <- function(..., na.group=FALSE, shortlabel, sep=', ') {
  # First, grab a copy of the call, which will be used to manufacture
  #  labels for unlabeled arguments
  # Then get the arguments as a list
  words <- as.character((match.call())[-1])
  allf <- list(...)
  # If there is only one argument, and it itself is a list, use
  #  it instead
  if(length(allf) == 1 && is.list(ttt <- unclass(allf[[1]]))) allf <- ttt
  nterms <- length(allf)
  
  # Keep the names of named args as their label, what was typed otherwise
  if (is.null(names(allf))) {
    argname <- words[1:nterms]
    if (missing(shortlabel))
      shortlabel <- all(sapply(allf, 
                               function(x) is.character(x) | inherits(x, 'factor')))
  }
  else {
    argname <- ifelse(names(allf) == '', words[1:nterms], names(allf))
    if (missing(shortlabel)) shortlabel <- FALSE
  }
  
  # If the arguments are not all the same length, stop now    
  # Mostly this is to stop calls with an improper object
  arglength <- sapply(allf, length)
  if (any(arglength != arglength[1])) 
    stop("all arguments must be the same length")
  if (!all(sapply(allf, is.atomic))) stop("all arguments must be vectors")
  
  # Process the first argument
  what <- allf[[1]]
  if(is.null(levels(what)))
    what <- factor(what)
  levs <- unclass(what) - 1
  wlab <- levels(what)
  if (na.group && any(is.na(what))){
    # add "NA" as a level
    levs[is.na(levs)] <- length(wlab)
    wlab <- c(wlab, "NA")
  }
  
  if (shortlabel) labs <- wlab
  else            labs <- paste(argname[1], wlab, sep='=')
  
  # Now march through the other variables, if any
  for (i in (1:nterms)[-1]) {
    what <- allf[[i]]
    if(is.null(levels(what)))
      what <- factor(what)
    wlab <- levels(what)
    wlev <- unclass(what) - 1
    if (na.group && any(is.na(wlev))){
      wlev[is.na(wlev)] <- length(wlab)
      wlab <- c(wlab, "NA")
    }
    if (!shortlabel) wlab <- format(paste(argname[i], wlab, sep='='))
    levs <- wlev + levs*(length(wlab))
    labs <- paste(rep(labs, rep(length(wlab), length(labs))),
                  rep(wlab, length(labs)), sep=sep)
  }
  levs <- levs + 1
  ulevs <- sort(unique(levs[!is.na(levs)]))
  levs <- match(levs, ulevs)
  labs <- labs[ulevs]
  
  factor(levs, labels=labs)
}
