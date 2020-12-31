#' Predicted hazard and survival curves from an \code{emfrail} object
#'
#' @importFrom stats .getXlevels delete.response drop.terms as.formula qnorm
#' @param object An \code{emfrail} fit object
#' @param newdata A data frame with the same variable names as those that appear in the \code{emfrail} formula, used to calculate the \code{lp} (optional).
#' @param lp A vector of linear predictor values at which to calculate the curves. Default is 0 (baseline).
#' @param strata The name of the strata (if applicable) for which the prediction should be made.
#' @param quantity Can be \code{"cumhaz"} and/or \code{"survival"}. The quantity to be calculated for the values of \code{lp}.
#' @param type Can be \code{"conditional"} and/or \code{"marginal"}. The type of the quantity to be calculated.
#' @param conf_int Can be \code{"regular"} and/or \code{"adjusted"}. The type of confidence interval to be calculated.
#' @param conf_level The width of the confidence intervals. By default, 95\% confidence intervals are calculated.
#' @param individual Logical. Are the observations in \code{newdata} from the same individual? See details.
#' @param ... Ignored
#'
#' @return The return value is a single data frame (if \code{lp} has length 1,
#' \code{newdata} has 1 row or \code{individual == TRUE}) or a list of data frames corresponding to each value of
#' \code{lp} or each row of \code{newdata} otherwise.
#' The names of the columns in the returned data frames are as follows: \code{time} represents the unique event time points
#' from the data set, \code{lp} is the value of the linear predictor (as specified in the input or as calculated from the lines of \code{newdata}).
#' By default, for each \code{lp} a data frame will contain the following columns: \code{cumhaz}, \code{survival},
#' \code{cumhaz_m}, \code{survival_m} for the cumulative hazard and survival, conditional and marginal, with corresponding confidence
#' bands. The naming of the columns is explained more in the Details section.
#'
#' @details
#' The function calculates predicted cumulative hazard and survival curves for given covariate
#' or linear predictor values; for the first, \code{newdata} must be specified and for the latter
#' \code{lp} must be specified. Each row of \code{newdata} or element of \code{lp} is consiered to be
#' a different subject, and the desired predictions are produced for each of them separately.
#'
#' In \code{newdata} two columns may be specified with the names \code{tstart} and \code{tstop}.
#' In this case, each subject is assumed to be at risk only during the times specified by these two values.
#' If the two are not specified, the predicted curves are produced for a subject that is at risk for the
#' whole follow-up time.
#'
#' A slightly different behaviour is observed if \code{individual == TRUE}. In this case, all the rows of
#' \code{newdata} are assumed to come from the same individual, and \code{tstart} and \code{tstop} must
#' be specified, and must not overlap. This may be used for describing subjects that
#' are not at risk during certain periods or subjects with time-dependent covariate values.
#'
#' The two "quantities" that can be returned are
#' named \code{cumhaz} and \code{survival}. If we denote each quantity with \code{q}, then the columns with the marginal estimates
#' are named \code{q_m}. The confidence intervals contain the name of the quantity (conditional or marginal) followed by \code{_l} or \code{_r} for
#' the lower and upper bound. The bounds calculated with the adjusted standard errors have the name of the regular bounds followed by
#' \code{_a}. For example, the adjusted lower bound for the marginal survival is in the column named \code{survival_m_l_a}.
#'
#' The \code{emfrail} only gives the Breslow estimates of the  baseline hazard \eqn{\lambda_0(t)} at the
#' event time points, conditional on the frailty. Let \eqn{\lambda(t)} be the baseline hazard for a linear predictor of interest.
#' The estimated conditional cumulative hazard is then
#' \eqn{\Lambda(t) = \sum_{s= 0}^t \lambda(s)}. The variance of \eqn{\Lambda(t)} can be calculated from the (maybe adjusted)
#' variance-covariance matrix.
#'
#' The conditional survival is obtained by the usual expression \eqn{S(t) = \exp(-\Lambda(t))}. The marginal survival
#' is given by
#' \deqn{\bar S(t) = E \left[\exp(-\Lambda(t)) \right] = \mathcal{L}(\Lambda(t)),}
#' i.e. the Laplace transform of the frailty distribution calculated in \eqn{\Lambda(t)}.
#'
#' The marginal hazard is obtained as \deqn{\bar \Lambda(t) = - \log \bar S(t).}
#'
#' The only standard errors that are available from \code{emfrail} are those for \eqn{\lambda_0(t)}. From this,
#' standard errors of \eqn{\log \Lambda(t)} may be calculated. On this scale, the symmetric confidence intervals are built, and then
#' moved to the desired scale.
#'
#' @note The linear predictor is taken as fixed, so the variability in the estimation of the regression coefficient is not taken into account.
#' Does not support left truncation (at the moment). That is because, if \code{individual == TRUE} and \code{tstart} and \code{tstop} are
#' specified, for the marginal estimates the distribution of the frailty is used to calculate the integral, and not
#' the distribution of the frailty given the truncation.
#'
#' For performance reasons, consider running with \code{conf_int = NULL}; the reason is that the \code{deltamethod} function that is used
#' to calculate the confidence intervals easily becomes slow when there is a large number of time points
#' for the cumulative hazard.
#'
#' @export
#'
#'
#' @examples
#' kidney$sex <- ifelse(kidney$sex == 1, "male", "female")
#' m1 <- emfrail(formula = Surv(time, status) ~  sex + age  + cluster(id),
#'               data =  kidney)
#'
#' # get all the possible prediction for the value 0 of the linear predictor
#' predict(m1, lp = 0)
#'
#' # get the cumulative hazards for two different values of the linear predictor
#' predict(m1, lp = c(0, 1), quantity = "cumhaz", conf_int = NULL)
#'
#' # get the cumulative hazards for a female and for a male, both aged 30
#' newdata1 <- data.frame(sex = c("female", "male"),
#'                        age = c(30, 30))
#'
#' predict(m1, newdata = newdata1, quantity = "cumhaz", conf_int = NULL)
#'
#' # get the cumulative hazards for an individual that changes
#' # sex from female to male at time 40.
#' newdata2 <- data.frame(sex = c("female", "male"),
#'                       age = c(30, 30),
#'                       tstart = c(0, 40),
#'                       tstop = c(40, Inf))
#'
#' predict(m1, newdata = newdata2,
#'         individual = TRUE,
#'         quantity = "cumhaz", conf_int = NULL)
#'
#' @seealso \code{\link{plot.emfrail}}, \code{\link{autoplot.emfrail}}
predict.emfrail <- function(object,
                            newdata = NULL,
                            lp = NULL,
                            strata = NULL,
                            quantity = c("cumhaz", "survival"),
                            type = c("conditional", "marginal"),
                            conf_int = NULL,
                            individual = FALSE,
                            conf_level = 0.95,
                            ...) {


  if(!is.null(quantity))
    if(any(!(quantity %in% c("cumhaz", "survival"))))
      stop("quantity misspecified")

  if(!is.null(type))
    if(any(!(type %in% c("conditional", "marginal"))))
      stop("type misspecified")

  # if(!is.null(conf_int))
  #   if(any(!(conf_int %in% c("regular", "adjusted"))))
  #     warning("conf_int misspecified")

  if(is.null(newdata) & is.null(lp)) stop("lp or newdata must be specified")



  if(is.list(object$hazard))
    if(length(object$hazard) > 1) {
      if(is.null(strata)) stop("strata must be specified")
      if(!(strata %in% names(object$hazard)))
        stop(paste0("strata must be one of: {", paste0(names(object$hazard), collapse = ", "), "}"))
      if(length(strata) > 1) stop("predict() only available for one strata")
  }



  if(!is.null(lp) & !is.null(newdata)) stop("specify either lp or newdata")

  if(!is.null(lp)) {
    if(isTRUE(individual)) stop("if lp is specified, individual must be FALSE")
    if(!is.numeric(lp)) stop("lp must be a numeric vector")
    tstart <- 0
    tstop <- Inf
  }

  if(!is.null(newdata)) {
    if(!inherits(newdata, "data.frame")) stop("newdata must be a data.frame")

    # browser()
    # here the point is ot take newdata and make it into a value of linear predictor
    # this is the hacky way of getting some stuff from the model fit
    mdata <- attr(object, "metadata")
    # [[1]] is terms
    mf <- model.frame(mdata[[1]], data = newdata, xlev = mdata[[2]])
    mm <- try(model.matrix(mdata[[1]], mf)[,-1,drop=FALSE])
    if(inherits(mm, "try-error")) stop("newdata probably misspecified")

    if(!is.null(strata))
      mm <- mm[1:nrow(mm), -grep("strata", dimnames(mm)[[2]])]

    # browser()
    lp <- as.numeric(mm %*% object$coef)

    if(!("tstart" %in% names(newdata))) {
      tstart = 0
      tstop = Inf
    } else
        if(!("tstop" %in% names(newdata)))
          stop("if tstart is specified then also tstop must be specified") else {
            if(any(newdata$tstop <= newdata$tstart)) stop("tstop must be larger than tstart")

            # check if there is an overlap
            if(individual == TRUE) {
              if(any(diff(as.numeric(as.matrix(newdata[c("tstart", "tstop")]))) < 0))
                stop("(tstart,tstop) must be ordered and without overlap")
            }
            tstart <- newdata$tstart
            tstop <- newdata$tstop
          }

    # at this point I should have a vector of lp
    # if individual == TRUE then I should also have tstart tstop for each lp
  }

  # there is a hack here so that the cumulative hazard actually starts from 0.
  # for this I add a fake time point before the start of the predicted curve.


  # names(object$tev)


  # select the correct strata
  if(!is.null(strata)) {
    tev <- object$tev[[strata]]
    hazard <- object$hazard[[strata]]
  } else {
    tev <- object$tev
    hazard <- object$hazard
  }


  delta_time <- min(diff(tev)) /2

  list_haz <- mapply(FUN = function(lp, haz, tstart, tstop, tev) {
            cbind(tev[tstart <= tev & tev < tstop],
                  lp,
                  haz[tstart <= tev & tev < tstop] * exp(lp))
    }, lp, list(hazard), tstart, tstop, list(tev), SIMPLIFY = FALSE)

  if(isTRUE(individual)) {
      res <- as.data.frame(do.call(rbind, list_haz))
      res <- rbind(c(min(res[,1] - delta_time), lp[1], 0), res)
      names(res) <- c("time", "lp", "cumhaz")
      res$cumhaz <- cumsum(res$cumhaz) # a bit dodgy
      attr(res, "bounds") <- "cumhaz"
      res <- list(res)
    } else
      res <- lapply(list_haz, function(x) {
                    res <- rbind(c(min(x[,1]- delta_time), x[1,2], 0), as.data.frame(x))
                    names(res) <- c("time", "lp", "cumhaz")
                    res$cumhaz <- cumsum(res$cumhaz)
                    attr(res, "bounds") <- "cumhaz"
              res})


  ncoef <- length(object$coef)


  # Here we start putting a bunch of standard errors and stuff inside
  # for each lp we will get a data frame with a "bounds" attribute.
  # this attributie is meant to keep track of which columns we have in the data frame

  # select the correct varH


  if(!is.null(strata)) {
    pos <- which(rep(names(object$hazard), sapply(object$hazard, length))== strata) + ncoef
  } else
    pos <- ncoef + 1:length(hazard)


  if(("regular" %in% conf_int) | ("adjusted" %in% conf_int)) {


    varH <-   object$var[pos, pos]
    varH_adj <- object$var_adj[pos, pos]

    x <- res[[1]]
    res <- lapply(res, function(x) {
      times_res <- match(x$time[2:length(x$time)], tev)

      loghaz <- log(hazard[times_res])

      # Now here I build a bunch of formulas - that is to get the confidence intervals with the delta mehtod
      xs <- lapply(seq_along(times_res), function(x) text1 <- paste0("x", x))
      for(i in 2:length(xs)) {
        xs[[i]] = paste0(xs[[i-1]], " + ", xs[[i]])
      }
      forms <- lapply(xs, function(x) as.formula(paste0("~log(", x, ")")))

      # attr(x, "bounds") <- "cumhaz"

      if("regular" %in% conf_int) {

        x$se_logH <- c(0, msm::deltamethod(g = forms,
                                      mean =hazard[times_res],
                                      cov = varH[times_res, times_res],
                                      ses = TRUE))

        attr(x, "bounds") <- c(attr(x, "bounds"), "cumhaz_l", "cumhaz_r")
        x$cumhaz_l <- pmax(0, exp(log(x$cumhaz) + qnorm((1 - conf_level) / 2) * x$se_logH))
        x$cumhaz_r <- exp(log(x$cumhaz) - qnorm((1 - conf_level) / 2) *x$se_logH)
      }

      if("adjusted" %in% conf_int) {
        x$se_logH_adj <- c(0, msm::deltamethod(g = forms,
                                          mean = hazard[times_res],
                                          cov = varH_adj[times_res, times_res],
                                          ses = TRUE))

        attr(x, "bounds") <- c(attr(x, "bounds"), "cumhaz_l_a", "cumhaz_r_a")
        x$cumhaz_l_a <- pmax(0, exp(log(x$cumhaz) + qnorm((1 - conf_level) / 2) *x$se_logH_adj))
        x$cumhaz_r_a <- exp(log(x$cumhaz) - qnorm((1 - conf_level) / 2)*x$se_logH_adj)
      }



      x
    })
  }




  est_dist <- object$distribution
  est_dist$frailtypar <- exp(object$logtheta)

  # the way to convert the cumulative hazard to all sort of quantities
  chz_to_surv <- function(x) exp(-x)
  surv_to_chz <- function(x) -1 * log(x)
  survival <- function(chz) exp(-chz)
  survival_m <- function(chz) laplace_transform(chz, est_dist)
  cumhaz_m <- function(chz) -1 * log(laplace_transform(chz, est_dist))

  scenarios <- expand.grid(quantity, type, conf_int)

  res <- lapply(res, function(x) {

    # bounds is made into a column because cbind does not keep attributes
    bounds <- attr(x, "bounds")
    if("survival" %in% quantity & "conditional" %in% type) {
      x <- cbind(x,
                 as.data.frame(lapply(x[bounds], survival),
                               col.names = sub("cumhaz", "survival", bounds)))
    }

    if("survival" %in% quantity & "marginal" %in% type) {
      x <- cbind(x,
                 as.data.frame(lapply(x[bounds], survival_m),
                               col.names = sub("cumhaz", "survival_m", bounds)))
    }

    if("cumhaz" %in% quantity & "marginal" %in% type) {
      x <- cbind(x,
                 as.data.frame(lapply(x[bounds], cumhaz_m),
                               col.names = sub("cumhaz", "cumhaz_m", bounds)))
    }

    if(!("cumhaz" %in% quantity & "conditional" %in% type)) {
      cols_to_keep <- which(!(names(x) %in% c(bounds, "se_logH", "se_logH_adj")))
      x <- x[,cols_to_keep]
    }
    x

  })


  if(length(res) == 1) res <- res[[1]]
  res
}
