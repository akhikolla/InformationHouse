#' Call STAN models
#'
#' Call STAN models. Called by \code{rwe_ps_powerp}.
#'
#' @param lst_data List of study data to be passed to STAN
#' @param stan_mdl STAN model including
#'    \describe{
#'     \item{powerps}{PS-power prior model for continuous outcomes}
#'     \item{powerpsbinary}{PS-power prior model for binary outcomes}
#'     \item{powerp}{Power prior model}
#' }
#'
#' @param chains STAN parameter. Number of Markov chainsm
#' @param iter STAN parameter. Number of iterations
#' @param warmup STAN parameter. Number of burnin.
#' @param control STAN parameter. See \code{rstan::stan} for details.
#' @param ... other options to call STAN sampling such as \code{thin},
#'     \code{algorithm}. See \code{rstan::sampling} for details.#'
#'
#' @return Result from STAN sampling
#'
#' @export
#'
rwe_stan <- function(lst_data,
                     stan_mdl = c("powerps", "powerpsbinary", "powerp"),
                     chains = 4, iter = 2000, warmup = 1000,
                     control = list(adapt_delta = 0.95), ...) {

    stan_mdl <- match.arg(stan_mdl)
    stan_rst <- rstan::sampling(stanmodels[[stan_mdl]],
                                data    = lst_data,
                                chains  = chains,
                                iter    = iter,
                                warmup  = warmup,
                                control = control,
                                ...)

    stan_rst
}


#' Get posterior samples based on PS-power prior approach
#'
#' Draw posterior samples of the parameters of interest for the PS-power prior
#' approach
#'
#' @inheritParams rwe_ps_dist
#' @inheritParams rwe_ps_borrow
#'
#' @param v_outcome Column name corresponding to the outcome
#' @param outcome_type Type of outcomes: \code{continuous} or \code{binary}
#' @param v_distance Vector of distances in PS distributions for each stratum.
#'     See \code{\link{rwe_ps_dist}}.
#' @param prior_type Whether treat power parameter as fixed (\code{fixed}) or
#'     fully Bayesian (\code{random})
#' @param seed Random seed
#' @param ... extra parameters for calling function \code{\link{rwe_stan}}
#'
#' @return A list with the following objects \describe{
#'     \item{post_theta}{Posterior samples of parameter of interest}
#'     \item{post_theta_stratum}{Posterior samples of parameter of interest in
#'     each stratum} \item{stan_rst}{Result from STAN sampling} }
#'
#' @examples
#' \donttest{
#'  dta_ps <- rwe_ps(ex_dta, v_covs = paste("V", 1:10, sep = ""),
#'                   v_grp = "group", cur_grp_level = "current")
#'  ps_dist <- rwe_ps_dist(dta_ps, metric = "ovl")
#'  post_smps <- rwe_ps_powerp(dta_ps, outcome_type = "binary",
#'                             total_borrow = 30,
#'                             v_distance = ps_dist$Dist[1:5],
#'                             v_outcome = "Y")}
#' @export
#'
rwe_ps_powerp <- function(data_withps, outcome_type = c("continuous", "binary"),
                          total_borrow = 0, v_distance = NULL,
                          prior_type = c("fixed", "random"),
                          v_outcome = "Y",  ..., seed = NULL) {

    stopifnot(inherits(data_withps,
                       what = get_rwe_class("DWITHPS")))

    type       <- match.arg(outcome_type)
    prior_type <- match.arg(prior_type)
    data       <- data_withps$data
    stopifnot(v_outcome %in% colnames(data))

    ## prepare stan data
    data    <- data[!is.na(data[["_strata_"]]), ]
    nstrata <- max(data[["_strata_"]]);
    stan_d  <- NULL;

    ## set random seed
    if (!is.null(seed))
        old_seed <- set.seed(seed)

    Y1      <- NULL;
    INX1    <- NULL;
    for (i in seq_len(nstrata)) {
        cur_01 <- get_cur_d(data, i, v_outcome)
        cur_d1 <- cur_01$cur_d1
        cur_d0 <- cur_01$cur_d0

        cur_n1 <- length(cur_d1)
        cur_d  <- c(N0    = length(cur_d0),
                    YBAR0 = mean(cur_d0),
                    SD0   = sd(cur_d0),
                    N1    = cur_n1,
                    YBAR1 = mean(cur_d1),
                    YSUM1 = sum(cur_d1));

        stan_d <- rbind(stan_d, cur_d)
        Y1     <- c(Y1, cur_d1)
        INX1   <- c(INX1, rep(i, length = cur_n1))
    }

    if (is.null(v_distance)) {
        v_distance <- rep(1 / nstrata, nstrata)
    } else {
        stopifnot(length(v_distance) == nstrata)
        v_distance <- v_distance / sum(v_distance)
    }

    lst_data  <- list(S     = nstrata,
                      A     = total_borrow,
                      RS    = v_distance,
                      FIXVS = as.numeric(prior_type == "fixed"),
                      N0    = stan_d[, "N0"],
                      N1    = stan_d[, "N1"],
                      YBAR0 = stan_d[, "YBAR0"],
                      SD0   = stan_d[, "SD0"])

    ## sampling
    rst_theta  <- NULL
    rst_thetas <- NULL
    if ("continuous" == type) {
        stan_mdl <- ifelse(1 == nstrata,  "powerp", "powerps")
        lst_data <- c(lst_data,
                      list(TN1  = length(Y1),
                           Y1   = Y1,
                           INX1 = INX1))
        rst_post <- rwe_stan(lst_data = lst_data, stan_mdl = stan_mdl, ...)
    } else {
        lst_data <- c(lst_data,
                      list(YBAR1 = as.numeric(stan_d[, "YBAR1"]),
                           YSUM1 = as.numeric(stan_d[, "YSUM1"])))

        if (1 < nstrata) {
            rst_post  <- rwe_stan(lst_data = lst_data,
                                  stan_mdl = "powerpsbinary",
                                 ...)
        } else {
            rst_post  <- NULL;
            rst_theta <- with(lst_data,
                              rbeta(2000, YSUM1 + A * YBAR0 + 0.5,
                                    N1 - YSUM1 + A * (1 - YBAR0) + 0.5))
        }
    }

    if (!is.null(rst_post)) {
        rst_theta  <- rstan::extract(rst_post, pars = "theta")$theta
        rst_thetas <- rstan::extract(rst_post, pars = "thetas")$thetas
    }

    ## reset random seed
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst <- list(post_theta         = rst_theta,
                post_theta_stratum = rst_thetas,
                stan_rst           = rst_post)

    class(rst) <- append(get_rwe_class("PPRST"),
                         class(rst))

    rst
}


#' Summary of PS-power prior Estimation
#'
#' Summarize the result in a class \code{RWE_POWERPRST} objects
#' generated by \code{\link{rwe_ps_cl}}.
#'
#' @param object A class \code{RWE_POWERPRST} object
#' @param ... Extra arguments
#'
#'
#' @return
#'
#' Summary of composite likelihood estimation
#'
#' @method summary RWE_POWERPRST
#'
#' @export
#'
#'
summary.RWE_POWERPRST <- function(object, ...) {
    f_overall <- function(dta) {
        c(mean(dta), var(dta))
    }

    overall <- f_overall(object$post_theta)
    post_stratum <- object$post_theta_stratum
    theta_strata <- NULL
    for (i in seq_len(ncol(post_stratum))) {
        theta_strata <- rbind(theta_strata,
                              c(i, f_overall(post_stratum[, i])))
    }
    colnames(theta_strata) <- c("Strata", "Theta", "Variance")

    ## return
    list(overall_mean      = overall[1],
         overall_variance  = overall[2],
         theta_by_stratum  = data.frame(theta_strata))
}
