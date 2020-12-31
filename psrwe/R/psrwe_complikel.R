#' PS-Integrated Composite Likelihood Estimation
#'
#' Estimate the mean of the outcome based on PS-integrated composite likelihood
#' approach. Variance is estimated by Jack-Knife method. Applies to the case
#' when there is only one external data source.
#'
#' @inheritParams rwe_ps_dist
#' @inheritParams rwe_ps_powerp
#'
#' @param v_borrow Vector of number of subjects to be borrowed from each stratum
#' @param ... Parameters for \code{rwe_cl}
#'
#' @return A dataframe with class name \code{RWE_CLRST}. It contains the
#'     composite estimation of the mean for each stratum as well as the
#'     jackknife estimation for each subject. The results should be further
#'     summarized by its S3 method \code{summary}.
#'
#' @examples
#' \donttest{
#' dta_ps <- rwe_ps(ex_dta,
#'                  v_covs = paste("V", 1:7, sep = ""),
#'                  v_grp = "Group",
#'                  cur_grp_level = "current",
#'                  nstrata = 5)
#' ps_dist <- rwe_ps_dist(dta_ps)
#' ps_borrow <- rwe_ps_borrow(total_borrow = 40, ps_dist)
#' rst_cl    <- rwe_ps_cl(dta_ps, v_borrow = ps_borrow)
#' summary(rst_cl)}
#'
#' @export
#'
rwe_ps_cl <- function(data_withps, v_borrow = 0, v_outcome = "Y", ...) {
    stopifnot(inherits(data_withps,
                       what = get_rwe_class("DWITHPS")))

    data <- data_withps$data
    stopifnot(v_outcome %in% colnames(data))

    data[["_id_"]] <- seq_len(nrow(data))

    ## prepare data
    data    <- data[!is.na(data[["_strata_"]]), ]
    nstrata <- max(data[["_strata_"]])

    if (length(v_borrow) != nstrata) {
        warning("Length of v_borrow is different from the number of strata")
        v_borrow <- rep(v_borrow, length.out = nstrata)
    }

    ## find mwle
    rst_theta <- NULL
    for (i in seq_len(nstrata)) {
        cur_01 <- get_cur_d(data, i, c(v_outcome, "_id_"))
        cur_d1 <- cur_01$cur_d1
        cur_d0 <- cur_01$cur_d0

        ns1    <- nrow(cur_d1)
        ns0    <- nrow(cur_d0)

        cur_borrow <- v_borrow[i]
        cur_theta  <- rwe_cl(dta_cur  = cur_d1[[v_outcome]],
                             dta_ext  = cur_d0[[v_outcome]],
                             n_borrow = cur_borrow,
                             ...)
        rst_theta  <- rbind(rst_theta,
                            c(i, NA, NA, ns1, ns0, cur_theta))

        ##jackknife
        for (j in 1:ns1) {
            cur_jk <- rwe_cl(dta_cur = cur_d1[-j, v_outcome],
                             dta_ext = cur_d0[[v_outcome]],
                             n_borrow = cur_borrow, ...)

            rst_theta  <- rbind(rst_theta,
                                c(i, 1, cur_d1[j, "_id_"],
                                  ns1 - 1, ns0, cur_jk))
        }

        if (ns0 > 0) {
            for (j in 1:ns0) {
                cur_jk <- rwe_cl(dta_cur = cur_d1[[v_outcome]],
                                 dta_ext = cur_d0[-j, v_outcome],
                                 n_borrow = cur_borrow, ...)

                rst_theta <- rbind(rst_theta,
                                   c(i, 0, cur_d0[j, "_id_"],
                                     ns1, ns0 - 1, cur_jk))
            }
        }
    }

    ## add information
    colnames(rst_theta) <- c("Strata", "Group", "ID", "N1", "N0", "Theta")
    rst_theta           <- data.frame(rst_theta)
    class(rst_theta)    <- append(get_rwe_class("CLRST"),
                                  class(rst_theta))


    ## return
    rst_theta
}


#' Composite Likelihood Estimation
#'
#' Estimate parameter of interest based composite likelihood for a single PS
#' stratum
#'
#' @inheritParams rwe_ps_powerp
#'
#' @param dta_cur Vector of outcome from a PS stratum in current study
#' @param dta_ext Vector of outcome from a PS stratum in external data source
#' @param n_borrow Number of subjects to be borrowed
#' @param equal_sd Boolean. whether sd is the same between the current study and
#'     external data source
#'
#' @return Maximum composite likelihood estimator of the mean
#'
#' @examples
#' x <- rnorm(100,  mean = 0, sd = 1)
#' y <- rnorm(1000, mean = 1, sd = 2)
#' rwe_cl(x, y, n_borrow = 20, equal_sd = FALSE)
#'
#' @export
#'
rwe_cl <- function(dta_cur, dta_ext, n_borrow = 0,
                   outcome_type = c("continuous", "binary"),
                   equal_sd = TRUE) {

    f_ll <- function(pars) {
        theta  <- pars[1]
        sig2_1 <- pars[2]
        sig2_0 <- pars[3]

        ll <- - n1 * log(sig2_1) / 2
        ll <- ll - n1 * mean((dta_cur - theta)^2) / 2 / sig2_1
        ll <- ll - n_borrow * log(sig2_0) / 2
        ll <- ll - n_borrow * mean((dta_ext - theta)^2) / 2 / sig2_0

        ll
    }

    f_gradient <- function(pars) {
        theta  <- pars[1]
        sig2_1 <- pars[2]
        sig2_0 <- pars[3]

        g <- numeric(length(pars))

        ## d logl / d theta
        g[1] <- n1 / sig2_1 * (mean(dta_cur) - theta) +
          n_borrow / sig2_0 * (mean(dta_ext) - theta)

        ## d logl / d sig2.1
        g[2] <- - n1 / 2 / sig2_1 +
          n1 * mean((dta_cur - theta)^2) / 2 / sig2_1 / sig2_1
        ## d logl / d sig2.0
        g[3] <- - n_borrow / 2 / sig2_0 +
          n_borrow * mean((dta_cur - theta)^2) / 2 / sig2_0 / sig2_0

        return(g)
    }


    type <- match.arg(outcome_type)
    n1   <- length(dta_cur)

    ## ignore external data
    if (0 == n_borrow) {
        ## placeholder
        dta_ext  <- dta_cur
        equal_sd <- TRUE;
    }

    init_theta <- (n1 / (n1 + n_borrow)) * mean(dta_cur) +
        (n_borrow / (n1 + n_borrow)) * mean(dta_ext)

    if (("continuous" == type & equal_sd) |
        "binary" == type) {
        rst <- init_theta
    } else {
        init_sig2_1 <- mean((dta_cur - init_theta)^2)
        init_sig2_0 <- mean((dta_ext - init_theta)^2)
        rst         <- optim(c(init_theta, init_sig2_1, init_sig2_0),
                             method  = "L-BFGS-B",
                             fn      = f_ll,
                             lower   = c(-Inf, 1e-6, 1e-6),
                             upper   = rep(Inf, 3),
                             control = list(fnscale = -1))$par[1];
    }

    rst
}

#' Summary of Composite Likelihood Estimation
#'
#' Summarize the result in a class \code{RWE_CLRST} objects
#' generated by \code{\link{rwe_ps_cl}}.
#'
#' @param object A class \code{RWE_CLRST} object
#' @param ... Extra arguments
#'
#'
#' @return
#'
#' Summary of composite likelihood estimation
#'
#' @method summary RWE_CLRST
#'
#' @export
#'
#'
summary.RWE_CLRST <- function(object, ...) {
    f_overall <- function(dta) {
        overall_theta <- dta$N1 * dta$Theta
        overall_theta <- sum(overall_theta) / sum(dta$N1)
    }

    f_jk <- function(jk_all, overall_mean) {
        var_theta <- (length(jk_all) - 1) / length(jk_all)
        var_theta <- var_theta * sum((jk_all - overall_mean)^2)
        var_theta
    }

    ## stratum specific theta
    theta_strata <- object %>%
      dplyr::filter(is.na(.data$Group)) %>%
      dplyr::select(-.data$Group, -.data$ID) %>%
      arrange(.data$Strata)

    ## stratum specific jknife variance
    var_strata <- NULL
    for (i in seq_len(nrow(theta_strata))) {
        cur_s  <- theta_strata[i, ]
        cur_jk <- object %>%
          dplyr::filter(!is.na(.data$Group) &
                        .data$Strata == cur_s[1, "Strata"])

        cur_var <- f_jk(cur_jk$Theta, cur_s[1, "Theta"])
        var_strata <- c(var_strata, cur_var)
    }
    theta_strata$Variance <- var_strata

    ## overall mean
    overall_mean <- f_overall(theta_strata)

    ## overall jackknife variance
    jk_all <- NULL
    for (i in seq_len(nrow(object))) {
        cur_i <- object[i, ]
        if (is.na(cur_i["Group"]))
            next

        cur_d <- theta_strata %>%
          dplyr::filter(.data$Strata != cur_i[1, "Strata"]) %>%
          dplyr::select(-.data$Variance)


        cur_d <- rbind(cur_d,
                       cur_i[1, c("Strata", "N1", "N0", "Theta")])

        cur_jk <- f_overall(cur_d)
        jk_all <- c(jk_all, cur_jk)
    }

    var_theta <- f_jk(jk_all, overall_mean)

    ## return
    list(overall_mean       = overall_mean,
         jackknife_variance = var_theta,
         theta_by_stratum   = theta_strata)
}


#' PS-Integrated Composite Likelihood Estimation for Randomized Study
#'
#' Estimate the treatment effect based on PS-integrated composite likelihood
#' approach. Subjects from the external data source augment an arm, usually the
#' control arm.
#'
#' @inheritParams rwe_ps
#' @inheritParams rwe_ps_cl
#' @inheritParams rwe_ps_borrow
#'
#' @param ... Parameters for \code{\link{rwe_ps_cl}}
#'
#' @return A list with estimated treatment effect and the variance, as well as
#'     the estimated mean for the treatment and control arms
#'
#' @export
#'
rwe_ps_cl2arm <- function(data_withps, v_arm = "Arm", trt_arm_level = 1,
                          total_borrow = 0, ...) {

    stopifnot(inherits(data_withps,
                       what = get_rwe_class("DWITHPS")))
    data <- data_withps$data

    ## treatment arm
    stopifnot(v_arm %in% colnames(data))
    inx_trt <- trt_arm_level == data[[v_arm]]
    inx_g1  <- 1 == data[["_grp_"]]

    ## treatment arm
    dta_trt <- data[inx_trt & inx_g1, ]

    ## add fake data to avoid warnings about no subjects found in the external
    ## data source
    fake_dta            <- dta_trt
    fake_dta[["_grp_"]] <- 0

    dta_trt_withps      <- data_withps
    dta_trt_withps$data <- rbind(dta_trt, fake_dta)

    est_trt <- rwe_ps_cl(dta_trt_withps,
                         v_borrow = rep(0, data_withps$nstrata),
                         ...)

    sum_trt <- summary(est_trt)

    ## control arm
    dta_ctl <- data[!inx_trt, ]
    dta_ctl_withps <- data_withps
    dta_ctl_withps$data <- dta_ctl

    ps_dist   <- rwe_ps_dist(dta_ctl_withps)
    ps_borrow <- rwe_ps_borrow(total_borrow, ps_dist)
    est_ctl   <- rwe_ps_cl(dta_ctl_withps,
                           v_borrow = ps_borrow, ...)

    sum_ctl <- summary(est_ctl)

    ## treatment effect
    mu  <- sum_trt$theta_by_stratum$Theta    - sum_ctl$theta_by_stratum$Theta
    vs  <- sum_trt$theta_by_stratum$Variance - sum_ctl$theta_by_stratum$Variance
    nct <- sum_trt$theta_by_stratum$N1       + sum_ctl$theta_by_stratum$N1

    ws   <- nct / sum(nct)
    mu   <- sum(ws * mu)
    vs   <- sum(ws^2 * vs)

    ## return
    list(treatment = sum_trt,
         control   = sum_ctl,
         effect    = list(Estimate = mu, Variance = vs))
}
