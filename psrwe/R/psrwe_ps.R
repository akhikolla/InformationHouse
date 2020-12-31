#' Get propensity scores
#'
#' Calculate propensity scores using logistic regression or random forest model
#'
#' @param data Dataframe with group assignment and covariates
#' @param ps_fml Propensity score formula. If null, all covariates will be
#'     included in the PS model in a linear form.
#' @param v_grp Column name corresponding to group assignment
#' @param cur_grp_level Group level for the current study. Default 1.
#'
#' @param v_covs Vector of column names corresponding to covariates
#' @param v_arm Column name corresponding treatment vs. control. Ignored for
#'     single arm studies. two-arm randomized studies.
#' @param trt_arm_level Arm level for the treatment arm. Ignored for single-arm
#'     studies.
#' @param nstrata Number of PS strata to be created
#' @param ... Parameters to get propensity scores. \itemize{
#'     \item{ps_method}{\code{logistic} forlogistic regression or
#'     \code{randomforest} for randomforest) \item{...}{Parameters for
#'     \code{randomForest}}} }
#'
#' @return A class \code{RWE_DWITHPS} list with items:
#'
#' \itemize{
#' \item{data}{Original data with column \code{_ps_} for estimated PS scores
#' and \code{_strata_} for PS stratum added}
#' \item{ps_fml}{PS formula for estimated PS scores}
#' \item{nstrata}{Number of strata}
#' }
#'
#' @examples
#'
#' rwe_ps(ex_dta, v_covs = paste("V", 1:7, sep = ""), v_grp = "Group",
#'        cur_grp_level = "current")
#'
#' @export
#'
rwe_ps <- function(data, ps_fml = NULL,
                   v_covs  = "V1",
                   v_grp   = "Group",
                   cur_grp_level = 1,
                   v_arm   = NULL,
                   trt_arm_level = 1,
                   nstrata = 5, ...) {

    dnames <- colnames(data)
    stopifnot(v_grp %in% dnames)

    ## generate formula
    if (is.null(ps_fml))
        ps_fml <- as.formula(paste(v_grp, "~",
                                   paste(v_covs, collapse = "+"),
                                   sep = ""))

    ## current study index will be kept in the results
    d1_inx   <- cur_grp_level == data[[v_grp]]
    keep_inx <- which(d1_inx)
    stopifnot(length(keep_inx) > 0)

    ## for 2-arm studies only
    if (!is.null(v_arm))
        d1_inx <- d1_inx & trt_arm_level == data[[v_arm]]

    ## get ps
    all_ps  <- get_ps(data, ps_fml = ps_fml, ...)
    d1_ps   <- all_ps[which(d1_inx)]

    ## add columns to data
    grp     <- rep(1, nrow(data));
    grp[which(data[[v_grp]] != cur_grp_level)] <- 0

    data[["_ps_"]]   <- all_ps
    data[["_grp_"]]  <- grp
    if (!is.null(v_arm))
        data[["_arm_"]]  <- data[[v_arm]]

    ## stratification
    if (nstrata > 0) {
        strata  <- rwe_cut(d1_ps, all_ps,
                           breaks = nstrata,
                           keep_inx = keep_inx)

        data[["_strata_"]] <- strata;
    }

    ## return
    rst <- list(data    = data,
                ps_fml  = ps_fml,
                nstrata = nstrata)

    class(rst) <- get_rwe_class("DWITHPS")
    rst
}


#' Summarize PS scores
#'
#' Get number of subjects and the distances of PS distributions for each PS
#' strata
#'
#' @inheritParams rwe_ps
#'
#' @param data_withps A class \code{RWE_DWITHPS} list. See \code{\link{rwe_ps}}.
#' @param min_n0 threshold for number of external subjects, below which the
#'     external data in the current stratum will be ignored by setting the PS
#'     distance to 0. Default value 10.
#' @param ... Parameters for \code{get_distance}, e.g., \code{metric} with
#'     options such as overlapping area (\code{ovl}).
#'
#' @return A class \code{RWE_PSDIST} dataframe with columns
#'   \itemize{
#' \item{Strata}{Index of stratum. 0 represents the overall information}
#' \item{N0,N1}{Number of subjects in group 0 and 1}
#' \item{N00, N10}{Number of arm 0 subjects in group 0 and 1, when arm exists}
#' \item{N01, N11}{Number of arm 1 subjects in group 0 and 1, when arm exists}
#' \item{Dist} Distance}
#'
#'
#' @examples
#'
#' dta_ps <- rwe_ps(ex_dta,
#'                  v_covs = paste("V", 1:7, sep = ""),
#'                  v_grp = "Group",
#'                  cur_grp_level = "current")
#'
#'
#' rwe_ps_dist(dta_ps, metric = "ovl")
#'
#'
#' @export
#'
rwe_ps_dist <- function(data_withps, min_n0 = 10,
                        trt_arm_level = NULL,
                        ...) {

    f_narm <- function(inx, dataps) {
        if (is.null(dataps[["_arm_"]]))
            return(c(length(inx), 0, 0))
        n0 <- length(which(0 == dataps[inx, "_arm_"]))
        n1 <- length(which(1 == dataps[inx, "_arm_"]))
        c(length(inx), n0, n1)
    }

    stopifnot(inherits(data_withps,
                       what = get_rwe_class("DWITHPS")))

    dataps   <- data_withps$data;
    nstrata  <- data_withps$nstrata;
    rst      <- NULL;
    for (i in seq_len(nstrata)) {
        inx_ps0 <- i == dataps[["_strata_"]] & 0 == dataps[["_grp_"]]
        inx_ps1 <- i == dataps[["_strata_"]] & 1 == dataps[["_grp_"]];
        n0_01   <- f_narm(which(inx_ps0), dataps)
        n1_01   <- f_narm(which(inx_ps1), dataps)

        if (!is.null(trt_arm_level) &
            !is.null(dataps[["_arm_"]])) {
            inx_ps0 <- inx_ps0 & trt_arm_level == dataps[["_arm_"]]
            inx_ps1 <- inx_ps1 & trt_arm_level == dataps[["_arm_"]]
        }

        ps0 <- dataps[which(inx_ps0), "_ps_"];
        ps1 <- dataps[which(inx_ps1), "_ps_"];

        if (0 == length(ps0) | 0 == length(ps1))
            warning("No samples in strata");

        if (any(is.na(c(ps0, ps1))))
            warning("NA found in propensity scores in a strata");

        if (length(ps0) < min_n0) {
            warning("Not enough data in the external data
                     in the current stratum.
                     External data ignored.");
            cur_dist <- 0;
        } else {
            cur_dist <- get_distance(ps0, ps1, ...)
        }

        rst <- rbind(rst, c(i, n0_01, n1_01, cur_dist))
    }

    ## overall
    inx_tot_ps0 <- which(0 == dataps[["_grp_"]])
    inx_tot_ps1 <- which(1 == dataps[["_grp_"]])
    n0_tot_01   <- f_narm(inx_tot_ps0, dataps)
    n1_tot_01   <- f_narm(inx_tot_ps1, dataps)

    ps0         <- dataps[inx_tot_ps0, "_ps_"];
    ps1         <- dataps[inx_tot_ps1, "_ps_"];
    all_dist    <- get_distance(ps0, ps1, ...)
    rst         <- rbind(rst, c(0, n0_tot_01, n1_tot_01, all_dist))

    ##return
    colnames(rst) <- c("Strata",
                       "N0", "N00", "N01",
                       "N1", "N10", "N11",
                       "Dist")
    rst           <- data.frame(rst)
    class(rst)    <- append(get_rwe_class("PSDIST"),
                            class(rst))

    rst
}


#' PS matching
#'
#' Match patients in external data source with patients in current study based
#' on PS using nearest neighbor method
#'
#' @inheritParams rwe_ps
#'
#' @param dta_cur current study dataframe
#' @param dta_ext external data source dataframe
#' @param ratio matching ratio
#'
#' @return A matrix with two columns. Column \code{pid}: id in dta_cur; Column
#'     \code{match_id}: matched id from dta_ext.
#'
#' @examples
#' dta_cur <- ex_dta[which(ex_dta$Group == "current"), ]
#' dta_ext <- ex_dta[which(ex_dta$Group != "current"), ]
#' rwe_ps_match(dta_cur, dta_ext, v_covs = paste("V", 1:7, sep = ""),
#'              ratio = 2)
#'
#' @export
#'
rwe_ps_match <- function(dta_cur, dta_ext, ratio = 3, ps_fml = NULL,
                         v_covs  = "V1") {

    n_cur <- nrow(dta_cur)
    n_ext <- nrow(dta_ext)
    ratio <- min(ratio, floor(n_ext / n_cur))

    dta_cur$grp_tmp <- 1
    dta_ext$grp_tmp <- 0
    dta_ext         <- dta_ext[, colnames(dta_cur)]

    dta    <- rbind(dta_cur, dta_ext)
    dta_ps <- rwe_ps(data = dta,
                     v_grp = "grp_tmp",
                     v_covs = v_covs,
                     nstrata = 0)

    ps        <- dta_ps$data[["_ps_"]]
    target    <- ps[1 : n_cur]
    candidate <- ps[- (1 : n_cur)]

    ## match
    rst <- c_match(target, candidate, ratio = ratio)[]
    rst <- rst[1:(ratio * n_cur)] + 1

    ## return
    cbind(pid      = rep(1:n_cur, each = ratio),
          match_id = rst)
}

#' Get number of subjects borrowed from each statum
#'
#' Based on PS distances, split the total number of subjects to be borrowed from
#' the external data source to each stratum
#'
#' @param total_borrow Total number of subjects to be borrowed
#' @param dta_ps_dist A class \code{RWE_PSDIST}
#' @param ... Other options
#'
#' @return A vector of number of subjects to be borrowed from each stratum
#'
#' @examples
#'
#' dta_ps <- rwe_ps(ex_dta,
#'                  v_covs = paste("V", 1:7, sep = ""),
#'                  v_grp = "Group",
#'                  cur_grp_level = "current")
#'
#' ps_dist <- rwe_ps_dist(dta_ps, metric = "ovl")
#' ps_borrow <- rwe_ps_borrow(total_borrow = 100, ps_dist)
#'
#' @export
#'
#'
rwe_ps_borrow <- function(total_borrow, dta_ps_dist, ...) {

    stopifnot(inherits(dta_ps_dist,
                       what = get_rwe_class("PSDIST")))

    ps_dist <- dta_ps_dist[seq_len(nrow(dta_ps_dist) - 1), ]
    rst     <- get_aborrow(total_borrow,
                           rs = ps_dist$Dist,
                           ns0 = ps_dist$N0, ...)

    rst
}


#' Plot PS distributions
#'
#' S3 method for visualizing PS adjustment
#'
#' @param x Class \code{RWE_DWITHPS} created by \code{\link{rwe_ps}}
#' @param plot_type Types of plots. \itemize{\item{ps}{PS density plot}
#'     \item{balance}{Covariate balance plot}}
#' @param ... Additional parameter for the plot
#'
#' @method plot RWE_DWITHPS
#'
#' @export
#'
plot.RWE_DWITHPS <- function(x, plot_type = c("ps", "balance"), ...) {
    type <- match.arg(plot_type)
    switch(type,
           ps = plot_ps(x, ...),
           balance = plot_balance(x, ...))
}


#' Create strata
#'
#' Cut a sequence of numbers into bins.
#'
#' @description
#'
#' The cut points are chosen such that there will with equal numbers in each bin
#' for \code{x}. By default, values of \code{y} that are outside the range of
#' \code{x} will be excluded from the bins, unless they are in the
#' \code{keep_inx}.
#'
#' @param x Vector of values based on which cut points will be determined
#' @param y Vector of values to be cut, default to be the same as \code{x}
#' @param breaks Number of cut points
#' @param keep_inx Indices of y that will be categorized as 1 or the largest bin
#'     even if their values are out of range of x, i.e. the y's that will not be
#'     trimmed
#'
#' @return A vector of stratum assignment for \code{y}. The y's that are outside
#'     the range of \code{x} and not in \code{keep_inx} are assigned \code{NA}
#'     in the result.
#' @examples
#'
#' x <- rnorm(100,  mean = 0, sd = 1)
#' y <- rnorm(1000, mean = 1, sd = 2)
#' rwe_cut(x, y, breaks = 5)
#'
#' @export
#'
#'
rwe_cut <- function(x, y = x, breaks = 5, keep_inx = NULL) {
    cuts <- quantile(x, seq(0, 1, length = breaks + 1))
    cuts[1] <- cuts[1] - 0.001
    rst <- rep(NA, length(y))
    for (i in 2:length(cuts)) {
        inx <- which(y > cuts[i - 1] & y <= cuts[i])
        rst[inx] <- i - 1
    }

    if (!is.null(keep_inx)) {
        inx <- which(y[keep_inx] <= cuts[1])
        if (0 < length(inx)) {
            rst[keep_inx[inx]] <- 1
        }

        inx <- which(y[keep_inx] > cuts[length(cuts)])
        if (0 < length(inx)) {
            rst[keep_inx[inx]] <- length(cuts) - 1
        }
    }

    rst
}
