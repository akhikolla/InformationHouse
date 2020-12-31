#' @export

kcpRS_workflow.default <- function(data,
                                   RS_funs = c("runMean", "runVar", "runAR", "runCorr"),
                                   wsize = 25,
                                   nperm = 1000,
                                   Kmax = 10,
                                   alpha = .05,
                                   varTest = FALSE,
                                   bcorr = TRUE,
                                   ncpu = 1
                                   ) {
  if (ncpu<=detectCores()){

    rm <- ifelse("runMean" %in% RS_funs, 1, 0)
    rv <- ifelse("runVar" %in% RS_funs, 1, 0)
    ra <- ifelse("runAR" %in% RS_funs, 1, 0)
    rc <- ifelse("runCorr" %in% RS_funs, 1, 0)

    kcp_mean <- NULL
    kcp_var <- NULL
    kcp_corr <- NULL
    kcp_AR <- NULL

    #check which tests are to be performed and if correction is asked
    ntest <- rm + rv + ra + rc #no. of tests
    alpha_per_test <-
      ifelse(isTRUE(bcorr), alpha / ntest, alpha)  #alpha per RS test


    #Running means
    if (rm == 1) {
      kcp_mean <- kcpRS(
        data,
        RS_fun = runMean,
        RS_name = "Mean",
        wsize,
        nperm,
        Kmax,
        alpha = alpha_per_test,
        varTest,
        ncpu
      )
    }
    ncp_mean <- length(kcp_mean$changePoints)
    if (rm == 1 &
        ncp_mean > 0 &
        (rv + ra + rc) > 0) {
      #if there is a mean change point and further tests are requested
      cps <- as.numeric(kcp_mean$changePoints)
      nv <- ncol(data)
      N <- nrow(data)
      bounds <- c(1, cps, N)
      nbounds <- length(bounds)
      dat_centered <- matrix(0, nrow <- N, ncol = nv)

      for (v in 1:nv) {
        for (k in 2:nbounds) {
          mean_temp <- mean(data[bounds[k - 1]:(bounds[k] - 1), v])
          dat_centered[bounds[k - 1]:(bounds[k] - 1), v] <-
            data[bounds[k - 1]:(bounds[k] - 1), v] - mean_temp
        }
      }

      dat_centered <- as.data.frame(dat_centered)
      colnames(dat_centered) <- colnames(data)
      data <- dat_centered
    }



    #Running var
    if (rv == 1) {
      kcp_var = kcpRS(
        data,
        RS_fun = runVar,
        RS_name = "Variance",
        wsize,
        nperm,
        Kmax,
        alpha = alpha_per_test,
        varTest,
        ncpu
      )
    }


    #Running AR
    if (ra == 1) {
      kcp_AR = kcpRS(
        data,
        RS_fun = runAR,
        RS_name = "Autocorrelation",
        wsize,
        nperm,
        Kmax,
        alpha = alpha_per_test,
        varTest,
        ncpu
      )
    }

    #Running corr
    if (rc == 1) {
      kcp_corr = kcpRS(
        data,
        RS_fun = runCorr,
        RS_name = "Correlation",
        wsize,
        nperm,
        Kmax,
        alpha = alpha_per_test,
        varTest,
        ncpu
      )
    }


    output <- list(
      "kcpMean" = kcp_mean,
      "kcpVar" = kcp_var,
      "kcpAR" = kcp_AR,
      "kcpCorr" = kcp_corr
    )

    class(output) <- "kcpRS_workflow"

    return(output)
  }
}
