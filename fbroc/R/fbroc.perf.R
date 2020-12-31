#' @export
perf.roc <- function(roc, ...) {
  .Deprecated("perf")
  perf(roc, ...)
}

#' Calculate performance for bootstrapped ROC curve
#'
#' Calculates different performance metric for ROC curves based on the bootstrap
#' results saved in an object of class \code{fbroc.roc}. Confidence intervals
#' are included.
#'
#' @param roc An object of class \code{fbroc.roc}.
#' @param metric A performance metric. Select "auc" for the AUC, "partial.auc" for the partial AUC, 
#' "tpr" for the TPR at a fixed FPR and "fpr" for the FPR at a fixed TPR.
#' @param conf.level The confidence level of the confidence interval.
#' @param fpr The fixed FPR at which the TPR is to be evaluated when \code{tpr} is selected as metric.
#' If partial AUC is investigated, then an FPR interval over which the partial area is to be calculated.
#' @param tpr The fixed TPR at which the FPR is to be evaluated when \code{fpr} is selected as metric.
#' If partial AUC is investigated, then an TPR interval over which the partial area is to be calculated.
#' @param correct.partial.auc Corrects partial AUC for easier interpretation using McClish correction.
#' Details are given below. Defaults to TRUE.
#' @param show.partial.auc.warning Whether to give warnings for partial AUCs below 0.5. Defaults to
#' true.
#' @param ... Further arguments, that are not used at this time.
#' @return A list of class \code{fbroc.perf}, containing the elements:
#' \item{Observed.Performance}{The observed performance.}
#' \item{CI.Performance}{Quantile based confidence interval for the performance.}
#' \item{conf.level}{Confidence level of the confidence interval.}
#' \item{metric}{Used performance metric.}
#' \item{params}{Parameters used to further specifiy metric, e.g. fixed TPR.}
#' \item{n.boot}{Number of bootstrap replicates used.}
#' \item{boot.results}{Performance in each bootstrap replicate.}
#' @seealso \code{\link{boot.roc}}, \code{\link{print.fbroc.perf}}, 
#'   \code{\link{plot.fbroc.perf}}
#' @template partial.auc.doc 
#' @examples
#' y <- rep(c(TRUE, FALSE), each = 100)
#' x <- rnorm(200) + y
#' result.boot <- boot.roc(x, y, n.boot = 100)
#' perf(result.boot, "auc")
#' perf(result.boot, "auc", conf.level = 0.99)
#' perf(result.boot, "partial.auc", fpr = c(0, 0.25), show.partial.auc.warning = FALSE)
#' @export
perf.fbroc.roc <- function(roc, metric = "auc", conf.level = 0.95, tpr = NULL, fpr = NULL, 
                           correct.partial.auc = TRUE, show.partial.auc.warning = TRUE, ...) {

  # start with data validation
  if (!is(roc, "fbroc.roc"))
    stop("roc must be of class fbroc.roc")
  if (length(metric) != 1 | class(metric) != "character")
    stop("metric must be character")
  if (!(metric %in% c("auc", "tpr", "fpr", "partial.auc")))
    stop(paste(metric,"is not a valid performance metric"))
  
  if (metric == "auc") {
    metric.text = "AUC"
    metric.number <- as.integer(0)
    param.vec <- 0
  }
  if (metric == "tpr") {
    param.vec = validate.single.numeric(fpr, "FPR")
    metric.text <- paste("TPR at a fixed FPR of", round(param.vec, 3))
    metric.number <- as.integer(1)
  }
  if (metric == "fpr") {
    param.vec = validate.single.numeric(tpr, "TPR")
    metric.text <- paste("FPR at a fixed TPR of", round(param.vec, 3))
    metric.number <- as.integer(2)
  }
  if (metric == "partial.auc") {
    if (is.null(tpr) & is.null(fpr)) stop("Either FPR or TPR must be specified")
    if (is.null(fpr)) {
      tpr.area = TRUE
      param <- sort(tpr)
      metric.number <- as.integer(4)
      metric.text <- paste("Partial AUC over TPR range [", round(param[1], 2),
                           ",", round(param[2], 2), "]", sep = "")
    }
    else {
      tpr.area = FALSE
      param <- sort(fpr)
      metric.number <- as.integer(3)
      metric.text <- paste("Partial AUC over FPR range [", round(param[1], 2),
                           ",", round(param[2], 2), "]", sep = "")
    }
    if (correct.partial.auc) metric.text <- paste("Corrected", metric.text) 
    if (length(param) != 2) stop("Interval must be given as a vector of length 2")
    if ((min(param) < 0) | (max(param) > 1)) stop("Interval values must be between 0 and 1")
    if (param[1] == param[2]) stop("Interval has width zero!")
    param.vec <- param
  }
  
  # call C++ to calculate actual results
  tpr.m <- matrix(roc$roc$TPR, nrow = 1)
  fpr.m <- matrix(roc$roc$FPR, nrow = 1)
  
  observed.perf <- get_cached_perf(tpr.m, fpr.m, param.vec, metric.number)
  #print(observed.perf)
  if (roc$use.cache) {
    perf.boot <- get_cached_perf(roc$boot.tpr, roc$boot.fpr, param.vec, metric.number)
  } else {
    perf.boot <- get_uncached_perf(roc$predictions,
                                   as.integer(roc$true.classes),
                                   param.vec,
                                   roc$n.boot,
                                   metric.number)
  }
  
  if (correct.partial.auc & metric == "partial.auc") {
    observed.perf <- partial.auc.index(observed.perf, param.vec, show.partial.auc.warning, tpr.area)
    perf.boot <- partial.auc.index(perf.boot, param.vec, show.partial.auc.warning, tpr.area)
  }
  
  # Quantile based confidence interval
  alpha <- 0.5 * (1 - conf.level)
  alpha.levels <- c(alpha, 1 - alpha) 
  ci <- as.numeric(quantile(perf.boot, alpha.levels))
  names(ci) <- NULL
  
  perf <- list(Observed.Performance = observed.perf,
               CI.Performance = ci,
               conf.level = conf.level,
               metric = metric.text,
               params = param.vec,
               n.boot = roc$n.boot,
               boot.results = perf.boot
  )
  class(perf) <- append(class(perf), "fbroc.perf")
  return(perf)
}