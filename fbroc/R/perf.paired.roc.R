#' Calculate performance for paired bootstrapped ROC curves
#' 
#' For a given metric this calculates the difference in performance between two paired predictors
#' stored in an object of class \code{fbroc.paired.roc} in addition to their individual performance.
#' 
#' @param roc An object of class \code{fbroc.paired.roc}.
#' @inheritParams perf.fbroc.roc
#' @export
#' @template partial.auc.doc 
#' @examples
#' data(roc.examples)
#' example <- boot.paired.roc(roc.examples$Cont.Pred, roc.examples$Cont.Pred.Outlier,
#'                                roc.examples$True.Class, n.boot = 100)
#' perf(example, metric = "auc")   
#' # Get difference in TPR at a FPR of 20%   
#' perf(example, metric = "tpr", fpr = 0.2)    
#' perf(example, metric = "partial.auc", fpr = c(0, 0.25), 
#'      show.partial.auc.warning = FALSE)                       
perf.fbroc.paired.roc <- function(roc, metric = "auc", conf.level = 0.95, tpr = NULL, fpr = NULL, 
                                  correct.partial.auc = TRUE, show.partial.auc.warning = TRUE, ...){
  # start with data validation
  if (!is(roc, "fbroc.paired.roc"))
    stop("roc must be of class fbroc.paired.roc")
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
      param <- sort(tpr)
      tpr.area <- TRUE
      metric.number <- as.integer(4)
      metric.text <- paste("Partial AUC over TPR range [", round(param[1], 2),
                           ",", round(param[2], 2), "]", sep = "")
    }
    else {
      param <- sort(fpr)
      tpr.area <- FALSE
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
  tpr.m1 <- matrix(roc$roc1$TPR, nrow = 1)
  fpr.m1 <- matrix(roc$roc1$FPR, nrow = 1)
  tpr.m2 <- matrix(roc$roc2$TPR, nrow = 1)
  fpr.m2 <- matrix(roc$roc2$FPR, nrow = 1)
  
  observed.perf1 <- get_cached_perf(tpr.m1, fpr.m1, param.vec, metric.number)
  observed.perf2 <- get_cached_perf(tpr.m2, fpr.m2, param.vec, metric.number)
  if (roc$use.cache) {
    perf.boot.list <- get_cached_perf_paired(roc$boot.tpr1, roc$boot.fpr1, 
                                             roc$boot.tpr2, roc$boot.fpr2,
                                             param.vec, metric.number)
  } else {
    perf.boot.list <- get_uncached_perf_paired(roc$predictions1,
                                               roc$predictions2,
                                               as.integer(roc$true.classes),
                                               param.vec,
                                               roc$n.boot,
                                               metric.number)
  }
  
  if (correct.partial.auc & metric == "partial.auc") {
    observed.perf1 <- partial.auc.index(observed.perf1, param.vec, show.partial.auc.warning, 
                                        tpr.area)
    observed.perf2 <- partial.auc.index(observed.perf2, param.vec, show.partial.auc.warning, 
                                        tpr.area)
    perf.boot.list[[1]] <- partial.auc.index(perf.boot.list[[1]], param.vec, 
                                             show.partial.auc.warning, tpr.area)
    perf.boot.list[[2]] <- partial.auc.index(perf.boot.list[[2]], param.vec, 
                                             show.partial.auc.warning, tpr.area)
  }
  
  # Quantile based confidence interval
  alpha <- 0.5 * (1 - conf.level)
  alpha.levels <- c(alpha, 1 - alpha) 
  ci1 <- as.numeric(quantile(perf.boot.list[[1]], alpha.levels))
  names(ci1) <- NULL
  ci2 <- as.numeric(quantile(perf.boot.list[[2]], alpha.levels))
  names(ci2) <- NULL
  ci.diff <- as.numeric(quantile(perf.boot.list[[1]] - perf.boot.list[[2]], alpha.levels))
  names(ci.diff) <- NULL
  cor.pred <- cor(perf.boot.list[[1]], perf.boot.list[[2]])
  
  perf <- list(Observed.Performance.Predictor1 = observed.perf1,
               CI.Performance.Predictor1 = ci1,
               Observed.Performance.Predictor2 = observed.perf2,
               CI.Performance.Predictor2 = ci2,
               Observed.Difference = observed.perf1 - observed.perf2,
               CI.Performance.Difference = ci.diff,
               conf.level = conf.level,
               Cor = cor.pred,
               metric = metric.text,
               params = param.vec,
               n.boot = roc$n.boot,
               boot.results.pred1 = perf.boot.list[[1]],
               boot.results.pred2 = perf.boot.list[[2]]
  )
  class(perf) <- append(class(perf), "fbroc.perf.paired")
  return(perf)
}