#' @rdname kcpRS
#' @param object An object of the type produced by \code{kcpRS_workflow}
#' @export


summary.kcpRS<-function(object,...){

  cat("\n")
  cat("SETTINGS:","\n")
  cat("    Running statistics monitored:", object$RS_name, "\n")
  cat("    Number of running", paste0(object$RS_name,"s monitored:"), ncol(object$RS), "\n")
  cat("    Selected window size:", object$wsize , "\n")
  cat("    Number of time windows:", nrow(object$RS) , "\n")
  cat("    Selected maximum number of change points:", ncol(object$CPs_given_K)-2 , "\n")
  cat("\n")

  cat("    Permutation test:", ifelse(object$nperm>0, "Yes", "No") , "\n")
  cat("        Number of permuted data sets used:", ifelse(object$nperm>0, object$nperm, NA) , "\n")

  cat("\n")
  cat("\n")
  cat("OUTPUT:")
  print(object)
}

