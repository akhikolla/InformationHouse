#' @rdname kcpRS_workflow
#' @param object An object of the type produced by \code{kcpRS_workflow}
#' @export

summary.kcpRS_workflow<-function(object,...){

  res_kcpMean=object$kcpMean
  res_kcpVar=object$kcpVar
  res_kcpAR=object$kcpAR
  res_kcpCorr=object$kcpCorr

  RMean=ifelse(class(res_kcpMean)=="kcpRS",1,0)
  RVar=ifelse(class(res_kcpVar)=="kcpRS",1,0)
  RAR=ifelse(class(res_kcpAR)=="kcpRS",1,0)
  RCorr=ifelse(class(res_kcpCorr)=="kcpRS",1,0)

  ntests=RMean+RVar+RAR+RCorr

  if (ntests==0){cat("No running statistic selected.","\n")}

  names=NULL
  if (ntests>0){
    if(RMean==1){names=c(names,as.character(res_kcpMean$RS_name))}
    if(RVar==1){names=c(names,as.character(res_kcpVar$RS_name))}
    if(RAR==1){names=c(names,as.character(res_kcpAR$RS_name))}
    if(RCorr==1){names=c(names,as.character(res_kcpCorr$RS_name))}
    }


  univ_RS=0
  while (univ_RS==0){                                    #check if there is a univariate stat monitores
                        if(RMean==1){x_for_sum=res_kcpMean}
                        univ_RS=1

                        if(RVar==1){x_for_sum=res_kcpVar}
                        univ_RS=1

                        if(RAR==1){x_for_sum=res_kcpAR}
                        univ_RS=1}

   if (univ_RS==0){x_for_sum=res_kcpCorr}


                        cat("\n")
                        cat("===============================================================================================","\n")
                        cat("SETTINGS:","\n")
                        cat("    Running statistics monitored:", "\n")
                                if (RMean==1){cat("       Mean", "\n")}
                                if (RVar==1){cat("       Variance", "\n")}
                                if (RAR==1){cat("       Autocorrelation", "\n")}
                                if (RCorr==1){cat("       Correlation", "\n")}

                        if (univ_RS==1){cat("    Number of variables monitored:", ncol(x_for_sum$RS), "\n")}
                        cat("    Number of correlations monitored:", ncol(res_kcpCorr$RS), "\n")


                        cat("    Selected window size:", x_for_sum$wsize , "\n")
                        cat("    Number of time windows:", nrow(x_for_sum$RS) , "\n")
                        cat("    Selected maximum number of change points:", ncol(x_for_sum$CPs_given_K)-2 , "\n")
                        cat("\n")

                        cat("    Permutation test:", ifelse(x_for_sum$nperm>0, "Yes", "No") , "\n")
                        cat("        Number of permuted data sets used:", ifelse(x_for_sum$nperm>0, x_for_sum$nperm, NA) , "\n")

                        cat("\n")
                        cat("\n")
                        cat("===============================================================================================","\n")
                        cat("OUTPUT:")
                        cat("\n")
  print(object)
}
