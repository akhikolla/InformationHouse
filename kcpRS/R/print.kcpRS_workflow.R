
#' @rdname kcpRS_workflow
#' @export

print.kcpRS_workflow<-function(x,...){
  res_kcpMean=x$kcpMean
  res_kcpVar=x$kcpVar
  res_kcpAR=x$kcpAR
  res_kcpCorr=x$kcpCorr

  RMean=ifelse(class(res_kcpMean)=="kcpRS",1,0)
  RVar=ifelse(class(res_kcpVar)=="kcpRS",1,0)
  RAR=ifelse(class(res_kcpAR)=="kcpRS",1,0)
  RCorr=ifelse(class(res_kcpCorr)=="kcpRS",1,0)

  ntests=RMean+RVar+RAR+RCorr

  if (ntests==0){warning("No running statistic selected.","\n")}

  if (ntests>0){
            if(RMean==1){
              cat("\n")
              cat("    KCP-Mean:","\n")
              print(res_kcpMean,kcp_details=FALSE) #kcp_details is set to FALSE: suppress other outputs.
              cat("    ===============================================================================================","\n")
              }
            if(RVar==1){
              cat("\n")
              cat("    KCP-Var:","\n")
              print(res_kcpVar,kcp_details=FALSE) #kcp_details is set to FALSE: suppress other outputs.
              cat("    ===============================================================================================","\n")
              }
            if(RAR==1){
              cat("\n")
              cat("    KCP-AR:","\n")
              print(res_kcpAR,kcp_details=FALSE) #kcp_details is set to FALSE: suppress other outputs.
              cat("    ===============================================================================================","\n")
              }
            if(RCorr==1){
              cat("\n")
              cat("    KCP-Corr:","\n")
              print(res_kcpCorr,kcp_details=FALSE) #kcp_details is set to FALSE: suppress other outputs.
              cat("    ===============================================================================================","\n")
            }
    }
}
