#' @rdname kcpRS_workflow
#' @param x An object of the type produced by \code{kcpRS_workflow}
#' @param ... Further plotting arguments
#' @importFrom graphics plot
#' @importFrom RColorBrewer brewer.pal
#' @export

plot.kcpRS_workflow<-function(x,...){


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

    par(mfrow=c(ntests,1))
      if(RMean==1){plot(res_kcpMean,title="Running Means")}
      if(RVar==1){plot(res_kcpVar,title="Running Variances")}
      if(RAR==1){plot(res_kcpAR,title="Running Autocorrelations")}
      if(RCorr==1){plot(res_kcpCorr,title="Running Correlations")}
  }

}
