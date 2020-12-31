
#Function to perform igcda algorithm for selected MNAR values and slsa algorithm or mle algorithm for MCAR values

impute.mix <-function (tab, prob.MCAR, threshold, conditions, repbio=NULL, reptech=NULL, methodMCAR="mle",
                       nknn=15, weight=1, selec="all", ind.comp=1, progress.bar=TRUE, q=0.95,
                       ncp.max=5,maxiter = 10, ntree = 100, variablewise = FALSE,
                       decreasing = FALSE, verbose = FALSE,
                       mtry = floor(sqrt(ncol(tab))), replace = TRUE,
                       classwt = NULL, cutoff = NULL, strata = NULL,
                       sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                       xtrue = NA, parallelize = c('no', 'variables', 'forests'),
                       methodMNAR="igcda", q.min = 0.025, q.norm = 3, eps = 0, distribution = "unif",
                       param1 = 3, param2 = 1, R.q.min=1){

  if (is.null(repbio)){repbio=as.factor(1:length(conditions));}
  if (is.null(reptech)){reptech=as.factor(1:length(conditions));}

  tab.mvs=tab;

  #Random draw of MCAR and MNAR values
  l.MCAR=matrix(0,nrow(tab),ncol(tab));
  l.MCAR[which(prob.MCAR>threshold)]=1;
  l.MCAR[which(!is.na(tab))]=0;

  if (methodMCAR == "mle") {
    tab.imp = impute.mle(tab=tab.mvs, conditions = conditions);
  }else {
    if (methodMCAR == "rf") {
      tab.imp = impute.RF(tab=tab.mvs, conditions = conditions,
                                             maxiter = maxiter, ntree = ntree, variablewise = variablewise,
                                             decreasing = decreasing, verbose = verbose,
                                             mtry = mtry, replace = replace,
                                             classwt = classwt, cutoff = cutoff, strata = strata,
                                             sampsize = sampsize, nodesize = nodesize, maxnodes = maxnodes,
                                             xtrue = xtrue, parallelize = parallelize);
    }else {
      if (methodMCAR == "pca") {
        tab.imp = impute.PCA(tab=tab.mvs, conditions = conditions,
                                                ncp.max=ncp.max);
      }else {
        tab.imp = impute.slsa(tab=tab.mvs, conditions = conditions, repbio = repbio,
                                                 reptech = reptech, nknn = nknn, weight = weight,
                                                 selec = selec, progress.bar = FALSE, ind.comp = ind.comp);
      }
    }
  }

  #Impute MCAR values
  tab.mvs[which(l.MCAR==1)]=tab.imp[which(l.MCAR==1)];

  #Impute remaining MNAR values

  #Replace remaining missing values using MNAR-devoted algorithms
  if (methodMNAR == "igcda"){
    tab.mvs.imp = impute.igcda(tab = tab.mvs, tab.imp = tab.imp,
                               conditions = conditions, q = q)
  }else{
    tab.mvs.imp = impute.pa(tab = tab.mvs, conditions = conditions,
                            q.min = q.min, q.norm = q.norm, eps = eps,
                            distribution = distribution, param1 = param1,
                            param2 = param2, R.q.min=R.q.min)
  }

  return(tab.mvs.imp);
}



