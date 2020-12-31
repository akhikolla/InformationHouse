
impute.mi=function(tab, conditions, repbio=NULL, reptech=NULL, nb.iter=3, nknn=15, selec=1000, siz=900,
                   weight=1, ind.comp=1, progress.bar=TRUE, x.step.mod=300, x.step.pi=300, nb.rei=100,
                   q=0.95, methodMCAR="mle",
                   ncp.max=5,maxiter = 10, ntree = 100, variablewise = FALSE,
                   decreasing = FALSE, verbose = FALSE,
                   mtry = floor(sqrt(ncol(tab))), replace = TRUE,
                   classwt = NULL, cutoff = NULL, strata = NULL,
                   sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                   xtrue = NA, parallelize = c('no', 'variables', 'forests'),
                   methodMNAR="igcda",q.min = 0.025, q.norm = 3, eps = 0,
                   distribution = "unif", param1 = 3, param2 = 1, R.q.min=1){

  if (is.null(repbio)){repbio=as.factor(1:length(conditions));}
  if (is.null(reptech)){reptech=as.factor(1:length(conditions));}

  if (progress.bar==TRUE){cat(paste("\n 1/ Initial imputation under the MCAR assumption... \n  "));}

  #Imputation of missing values with the slsa algorithm
  if (methodMCAR == "mle") {
    tab.imp = impute.mle(tab=tab, conditions = conditions);
  }else {
    if (methodMCAR == "rf") {
      tab.imp = impute.RF(tab=tab, conditions = conditions,
                          maxiter = maxiter, ntree = ntree, variablewise = variablewise,
                          decreasing = decreasing, verbose = verbose,
                          mtry = mtry, replace = replace,
                          classwt = classwt, cutoff = cutoff, strata = strata,
                          sampsize = sampsize, nodesize = nodesize, maxnodes = maxnodes,
                          xtrue = xtrue, parallelize = parallelize);
    }else {
      if (methodMCAR == "pca") {
        tab.imp = impute.PCA(tab=tab, conditions = conditions,
                             ncp.max=ncp.max);
      }else {
        tab.imp = impute.slsa(tab=tab, conditions = conditions, repbio = repbio,
                              reptech = reptech, nknn = nknn, weight = weight,
                              selec = selec, progress.bar = FALSE, ind.comp = ind.comp);
      }
    }
  }

  if (progress.bar==TRUE){cat(paste("\n\n 2/ Estimation of the mixture model in each sample... \n  "));}

  #Estimation of the mixture model
  res=estim.mix(tab=tab, tab.imp=tab.imp, conditions=conditions, x.step.mod=x.step.mod, x.step.pi=x.step.pi, nb.rei=nb.rei);

  if (progress.bar==TRUE){cat(paste("\n 3/ Estimation of the probabilities each missing value is MCAR... \n  "));}

  #Computing probabilities to be MCAR
  born=estim.bound(tab=tab,conditions=conditions,q=q);
  proba=prob.mcar.tab(born$tab.upper,res);

  if (progress.bar==TRUE){cat(paste("\n 4/ Multiple imputation strategy with mi.mix... \n  "));}

  #Multiple imputation strategy
  data.mi=mi.mix(tab=tab, tab.imp=tab.imp, prob.MCAR=proba, conditions=conditions, repbio=repbio, reptech=reptech,
                 nb.iter=nb.iter, nknn=nknn, weight=weight, selec=selec, siz=siz, ind.comp=ind.comp,
                 methodMCAR=methodMCAR, q=q, progress.bar=progress.bar, ncp.max=ncp.max,
                 maxiter = maxiter, ntree = ntree, variablewise = variablewise,
                 decreasing = decreasing, verbose = verbose,
                 mtry = mtry, replace = replace,
                 classwt = classwt, cutoff = cutoff, strata = strata,
                 sampsize = sampsize, nodesize = nodesize, maxnodes = maxnodes,
                 xtrue = xtrue, parallelize = parallelize,
                 methodMNAR=methodMNAR,q.min = q.min, q.norm =q.norm, eps = eps,
                 distribution = distribution, param1 = param1, param2 = param2, R.q.min=R.q.min);

  if (progress.bar==TRUE){cat(paste("\n\n 5/ Imputation of rows with only missing values in a condition with impute.pa... \n  "));}

  data.final=impute.pa(tab=data.mi, conditions=conditions, q.min=q.min, q.norm=q.norm, eps=eps);

  return(data.final$tab.imp)
}







