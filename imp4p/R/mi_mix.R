
####################################
#
#Function to perform multiple imputation
#
####################################

mi.mix=function(tab, tab.imp, prob.MCAR, conditions, repbio = NULL, reptech = NULL, nb.iter = 3,
                nknn = 15, weight = 1, selec = "all", siz = 500, ind.comp = 1, methodMCAR = "mle",
                q = 0.95, progress.bar = TRUE, details= FALSE, ncp.max=5,
                maxiter = 10, ntree = 100, variablewise = FALSE,
                decreasing = FALSE, verbose = FALSE,
                mtry = floor(sqrt(ncol(tab))), replace = TRUE,
                classwt = NULL, cutoff = NULL, strata = NULL,
                sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                xtrue = NA, parallelize = c('no', 'variables', 'forests'),
                methodMNAR="igcda",q.min = 0.025, q.norm = 3, eps = 0,
                distribution = "unif", param1 = 3, param2 = 1, R.q.min=1){

  if (is.null(repbio)) {repbio = as.factor(1:length(conditions));}
  if (is.null(reptech)) {reptech = as.factor(1:length(conditions));}
  if (selec>=nrow(tab)){selec=nrow(tab)-1;}
  if (selec == "all") {
    if (siz > nrow(tab)) {siz = nrow(tab) - 1;}
  }
  else {
    if (siz > selec) {siz = selec - 1;}
  }

  data_imp = array(NA, dim = c(nrow(tab), ncol(tab), nb.iter))
  l.NA = matrix(0, nrow(tab), ncol(tab))
  if (progress.bar == TRUE) {cat(paste("\n Iterations: \n"));}
  iter = 1
  while ((iter <= nb.iter)) {
    if (progress.bar == TRUE) {
      cat(paste("\n", iter, "/", nb.iter, " - "))
    }
    tab.mvs = as.matrix(tab)

    #matrix of missing values
    l.NA[which(is.na(tab.mvs))] = 1

    #matrix of potential MCAR values
    l.MCAR = matrix(mapply(prob.MCAR, FUN = rbinom, size = 1, n = 1), nrow(prob.MCAR), ncol(prob.MCAR))

    #Missing and MCAR
    l.MCAR = l.NA * l.MCAR

    #Replace MCAR values by tab.imp values
    tab.mvs[which(l.MCAR == 1)] = tab.imp[which(l.MCAR == 1)]

    data_imp[, , iter] = tab.mvs

    #Replace remaining missing values using MNAR-devoted algorithms
    if (methodMNAR == "igcda"){
        tab.mvs.imp = impute.igcda(tab = tab.mvs, tab.imp = tab.imp,
                               conditions = conditions, q = q)
    }else{
      tab.mvs.imp = impute.pa(tab = tab.mvs, conditions = conditions,
                              q.min = q.min, q.norm = q.norm, eps = eps,
                              distribution = distribution, param1 = param1,
                              param2 = param2, R.q.min=R.q.min)
      tab.mvs.imp = tab.mvs.imp$tab.imp
    }

    if (progress.bar == TRUE) {
      cat(paste("Imputation MNAR OK - \n"))
    }

    rna = rowSums(l.NA)
    rna = which(rna > 0)

    if (progress.bar == TRUE) {
      cat(paste("Imputation MCAR in progress - \n"))
    }

    for (i in 1:length(rna)) {
      tab.mod = tab.mvs.imp
      tab.mod.imp = tab.mod
      tab.mod[rna[i], ] = tab.mvs[rna[i], ]

      lab = (1:nrow(tab.mod))[-rna[i]]
      sel = selec
      if (selec == "all") {sel = nrow(tab.mod) - 1;}
      list.select = sample(lab, size = max(sel, min(siz,nrow(tab.mod) - 1)), replace = FALSE)
      list.select = c(list.select, rna[i])

      if (methodMCAR == "mle") {
        tab.mod.imp[list.select, ] = impute.mle(tab.mod[list.select,], conditions = conditions);
      }else{
        if (methodMCAR == "rf") {
          tab.mod.imp[list.select, ] = impute.RF(tab.mod[list.select,], conditions = conditions,
                                                 maxiter = maxiter, ntree = ntree, variablewise = variablewise,
                                                 decreasing = decreasing, verbose = verbose,
                                                 mtry = mtry, replace = replace,
                                                 classwt = classwt, cutoff = cutoff, strata = strata,
                                                 sampsize = sampsize, nodesize = nodesize, maxnodes = maxnodes,
                                                 xtrue = xtrue, parallelize = parallelize);
        }else {
          if (methodMCAR == "pca") {
            tab.mod.imp[list.select, ] = impute.PCA(tab.mod[list.select,], conditions = conditions,
                                                    ncp.max=ncp.max);
          }else {
            tab.mod.imp[list.select, ] = impute.slsa(tab.mod[list.select,], conditions = conditions, repbio = repbio,
                                                 reptech = reptech, nknn = nknn, weight = weight,
                                                 selec = selec, progress.bar = FALSE, ind.comp = ind.comp);
          }
        }
      }

      data_imp[rna[i], , iter] = tab.mod.imp[list.select,][length(list.select), ];
      if (progress.bar == TRUE) {
        if (i%/%floor(0.01 * length(rna)) == i/floor(0.01 *length(rna))) {
          cat(c((i%/%(0.01 * length(rna)))), "% - ")
        }
      }
    }
    iter = iter + 1
  }
  if (details == TRUE){
      data_fin = list(
      imputed.matrix=apply(data_imp, 1:2, mean),
      sd.imputed.matrix=apply(data_imp, 1:2, sd),
      all.imputed.matrices=data_imp);
  }else{
      data_fin = apply(data_imp, 1:2, mean)
  }
  return(data_fin)
}
