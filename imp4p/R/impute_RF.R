
impute.RF=function (tab, conditions,
                    maxiter = 10, ntree = 100, variablewise = FALSE,
                    decreasing = FALSE, verbose = FALSE,
                    mtry = floor(sqrt(ncol(tab))), replace = TRUE,
                    classwt = NULL, cutoff = NULL, strata = NULL,
                    sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                    xtrue = NA, parallelize = c('no', 'variables', 'forests')) {

  tab_imp=as.matrix(tab);

  new_tab=NULL
  new_conditions=NULL
  index=NULL
  for (j in 1:length(levels(conditions))){
    index=c(index,which(conditions==levels(conditions)[j]))
    new_tab=cbind(new_tab,tab_imp[,which(conditions==levels(conditions)[j])])
    new_conditions=c(new_conditions,conditions[which(conditions==levels(conditions)[j])])
  }

  tab_imp=new_tab
  conditions=new_conditions
  conditions=factor(as.character(conditions),levels=as.character(unique(conditions)));

  nb_cond=length(levels(conditions));
  nb_rep=rep(0,nb_cond);
  k=1;

  for (n in 1:nb_cond){
    nb_rep[n]=sum((conditions==levels(conditions)[n]));
    xincomplete=as.matrix(tab[,(k:(k+nb_rep[n]-1))]);
    nbna=fast_apply_nb_na(xincomplete,1);
    if (sum(nbna)>0){
      xincomplete1=xincomplete[which(nbna!=nb_rep[n]),];
      nbna2=fast_apply_nb_na(xincomplete1,1);
      if (sum(nbna2)>0){
        rngseed(1234567);
        xcomplete1 <- missForest(xincomplete1,
                             maxiter = maxiter, ntree = ntree, variablewise = variablewise,
                             decreasing = decreasing, verbose = verbose,
                             mtry = mtry, replace = replace,
                             classwt = classwt, cutoff = cutoff, strata = strata,
                             sampsize = sampsize, nodesize = nodesize, maxnodes = maxnodes,
                             xtrue = xtrue, parallelize = parallelize)

        tab_imp[which(nbna!=nb_rep[n]),(k:(k+nb_rep[n]-1))]=xcomplete1$ximp;
      }
    }

    k=k+nb_rep[n];
  }

  tab_imp[,index]=tab_imp

  return(tab_imp)
}
