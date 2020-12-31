
#impute missing values of each column of tab
#with a small value.
#The small value is randomly draws from a uniform distribution
#between the min of observed values in the column
#and the q% quantile of these observed values.

impute.pa=function (tab, conditions, q.min = 0.025, q.norm = 3, eps = 0, distribution = "unif", param1 = 3, param2 = 1, R.q.min=1){

  if (length(colnames(tab))==0){
    colnames(tab)=seq(1,ncol(tab),by=1)
  }

  tab_imp = as.matrix(tab)

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

  qu = apply(tab_imp, 2, quantile, na.rm = TRUE, q.min)

  nb_cond = length(levels(conditions))

  nb_rep = rep(0, nb_cond)

  k = 1

  j = 1

  param=NULL

  for (i in 1:nb_cond) {

    nb_rep[i] = sum((conditions == levels(conditions)[i]))

    sde = apply(tab_imp[, (k:(k + nb_rep[i] - 1))], 1, sd, na.rm = TRUE)

    while (j < (k + nb_rep[i])) {

      if (distribution == "unif") {

        param=rbind(param,data.frame(n = sum(is.na(tab_imp[,j])), min = qu[j] - eps - q.norm * median(sde,na.rm = TRUE), max = qu[j] - eps))

        rownames(param)[j]=colnames(tab)[j]

        tab_imp[which(is.na(tab_imp[, j])), j] = runif(n = sum(is.na(tab_imp[,j])), min = qu[j] - eps - q.norm * median(sde,na.rm = TRUE), max = qu[j] - eps)

      }

      else if (distribution == "beta") {

        param=rbind(param,data.frame(n = sum(is.na(tab_imp[,j])), min = qu[j] - eps - q.norm * median(sde,na.rm = TRUE), max = qu[j] - eps, param1 = param1, param2 = param2))

        rownames(param)[j]=colnames(tab)[j]

        tab_imp[which(is.na(tab_imp[, j])), j] = translatedRandomBeta(n = sum(is.na(tab_imp[,j])), min = qu[j] - eps - q.norm * median(sde,na.rm = TRUE), max = qu[j] - eps, param1 = param1, param2 = param2)

      }

      else if (distribution == "dirac") {

        param=rbind(param,data.frame(imputed.value=R.q.min*quantile(tab[,j], probs=q.min, na.rm=T)))

        rownames(param)[j]=colnames(tab)[j]

        tab_imp[which(is.na(tab_imp[, j])), j] = R.q.min*quantile(tab[,j], probs=q.min, na.rm=T)

      }

      j = j + 1

    }

    k = k + nb_rep[i]

  }

  tab_imp[,index]=tab_imp

  return(list(tab.imp=tab_imp,para=param))

}

