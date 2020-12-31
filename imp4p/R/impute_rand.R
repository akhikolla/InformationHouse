
impute.rand=function(tab,conditions){

  tab_imp=as.matrix(tab);

  nb_cond=length(levels(conditions));
  nb_rep=rep(0,nb_cond);
  k=1;

  for (n in 1:nb_cond){
  
    nb_rep[n]=sum((conditions==levels(conditions)[n]));

    md=suppressWarnings(apply(tab[,(k:(k+nb_rep[n]-1))],1,mean,na.rm=T));
    asd=suppressWarnings(apply(tab[,(k:(k+nb_rep[n]-1))],1,sd,na.rm=T));
    masd=quantile(asd,0.25,na.rm=T);

    for (i in 1:length(md)){
      tab_imp[i,(k:(k+nb_rep[n]-1))][which(is.na(tab[i,(k:(k+nb_rep[n]-1))]))]=md[i]+rnorm(sum(is.na(tab[i,(k:(k+nb_rep[n]-1))])),0,masd);
    }

    k=k+nb_rep[n];
  }
  
  tab_imp[is.nan(tab_imp)]=NA;

  return(tab_imp);
}
  

 


