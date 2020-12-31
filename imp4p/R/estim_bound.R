
#Determine lower and upper bounds for missing values

estim.bound=function(tab,conditions,q=0.95){

  nb_cond=length(levels(conditions));
  nb_rep=rep(0,nb_cond);
  k=1;
  lb=NULL;
  ub=NULL;
  for (n in 1:nb_cond){
    nb_rep[n]=sum((conditions==levels(conditions)[n]));
    nb.NA=apply(tab[,(k:(k+nb_rep[n]-1))],1,function(x){sum(is.na(x));});

    mat=apply(tab[(nb.NA!=nb_rep[n]),(k:(k+nb_rep[n]-1))],1,max,na.rm=TRUE);
    mit=apply(tab[(nb.NA!=nb_rep[n]),(k:(k+nb_rep[n]-1))],1,min,na.rm=TRUE);
    rmatmit=mat-mit;
    rangem=quantile(rmatmit,q);
    bi=mat-rangem;
    bs=mat;
    binf=rep(0,length(tab[,1]));
    bsup=rep(0,length(tab[,1]));
    binf[(nb.NA!=nb_rep[n])]=bi;
    bsup[(nb.NA!=nb_rep[n])]=bs;

    lb=cbind(lb,matrix(rep(binf,ncol(tab[,(k:(k+nb_rep[n]-1))])),nrow(tab[,(k:(k+nb_rep[n]-1))]),ncol(tab[,(k:(k+nb_rep[n]-1))])));
    ub=cbind(ub,matrix(rep(bsup,ncol(tab[,(k:(k+nb_rep[n]-1))])),nrow(tab[,(k:(k+nb_rep[n]-1))]),ncol(tab[,(k:(k+nb_rep[n]-1))])));

    k=k+nb_rep[n];
  }
  lb[lb==0]=NA;
  ub[ub==0]=NA;

  return(list(tab.lower=lb,tab.upper=ub));
}
