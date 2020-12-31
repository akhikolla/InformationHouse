
#Mod?lisation du processus des valeurs manquantes ? l'aide d'un mod?le logit
pi.mcar.probit=function(tab,conditions){
  
  nb_cond=length(levels(conditions))
  nb_rep=rep(0,nb_cond)
  for (i in 1:nb_cond){
    nb_rep[i]=sum((conditions==levels(conditions)[i]));
  }
  
  j=1
  tab2=tab
  moy=matrix(0,length(tab[,1]),nb_cond)
  pi_mcar=rep(0,length(tab[1,]))
  pi_mcar2=pi_mcar
  pi_na=pi_mcar
  a=rep(0,length(tab[1,]))
  b=a
  for (i in 1:nb_cond){
    k=j;
    
    moy[,i]=apply(tab[,k:(k+nb_rep[i]-1)],1,mean,na.rm=T);
    
    while (j<(k+nb_rep[i])){
      nna=sum(is.na(tab[,j]))
      pi_na[j]=nna/length(tab[,j])
      tab2[,j]=as.numeric(is.na(tab[,j]));
      myprobit <- glm(tab2[,j] ~ moy[,i], family=binomial(link="probit"));
      a[j]=myprobit$coefficients[1];
      b[j]=myprobit$coefficients[2];
      #pi_mcar[j]=pnorm(a[j],mean=0,sd=1)/pi_na[j];
      pi_mcar[j]=pnorm(a[j]+b[j]*max(moy[,i],na.rm=T),mean=0,sd=1)/pi_na[j];
      j=j+1;
    }
  }
  
  return(list(pi.mcar=pi_mcar,coef1=a,coef2=b))
}
  




