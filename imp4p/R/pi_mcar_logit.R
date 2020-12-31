
#Mod?lisation du processus des valeurs manquantes ? l'aide d'un mod?le logit
pi.mcar.logit=function(tab,conditions){

nb_cond=length(levels(conditions))
nb_rep=rep(0,nb_cond)
for (i in 1:nb_cond){
  nb_rep[i]=sum((conditions==levels(conditions)[i]));
}

j=1
tab2=tab
moyy=matrix(0,length(tab[,1]),nb_cond)
pi_mcar=rep(0,length(tab[1,]))
pi_mcar2=pi_mcar
pi_na=pi_mcar
a=rep(0,length(tab[1,]))
b=a
for (i in 1:nb_cond){
  k=j;
  
  moyy [,i]=apply(tab[,k:(k+nb_rep[i]-1)],1,mean,na.rm=T);
  
  while (j<(k+nb_rep[i])){
    nna=sum(is.na(tab[,j]))
    pi_na[j]=nna/length(tab[,j])
    tab2[,j]=as.numeric(is.na(tab[,j]))
    mylogit <- glm(tab2[,j] ~ moyy[,i], family = "binomial")
    a[j]=mylogit$coefficients[1]
    b[j]=mylogit$coefficients[2]
    #pi_mcar[j]=(exp(a[j])/(1+exp(a[j])))/pi_na[j]
    pi_mcar[j]=(exp(a[j]+b[j]*max(moyy[,i],na.rm=T))/(1+exp(a[j]+b[j]*max(moyy[,i],na.rm=T))))/pi_na[j]
    j=j+1;
  }
}

return(list(pi.mcar=pi_mcar,coef1=a,coef2=b))
}
  
  
  
  
  