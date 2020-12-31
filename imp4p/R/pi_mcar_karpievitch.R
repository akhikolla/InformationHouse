
pi.mcar.karpievitch=function(tab,conditions){
  nb_cond=length(levels(conditions))
  nb_rep=rep(0,nb_cond)
  for (i in 1:nb_cond){
    nb_rep[i]=sum((conditions==levels(conditions)[i]));
  }
  
  p.NA=matrix(0,length(tab[,1]),nb_cond)
  pi=rep(0,nb_cond)
  moy=p.NA
  k=1
  for (i in 1:nb_cond){
    #calcul de la proportion de NA dans chaque condition
    nb.NA<-apply(tab[,k:(k+nb_rep[i]-1)],1,function(x) length(which(is.na(x))==TRUE));
    p.NA[,i]=nb.NA/(nb_rep[i]);   
    #calcul du log2 moyen
    moy[,i]=apply(tab[,k:(k+nb_rep[i]-1)],1,mean,na.rm=T);
    #r?gression par splines cubiques de degr? 5
    reg=smooth.spline(p.NA[!is.na(moy[,i]),i]~moy[!is.na(moy[,i]),i],df=5)$y;
    pi[i]=reg[length(reg)]/(sum(nb.NA)/nb_rep[i]);
    k=k+nb_rep[i];
  }
  
  pi[pi<0]=0;
  pi[pi>1]=1;
  
  return(list(pi.mcar=pi,prop.na=p.NA,moy=moy)) 
}
 