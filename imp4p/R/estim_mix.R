
####################################
#
# Estimating the mixture model
#
####################################


estim.mix=function(tab, tab.imp, conditions, x.step.mod=150, x.step.pi=150, nb.rei=200){

  #require(Iso);
  min.x=min(tab.imp,na.rm=TRUE);
  max.x=max(tab.imp,na.rm=TRUE);

  x.min=min.x-(max.x-min.x)/2;
  x.max=max.x+(max.x-min.x)/2;

  #Function to minimize under Weibull assumption
  fr <- function(pi,x,pi_est,w,Ft) {
    K <- pi[1];
    a <- pi[2];
    d <- pi[3];
    f=K+(1-K)*exp(-a*(x)^d)/(1-Ft);
    return(sum(w*(pi_est-f)^2))
  }

  #Function to estimate the asymptotic variance of pi_mcar
  asy_v <- function(a,u,v) {
    delta=(a-1)*u/(a*u-1)
    ll=(a*(v-1)*u-v*u+1)/(1-a)
    g=(1-delta*v)*delta*v*(((1/u)+v/ll)^2)/((1-a*u)^2)
    h=((1/a)-1)/(v*(1-u))+(1-delta*v)*delta*v*((1/v)+u/ll)^2
    k=(1-delta*v)*((1-a)*ll^2)
    av=(1-a)*h/(a*(g*h-k^2))
    return(av)
  }

  #Number of conditions
  nb_cond=length(levels(conditions));

  #Initialisation number of replicates by condition
  nb_rep=rep(0,nb_cond);

  #Initialisation interval on which pi is estimated
  abs=matrix(0,x.step.pi,length(tab[1,]));

  #Initialisation number of missing values on each row in the condition
  dat_nbNA=matrix(0,length(tab[,1]),nb_cond);

  pi_abs=rep(0,nb_cond);

  #Proportion of missing values in each sample
  pi_na=rep(0,length(tab[1,]));

  #Proportion of MCAR values in each sample
  pi_mcar=rep(0,length(tab[1,]));

  MinRes=matrix(0,4,length(tab[1,]));

  sta=rep(0,length(tab[1,]));

  #Number of missing values on each row in the condition
  dat_nbNA=matrix(0,length(tab[,1]),nb_cond);

  PI_INIT=matrix(0,x.step.pi,length(tab[1,]));
  TREND_INIT=matrix(0,x.step.pi,length(tab[1,]));
  VAR_PI_INIT=matrix(0,x.step.pi,length(tab[1,]));

  ABSC=matrix(0,x.step.mod,length(tab[1,]));
  FTOT=matrix(0,x.step.mod,length(tab[1,]));
  FNA=matrix(0,x.step.mod,length(tab[1,]));
  FOBS=matrix(0,x.step.mod,length(tab[1,]));

  k=1;
  j=1;
  for (i in 1:nb_cond){
    #Number of replicates in each condition
    nb_rep[i]=sum((conditions==levels(conditions)[i]));
    #Number of missing values on each row in the condition
    dat_nbNA[,i]=apply(tab[,(k:(k+nb_rep[i]-1))],1, function(x) sum(is.na(x)));
    #Percentage of rows without observed values in the condition
    pi_abs[i]=sum(dat_nbNA[,i]==nb_rep[i])/sum(dat_nbNA[,i]>0);
    #Deletion of rows without data in the condition
    liste_retire=which(dat_nbNA[,i]==nb_rep[i]);
    liste_garde=which(dat_nbNA[,i]!=nb_rep[i]);
    #Data after deletion of rows without data in the condition
    tab2=as.matrix(tab[liste_garde,])
    ############
    #Roughly imputed data after deletion
    tab_imp2=as.matrix(tab.imp[liste_garde,]);

    ############

    while (j<(k+nb_rep[i])){
      #Percentage of missing values in sample j
      nna=sum(is.na(tab2[,j]));
      #Percentage of missing values in sample j (after deletion of empty rows in the condition)
      pi_na[j]=sum(is.na(tab2[,j]))/length(tab2[,j]);

      ############
      #Rough estimation of the distribution of missing values in sample j
      v_na=tab_imp2[is.na(tab2[,j]),j];
      F_na=ecdf(v_na);

      ############
      #Distribution of observed values in sample j
      F_obs=ecdf(tab2[which(!is.na(tab2[,j])),j]);

      ############
      #Determination of the interval on which is estimated pi^MLE(x)
      #F_na(x)>0 et F_obs(x)>0
      xmin=max(min(v_na),min(tab2[which(!is.na(tab2[,j])),j]));
      #F_na(x)<1 et F_obs(x)<1
      xmax=min(-sort(-v_na,na.last=T)[2],-sort(-tab2[which(!is.na(tab2[,j])),j],na.last=T)[2]);
      #initial interval
      absi=seq(xmin,xmax,length.out=x.step.pi)
      pi_init=NULL;
      for (x in absi){
        Ftotv=pi_na[j]*F_na(x)+(1-pi_na[j])*F_obs(x)
        pi_init_x=(1-F_na(x))/(1-Ftotv);
        pi_init=c(pi_init,pi_init_x);
      }
      #Determination of Mj
      kMj=1
      Mj=1
      while(pi_init[kMj]>=mean(pi_init)){Mj=absi[kMj];kMj=kMj+1;}
      abs[,j]=seq(Mj,xmax,length.out=x.step.pi);

      ############
      #Estimation of pi^MLE(x)
      varasy=NULL;
      pi_init=NULL;
      F_tot=NULL
      for (x in abs[,j]){

        ############
        #Initial estimator of the proportion
        one_minus_Ftot=pi_na[j]*(1-F_na(x))+(1-pi_na[j])*(1-F_obs(x))
        pi_init_x=(1-F_na(x))/(one_minus_Ftot);

        ############
        #Estimation of the asymptotic variance of the estimator
        pp=(1-F_obs(x))
        aa=pi_na[j]
        uu=pi_init_x

        if (uu>=1){uu=1-1e-8;}
        if (aa>=1){aa=1-1e-8;}
        if (pp>=1){pp=1-1e-8;}
        if (uu<=0){uu=1e-8;}
        if (aa<=0){aa=1e-8;}
        if (pp<=0){pp=1e-8;}

        pi_init=c(pi_init,uu)
        varasy=c(varasy,asy_v(aa,uu,pp))
        F_tot=c(F_tot,1-one_minus_Ftot)

      }
      varasy=varasy/length(tab2[which(!is.na(tab2[,j])),j])

      PI_INIT[,j]=pi_init;
      VAR_PI_INIT[,j]=varasy;

      #Minimization step
      #The minimization is performed nb.rei times to try to find a global minimum
      pim=matrix(0,nb.rei,4)
      for (nbit in 1:nb.rei){
        #First initialisation of the 3 parameters
        #pi_mcar is initialized between the min of pi_init - 0.1 and this same value + 0.1
        #alpha is initialized between 0 and 1
        #d is initialized between 0.5 and 15
        init=c(runif(1,max(c(0,min(pi_init)-0.1)),min(c(min(pi_init)+0.1,1))),runif(1,0,1),runif(1,0.5,15))
        nbtest=1;
        #Find a working starting point for the minimization
        while ((!inherits(try(optim(init, fr, gr=NULL, x=(abs[,j]-xmin) ,
                                    pi_est=pi_init, lower=c(max(c(0,min(pi_init)-0.1)),0,0.5),
                                    upper=c(min(c(min(pi_init)+0.1,1)),1,15), method="L-BFGS-B",
                                    w=(1/varasy)/(sum(1/varasy)), Ft=F_tot), TRUE), "try-error")==FALSE)&(nbtest<nb.rei)){
          init=c(runif(1,max(c(0,min(pi_init)-0.1)),min(c(min(pi_init)+0.1,1))),runif(1,0,1),runif(1,0.5,15));
          nbtest=nbtest+1;
        }
        #Minimization
        re=optim(init, fr, gr=NULL, x=(abs[,j]-xmin) ,
                 pi_est=pi_init, lower=c(0+1e-08,0,0.5),
                 upper=c(1-1e-08,1,15), method="L-BFGS-B",
                 w=(1/varasy)/(sum(1/varasy)), Ft=F_tot);
        pim[nbit,1]=re$par[1]
        pim[nbit,2]=re$par[2]
        pim[nbit,3]=re$par[3]
        pim[nbit,4]=re$value[1]
      }
      #Final minimizations from the three best starting points and adding a small variation
      re=optim(pim[order(pim[,4]),][1,1:3]+1e-3, fr, gr=NULL, x=(abs[,j]-xmin) ,
               pi_est=pi_init, lower=c(0+1e-08,0,0.5),
               upper=c(1-1e-08,1,10), method="L-BFGS-B",
               w=(1/varasy)/(sum(1/varasy)), Ft=F_tot);
      pim=rbind(pim,c(re$par[1],re$par[2],re$par[3],re$value));
      re=optim(pim[order(pim[,4]),][2,1:3]+1e-3, fr, gr=NULL, x=(abs[,j]-xmin) ,
               pi_est=pi_init, lower=c(0+1e-08,0,0.5),
               upper=c(1-1e-08,1,10), method="L-BFGS-B",
               w=(1/varasy)/(sum(1/varasy)), Ft=F_tot);
      pim=rbind(pim,c(re$par[1],re$par[2],re$par[3],re$value));
      re=optim(pim[order(pim[,4]),][3,1:3]+1e-3, fr, gr=NULL, x=(abs[,j]-xmin) ,
               pi_est=pi_init, lower=c(0+1e-08,0,0.5),
               upper=c(1-1e-08,1,10), method="L-BFGS-B",
               w=(1/varasy)/(sum(1/varasy)), Ft=F_tot);
      pim=rbind(pim,c(re$par[1],re$par[2],re$par[3],re$value));

      #Final estimation of the proportion of MCAR values: the best of all the previous minimizations
      pi_mcar[j]=pim[which.min(pim[,4]),1];

      MinRes[1:4,j]=pim[which.min(pim[,4]),1:4];
      ############
      #Find the normal distribution of complete values with
      #Regression between sufficiently high quantiles of observed values and normal quantiles
      trend=pim[which.min(pim[,4]),1]+((1-pim[which.min(pim[,4]),1])/(1-F_tot))*exp(-pim[which.min(pim[,4]),2]*(abs[,j]-xmin)^(pim[which.min(pim[,4]),3]));
      TREND_INIT[,j]=trend;
      pq=rep(0,length(trend))
      for (ii in 1:length(trend)){
        pq[ii]=pnorm(q=pi_init[ii],mean=trend[ii],sd=sqrt(varasy)[i])
      }

      eta=min(abs[which(pq>0.05),j])
      sta[j]=F_obs(eta);

      gamma=pi_na[j]*(1-pi_mcar[j])/(1-pi_mcar[j]*pi_na[j]);
      interv=gamma+(1-gamma)*seq(sta[j],1,length.out=100);
      q.normal = qnorm(interv, mean = 0, sd = 1);
      q.curr.sample = quantile(tab2[,j], probs = seq(sta[j],1,length.out=100), na.rm = T);
      q.normal = q.normal[1:(length(q.normal)-1)]
      q.curr.sample = q.curr.sample[1:(length(q.curr.sample)-1)]
      temp.QR = lm(q.curr.sample ~ q.normal);
      #mean and variance of the estimated normal density of complete values
      m = temp.QR$coefficients[1];
      v = (as.numeric(temp.QR$coefficients[2]))^2;

      ############
      #Final estimation of the mixture model on the interval precised in input
      ABSC[,j]=seq(x.min,x.max,length.out=x.step.mod);

      F_TOT=NULL
      F_NA=NULL
      F_OBS=NULL
      for (x in ABSC[,j]){
          F_TOT=c(F_TOT,pnorm(x,mean=m,sd=sqrt(v)));
          F_NA=c(F_NA,F_na(x));
          F_OBS=c(F_OBS,F_obs(x));
      }

      FTOT[,j]=F_TOT;
      FNA[,j]=F_NA;
      FOBS[,j]=F_OBS;

      j=j+1;
    }

    k=k+nb_rep[i];
  }

  return(list(abs.pi=abs,pi.init=PI_INIT,var.pi.init=VAR_PI_INIT,trend.pi.init=TREND_INIT,
              abs.mod=ABSC[,1],pi.na=pi_na,F.na=FNA,F.tot=FTOT,F.obs=FOBS,pi.mcar=pi_mcar,MinRes=MinRes));
}
