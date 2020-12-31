#Simulation of data sets
#with a hierarchical design

sim.data=function(nb.pept=15000,nb.miss=5000,pi.mcar=0.2,para=3,nb.cond=1,nb.repbio=3,nb.sample=3,m.c=25,sd.c=2,sd.rb=0.5,sd.r=0.2){

  tab=matrix(0,nb.pept,nb.sample*nb.cond*nb.repbio);
  pmnar=tab;
  moyenne=pmnar;
  moyenne.cond=matrix(0,nb.pept,nb.cond);
  moyenne.repbio=matrix(0,nb.pept,nb.repbio);

  nbMCAR=floor(nb.miss*pi.mcar);
  nbMNAR=floor(nb.miss*(1-pi.mcar));

  #Generating the complete values using a normal distribution
  xmin=rep(0,ncol(tab));
  xmax=xmin;
  for (i in 0:(nb.cond-1)){
    moyenne.cond[1:(nb.pept-nbMNAR),i+1]=rnorm(nb.pept-nbMNAR,mean=m.c,sd=sd.c);
    for (rb in 0:(nb.repbio-1)){
      moyenne.repbio[1:(nb.pept-nbMNAR),rb+1]=rnorm(nb.pept-nbMNAR,mean=0,sd=sd.rb);
      for (j in (nb.sample*(rb+i*nb.repbio)+1):(nb.sample*(rb+i*nb.repbio)+nb.sample)){
        tab[1:(nb.pept-nbMNAR),j]=rnorm(nb.pept-nbMNAR,mean=moyenne.cond[1:(nb.pept-nbMNAR),i+1]+moyenne.repbio[1:(nb.pept-nbMNAR),rb+1],sd=sd.r);
        xmin[j]=min(qnorm(0.01,mean=moyenne.cond[1:(nb.pept-nbMNAR),i+1]+moyenne.repbio[1:(nb.pept-nbMNAR),rb+1],sd=sd.r));
        xmax[j]=max(qnorm(0.99,mean=moyenne.cond[1:(nb.pept-nbMNAR),i+1]+moyenne.repbio[1:(nb.pept-nbMNAR),rb+1],sd=sd.r));
      }
    }
  }

  #Creating enough MNAR values
  for (i in 0:(nb.cond-1)){
      #1) generate random probability to be mnar
      p_mnar=runif(nbMNAR,0,1);
      moyenne.cond[(nb.pept-nbMNAR+1):nrow(tab),i+1]=min(xmin)+(max(xmax)-min(xmin))*(1-p_mnar)/para;
      for (rb in 0:(nb.repbio-1)){
          moyenne.repbio[(nb.pept-nbMNAR+1):nrow(tab),rb+1]=rnorm(nbMNAR,mean=0,sd=sd.rb);
          for (j in (nb.sample*(rb+i*nb.repbio)+1):(nb.sample*(rb+i*nb.repbio)+nb.sample)){
            tab[(nb.pept-nbMNAR+1):nrow(tab),j]=rnorm(nbMNAR,mean=moyenne.cond[(nb.pept-nbMNAR+1):nrow(tab),i+1]+moyenne.repbio[(nb.pept-nbMNAR+1):nrow(tab),rb+1],sd=sd.r);
          }
      }
  }


  #Generating the probabilities to be MNAR
  for (j in 1:ncol(tab)){
    for (i in 1:nrow(tab)){
      if((1-para*(tab[i,j]-xmin[j])/(xmax[j]-xmin[j]))>=0){
        pmnar[i,j]=1-para*(tab[i,j]-xmin[j])/(xmax[j]-xmin[j]);
      }else{pmnar[i,j]=1e-08;}
    }
  }

  tab.comp=tab;

  #Uniform drawing in the initial distribution to generate MCAR values
  listeMCAR=matrix(0,nbMCAR,nb.sample*nb.cond*nb.repbio);

  if (pi.mcar!=0){
    for (i in 0:(nb.cond-1)){
      for (rb in 0:(nb.repbio-1)){
        for (j in (nb.sample*(rb+i*nb.repbio)):(nb.sample*(rb+i*nb.repbio)+nb.sample-1)){
          liste=sample(1:(nb.pept-nbMNAR),size=nbMCAR);
          tab[liste,j+1]=NA;
          listeMCAR[,j+1]=sort(liste);
        }
      }
    }
  }


  # #If we do not have enough generated MNAR values corresponding to the chosen distribution in the original simulated dataset,
  # #we add some potential MNAR values in the dataset
  # avera=matrix(0,nb.pept,nb.sample*nb.cond*nb.repbio);
  # for (i in 0:(nb.cond-1)){
  #   list_max_notsmallpmnar=which(apply(pmnar[,(nb.sample*(i*nb.repbio)+1):(nb.sample*((nb.repbio-1)+i*nb.repbio)+nb.sample)],1,max)!=1e-08)
  #   list_MCAR_values=as.numeric(levels(as.factor(as.vector(listeMCAR[,(nb.sample*(i*nb.repbio)):(nb.sample*((nb.repbio-1)+i*nb.repbio)+nb.sample-1)+1]))));
  #   if (length(list_max_notsmallpmnar)<nbMNAR){
  #       avera=moyenne.cond[,i+1]+moyenne.repbio[,rb+1];
  #       maxmnar=max(avera[list_max_notsmallpmnar]);
  #       minmnar=min(avera[list_max_notsmallpmnar]);
  #       #list of new potential MNARs
  #       #nbcolonne=length((nb.sample*(rb+i*nb.repbio)+1):(nb.sample*(rb+i*nb.repbio)+nb.sample))
  #       x=which((avera>maxmnar))
  #       list_potential_new_mnar=x[!x%in%list_MCAR_values];
  #       list_new_mnar=list_potential_new_mnar[floor(runif(nbMNAR-length(list_max_notsmallpmnar),1,length(list_potential_new_mnar)))];
  #       for (k in list_new_mnar){
  #           moyenne.cond[k,i+1]=runif(1,minmnar,maxmnar);
  #           for (rb in 0:(nb.repbio-1)){
  #               moyenne.repbio[k,rb+1]=rnorm(1,mean=0,sd=sd.rb);
  #               for (j in (nb.sample*(rb+i*nb.repbio)+1):(nb.sample*(rb+i*nb.repbio)+nb.sample)){
  #                   tab[k,j]=rnorm(1,mean=moyenne.cond[k,i+1]+moyenne.repbio[k,rb+1],sd=sd.r);
  #               }
  #           }
  #       }
  #
  #   }
  # }
  #
  # #Recomputing of the probability to be MNAR
  # pmnar=matrix(0,nb.pept,nb.sample*nb.cond*nb.repbio);
  #
  # for (j in 1:ncol(tab)){
  #     xmin=min(tab[,j])
  #     xmax=max(tab[,j])
  #     for (i in 1:nrow(tab)){
  #     if (!is.na(tab[i,j])){
  #       if((1-para*(tab[i,j]-xmin)/(xmax-xmin))>=0){
  #         pmnar[i,j]=1-para*(tab[i,j]-xmin)/(xmax-xmin);
  #       }else{pmnar[i,j]=1e-08;}
  #     }
  #   }
  # }



  # #Uniform drawing in the initial distribution to generate MCAR values
  # listeMCAR=matrix(0,nbMCAR,nb.sample*nb.cond*nb.repbio);
  #
  # if (pi.mcar!=0){
  #   for (i in 0:(nb.cond-1)){
  #     for (rb in 0:(nb.repbio-1)){
  #       for (j in (nb.sample*(rb+i*nb.repbio)):(nb.sample*(rb+i*nb.repbio)+nb.sample-1)){
  #         liste=sample(1:nb.pept,size=nbMCAR);
  #         tab[liste,j+1]=NA;
  #         listeMCAR[,j+1]=sort(liste);
  #       }
  #     }
  #   }
  # }

  #Drawing to find MNAR values
  if (pi.mcar!=1){
    for (i in 0:(nb.cond-1)){
      for (rb in 0:(nb.repbio-1)){
        for (j in (nb.sample*(rb+i*nb.repbio)):(nb.sample*(rb+i*nb.repbio)+nb.sample-1)){

          d=1:nb.pept
          sa=d[-c(which(is.na(tab[,j+1])))];
          if (length(sa)>=nbMNAR){
            liste=sample(sa,size=nbMNAR,prob=pmnar[-c(which(is.na(tab[,j+1]))),j+1]);
          }else{
            warning(paste("\n The proportion of generated MCAR values is superior to the chosen proportion in sample ",j+1,".\n"));
            liste=sa;
          }

          tab[liste,j+1]=NA;

        }

      }

    }
  }


  #Delete rows with only missing values for each condition
  lMCAR=vector("list",ncol(listeMCAR))
  nMCAR=rep(0,ncol(tab));
  nNA=rep(0,ncol(tab));
  for (i in 0:(nb.cond-1)){
    colo=(nb.sample*(i*nb.repbio)+1):(nb.sample*((nb.repbio-1)+i*nb.repbio)+nb.sample);
    nb_na=apply(tab[,colo],1,function(x) sum(is.na(x)));
    ldel=which(nb_na==length(colo))
    if(length(ldel)!=0){
      tab1=tab[-ldel,];
      tab.comp1=tab.comp[-ldel,];
    }else{
      tab1=tab;
      tab.comp1=tab.comp;
    }
    nNA=apply(tab1,2,function(x) sum(is.na(x)));
    for (j in colo){
      ll1=listeMCAR[!listeMCAR[,j]%in%ldel,j];
      ll=NULL
      for (i in 1:length(ll1)){
        ll=c(ll,ll1[i]-sum(ldel<ll1[i]))
      }
      lMCAR[[j]]=ll;
      nMCAR[j]=length(lMCAR[[j]]);
    }
  }

  return(list(dat.obs=tab1, dat.comp=tab.comp1, list.MCAR=lMCAR, nMCAR=nMCAR, nNA=nNA, conditions=gen.cond(nb.cond,nb.repbio*nb.sample), repbio=gen.cond(nb.cond*nb.repbio,nb.sample)))
}

