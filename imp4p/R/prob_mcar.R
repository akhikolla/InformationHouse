
#Compute a vector of probabilities to be MCAR

prob.mcar=function(b.u,absc,pi.na,pi.mcar,F.tot,F.obs){
  prob=rep(0,length(b.u));
  for (i in 1:length(b.u)){

      if (!is.na(b.u[i])){
          if (!is.nan(b.u[i])){
              if (b.u[i]!=Inf){

                i.bsup=which.min(abs(absc-b.u[i]));
                if ((F.tot[i.bsup]!=0)){
                  prob[i]=pi.mcar*pi.na/(1-(1-pi.na)*(F.obs[i.bsup])/(F.tot[i.bsup]));
                }else{
                  prob[i]=pi.mcar*pi.na;
                }

            }
        }
    }

  }

  prob[prob>1]=1;
  prob[prob<0]=0;

  return(prob)
}

