#include <Rcpp.h>

using namespace Rcpp;



// [[Rcpp::export]]


SEXP B_est1(NumericVector time, NumericVector status, 
  NumericVector stime, 
  NumericVector stime1, NumericVector G, NumericVector X, 
  double tau_bet, 
  NumericMatrix epstheta, NumericMatrix Edot, double pdim){


  int n=time.size(), indik=0;
	
  int k=stime.size(), k_tmp=0;
	
  double muhat, C, D,  Bdot_tmp, tmpA=0, tmp11,tmp2,tmp3,tmp4;
	
  NumericVector B(k), dB(k),status_t(n), status_v(k), dN(n), Gc(n), Gt(n), 
mudot(1), 
    se(k), del_se(k), corA(k), stime_tmp(k), indik_v(k);
	
  NumericVector C_dot(pdim), D_dot(pdim), H(n), H_dot(n), resid(n), tmpF(k), 
F_st(k),
    tmpv1(n),Edot_vec(n) ;
	
  NumericMatrix H_resid(k,n), d_eps_B(k,n), eps_B(k,n), eps_B_thet(k,n), eps_B_tmp(k,n);
  
  NumericMatrix dB_dot(k,pdim), B_dot(k,pdim);

  
  double tot_risk=0, del=0.01;
  
  tot_risk=n*stime1[0];
  
  for (int j=1;j<n;j++){
    
    tot_risk+= (n-j)*(stime1[j]-stime1[j-1]);
 	
  }
  
  muhat=0;
  mudot[0]=1;
  
  for (int i=0;i<n;i++){
 	  
    muhat+=G[i]/n;
   	
    status_t[i]=(time[i]>=stime[0]) ? 1:0;
   	
    dN[i]=(time[i]==stime[0]) ? 1:0;
   	
    status_v[0]+=(time[i]>=stime[0]) ? 1:0;
 	
  }
 	
  double A=0, tmp1=0, beta=0;
 	
  for (int i=0;i<n;i++){
    
    Gc[i]=G[i];
		
    dB[0]+=Gc[i]*dN[i];
    
    A+=Gc[i]*status_t[i]*X[i];
		
    tmp1+=Gc[i]*status_t[i]*X[i]*X[i];
  
  }
 	
  for (int i=0;i<n;i++){
    
    H[i]=Gc[i]*status_t[i]/A;
		
    H_dot[i]=Gc[i]*status_t[i]*X[i]/A-Gc[i]*status_t[i]*tmp1/(A*A);
  
  }
    		
  
  dB[0]=dB[0]/A;
 	
  B[0]=dB[0];
 	
  corA[0]=A;
 	
  beta=beta+sum(status_t)*dB[0];
 	
  tmpF[0]=1-sum(H_dot*dN);
  
  F_st[0]=tmpF[0];
 	
  
  for (int j1=0;j1<pdim;j1++){
    
    tmp11=0;
    tmp2=0;
    tmp3=0;
    tmp4=0;
    
    for (int i=0;i<n;i++){
     	
      tmp11+=Edot(i,j1)*dN[i];
     	
      tmp2+=Gc[i]*status_t[i]*X[i];
      
      tmp3+=Gc[i]*dN[i];
      
      tmp4+=Edot(i,j1)*status_t[i]*X[i];
    
    }  	
    
    dB_dot(0,j1)=-tmp11/tmp2+tmp3*tmp4/(tmp2*tmp2);
    
    B_dot(0,j1)=dB_dot(0,j1);
  
  }
 	
  for (int i=0;i<n;i++){
		
    resid[i]=dN[i]-X[i]*status_t[i]*dB[0];
		
    H_resid(0,i)=H[i]*resid[i];
  
  }

  
  for (int j=1;j<k;j++){
    
    for (int i=0;i<n;i++){
    	
      status_t[i]=(time[i]>=stime[j]) ? 1:0;
    	
      status_v[j]+=(time[i]>=stime[j]) ? 1:0;
    	
      dN[i]=(time[i]==stime[j]) ? 1:0;
    	
      Gt[i]=Gc[i]*exp(B[j-1]*X[i]);
   	
    }
   	
    double A=0;tmp1=0;
   	
    for (int i=0;i<n;i++){
  		
      dB[j]+=Gc[i]*exp(B[j-1]*X[i])*dN[i];
  		
      A+=Gc[i]*status_t[i]*exp(B[j-1]*X[i])*X[i];
  		
      tmpA=(A<0) ? -A:A;
  		
      indik_v[j]=(tmpA<del) ? 0:1;
  		
      indik=(stime[j]<tau_bet) ? 1:0;
  		
      k_tmp=(stime[j]<tau_bet) ? j:k_tmp;
  		
      tmp1+=Gc[i]*status_t[i]*exp(B[j-1]*X[i])*X[i]*X[i];
		
    }
 	  
    for (int i=0;i<n;i++){
		  
      H[i]=Gc[i]*status_t[i]*exp(B[j-1]*X[i])/A;
  		
      H_dot[i]=Gc[i]*status_t[i]*exp(B[j-1]*X[i])*X[i]/A-Gc[i]*status_t[i]*exp(B[j-1]*X[i])*tmp1/(A*A);

      tmpv1[i]=exp(B[j-1]*X[i]);
		
    }
    		
 	  
    dB[j]=indik_v[j]*dB[j]/A;
   	
    B[j]=B[j-1]+dB[j];
   	
    stime_tmp[j]=(tmpA<del) ? stime_tmp[j-1]:stime[j];
   	
    corA[j]=A;
   	
    beta=beta+indik*sum(status_t)*dB[j];

    
    for (int i=0;i<n;i++){
  		
      resid[i]=dN[i]-X[i]*status_t[i]*dB[j];
  		
      H_resid(j,i)=H[i]*resid[i];
    
    }

    
    tmpF[j]=1-sum(H_dot*dN);
    
    F_st[j]=F_st[j-1]*tmpF[j];
    
    C=sum( Gc*tmpv1*dN);
    
    D=sum( Gc*tmpv1*status_t*X);
    
    
    for (int j1=0;j1<pdim;j1++){
      
      for (int i=0;i<n;i++){ 
        Edot_vec[i]=Edot(i,j1);
      }
   	    
      C_dot[j1]=sum(  (B_dot(j-1,j1)*Gc*X-Edot_vec)*tmpv1*dN );
        
      D_dot[j1]=sum(  (B_dot(j-1,j1)*Gc*X-Edot_vec)*tmpv1*X*status_t );
        
      dB_dot(j,j1)=indik_v[j]*(C_dot[j1]*D-C*D_dot[j1])/(D*D);
        
      B_dot(j,j1)=B_dot(j-1,j1)+dB_dot(j,j1);
      
    }
    
  }
                                 
    
  for (int i=0;i<n;i++){
      
    eps_B(0,i)=F_st[0]*H_resid(0,i);
   	  
    eps_B_tmp(0,i)=F_st[0]*H_resid(0,i);
   	  
    Bdot_tmp=0;
   	  
    for (int j1=0;j1<pdim;j1++){ 
      Bdot_tmp=Bdot_tmp+B_dot(0,j1)*epstheta(i,j1);
    }
 	    
    eps_B_thet(0,i)=eps_B(0,i)+Bdot_tmp;
    
  }
    
  for (int j=1;j<k;j++){
 	    
    for (int i=0;i<n;i++){
 	      
      d_eps_B(j,i)=F_st[j-1]*H_resid(j,i);
   	    
      eps_B_tmp(j,i)=eps_B_tmp(j-1,i)+indik_v[j]*d_eps_B(j,i);
   	    
      eps_B(j,i)=eps_B_tmp(j,i)/F_st[j];
   	    
      Bdot_tmp=0;
   	    
      for (int j1=0;j1<pdim;j1++){ 
        Bdot_tmp=Bdot_tmp+B_dot(j,j1)*epstheta(i,j1);
      }
 	      
      eps_B_thet(j,i)=eps_B(j,i)+Bdot_tmp;
 	    
    }
    
  }
    
  for (int j=0;j<k;j++){
 	    
    for (int i=0;i<n;i++){
 	      
      se[j]+=eps_B_thet(j,i)*eps_B_thet(j,i);
 	  
    }
 	  
    se[j]=sqrt(se[j]);
  
  }
  
  for (int j=1;j<k;j++){	
    del_se[j]=(se[j]-se[j-1])/se[j-1];
  }
  
  tot_risk=n*stime1[0];
  
  for (int j=1;j<n;j++){
    
    indik=(stime1[j]<=stime[k_tmp]) ? 1:0;
 	  
    tot_risk+= (n-j)*(stime1[j]-stime1[j-1])*indik;
 	
  }
 	
  beta=beta/tot_risk;

  
  return List::create(
Named("stime") = stime,
 Named("B")= B, Named("se")= se,
    
    Named("del_se")= del_se,
 Named("beta")= beta,
 Named("Bdot")= B_dot,
    
    Named("dBdot")= dB_dot,
 Named("tmp11")= tmp11,
 Named("tmp2")= tmp2,
  
    Named("tmp3")= tmp3,
 Named("tmp4")= tmp4,
 Named("eps_B")= eps_B_thet,
    
    Named("corA")= corA, Named("at_risk")= status_v,
 Named("tot_risk")= tot_risk
);
    

}

