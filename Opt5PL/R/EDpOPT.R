EDpOPT<-function(LB,UB,P,EDp,grid=0.01,r=30,epsilon=.001,N_dose=FALSE,log_scale=TRUE)
{
      I5<-10^-10*diag(5)
      change<-ifelse(N_dose==FALSE,0,1)
      change1<-ifelse(log_scale==FALSE,0,1)

#Required arguments
	k1<-length(P)
    	p<-EDp 
	e1<-10^-7
	e2<-epsilon
	n<-1
	p1<-1
    	it<-1
	nit<-r
	gr<-grid
      lb<-LB
      LB<-ifelse(change==0,log(LB),log(-UB))
      UB<-ifelse(change==0,log(UB),log(-lb))

#Setting parameter values
      P[3]<-ifelse(change==0,log(P[3]),log(-P[3]))
      T<-P	

#Initial design points with weights	
	X<-c(LB,LB+(UB-LB)/4,LB+2*(UB-LB)/4,LB+3*(UB-LB)/4,UB)
	W<-rep(1/length(X),length(X)-1)

# Run V-algorithm to get initial design points
      
      M<-upinfor(W,T,X,k1)
      while(n<nit){
            x<-seq(LB,UB,gr)
            ds<-rep(0,length(x))
            if(exp(UB)<999) {
                 inv<-ginv(M)
            } else {
                 inv<-Inv(M,I5)
            }
            for (i in 1:length(x)) {
                 ds[i]<-DS1(T,x[i],inv,p,k1)
            }
            newX<-x[which.max(ds)]
            newds<-max(ds)
            an<-1/(n+1)
            newM<-(1-an)*M+an*f(T,newX,k1)%*%t(f(T,newX,k1))
            M<-newM
            X<-c(X,newX)
            n<-n+1
      }
      r<-length(X)
      X<-unique(X[(r-k1):r])
	X<-sort(X,decreasing=FALSE)

#Searching optimal design using the initial design selected
	cat(format("Computing the difference between the sensitivity function and the upper bound", width=80),"\n")
      while(p1>e2) {
            x<-seq(LB,UB,gr)
            ds<-rep(0,length(x))
            D<-S_weight(X,T,e1,c_weight,p,k1,UB,I5)
            X<-D[1,]
            k<-length(X)
            W<-D[2,1:k-1]
            M<-upinfor(W,T,X,k1)
            if(exp(UB)<999) {
                 inv<-ginv(M)
            } else {
                 inv<-Inv(M,I5)
            }
            for (i in 1:length(x)) {
                  ds[i]<-DS1(T,x[i],inv,p,k1)
            }
            newX<-x[which.max(ds)]
            newds<-max(ds)
            X<-c(X,newX)
            X<-unique(sort(X,decreasing=FALSE))
            newp<-abs(newds-1)
            if(abs(newp-p1)<.0000001) newp<-10^-20
            if(it>20) newp<-10^-20
            p1<-newp
            it<-it+1
		cat(p1,"\n")
      }

#Verification of the c-optimal design	
      X<-D[1,]
      W<-D[2,1:length(X)-1]
      x<-seq(LB,UB,gr)
      ds<-rep(0,length(x))
      M<-upinfor(W,T,X,k1)
      if(exp(UB)<999) {
            inv<-ginv(M)
      } else {
            inv<-Inv(M,I5)
      }

      for (i in 1:length(x)) 
            ds[i]<-DS1(T,x[i],inv,p,k1)

      if(change==0) {
            x<-exp(x)
            Dose<-round(exp(D[1,]),2)
      } else {
            x<--exp(x)
            Dose<--round(exp(D[1,]),2)
      }
	Weight<-round(D[2,],3)    
	D<-rbind(Dose,Weight)
      if(change1==1) {
	      plot(x,ds,log="x",cex=.3,ylab="Sensitivity function",xlab="dose")
      } else {
            plot(x,ds,cex=.3,ylab="Sensitivity function",xlab="dose")
      }
#Print optimal design rescaled on original dose level
      L<-list()
      L[[1]]<-D
	names(L)<-"c-optimal design"
      return(L)
}





