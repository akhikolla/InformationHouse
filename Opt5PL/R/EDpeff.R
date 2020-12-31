EDpeff<-function(weight,dose,P,EDp,LB,UB,r=30,epsilon=.001,grid=.01, N_dose=FALSE)
{
#checking negative dose levels
      change<-ifelse(N_dose==F,0,1)

#Required arguments
	k1<-length(P)
    	#EDp(p percentage)
    	p<-EDp
	e1<-10^-7
	e2<-epsilon
	n<-1
	p1<-1
    	it<-1
	nit<-r
	gr<-grid
      I5<-10^-10*diag(5)
      lb<-LB
      LB<-ifelse(change==0,log(LB),log(-UB))
      UB<-ifelse(change==0,log(UB),log(-lb))

#Setting parameter values
      P[3]<-ifelse(change==0,log(P[3]),log(-P[3]))
      T<-P

#Initial design points with weights
	X<-c(LB,LB+(UB-LB)/4,LB+2*(UB-LB)/4,LB+3*(UB-LB)/4,UB)
	W<-rep(1/length(X),length(X)-1)

####Required functions to run the algorithm#######################

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
	X<-sort(X,decreasing=F)

#Searching optimal design using the initial design selected
	while(p1>e2) {
            x<-seq(LB,UB,gr)
            ds<-rep(0,length(x))
            D<-S_weight(X,T,e1,c_weight,p,k1,UB,I5)
            X<-D[1,]
            k<-length(X)
            W<-D[2,1:k-1]
            if(exp(UB)<999) {
                  inv<-ginv(upinfor(W,T,X,k1))
            } else {
                  inv<-Inv(upinfor(W,T,X,k1),I5)
            }
            for (i in 1:length(x)) {
                  ds[i]<-DS1(T,x[i],inv,p,k1)
            }
            newX<-x[which.max(ds)]
            newds<-max(ds)
            X<-c(X,newX)
            X<-unique(sort(X,decreasing=F))
            newp<-abs(newds-1)
            if(abs(newp-p1)<.0000001) newp<-10^-20
            if(it>20) newp<-10^-20
            p1<-newp
            it<-it+1
      }
	Cpoints<-D[1,]
	#Clength=length(Cpoints)
	Cweight<-D[2,1:(length(Cpoints)-1)]
	if (change==0) {
            Dose<-round(exp(D[1,]),2)
      } else {
		Dose<--round(exp(D[1,]),2)
      }
	Weight<-round(D[2,],3)

#C-optimal design
	Copt<-rbind(Dose,Weight)

#Design given
      k<-length(dose)
      for (i in 1:k)
            if(dose[i]==0) dose[i]<-.0001
      if(change==1) dose<--dose
      dose<-log(dose)
      if(exp(UB)<999) {
            inv<-ginv(upinfor(weight,T,dose,k1))
	      inv_min<-ginv(upinfor(Cweight,T,Cpoints,k1))
      } else {
            inv<-Inv(upinfor(weight,T,dose,k1),I5)
            inv_min<-Inv(upinfor(Cweight,T,Cpoints,k1),I5)
      }
	eff<-(t(g(T,p))%*%inv_min%*%g(T,p))/(t(g(T,p))%*%inv%*%g(T,p))
	size<-100*(1/eff-1)
	eff<-round(eff,2)
	size<-round(size,2)
	cat(format("c-optimal design", width=50),"\n")
      print(Copt)
      cat(format("c-efficiency", width=50),"\n")
      print(eff)
	cat(format("%More Samples Needed", width=50),"\n")
	print(size)
}

