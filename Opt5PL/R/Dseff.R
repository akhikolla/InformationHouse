Dseff<-function(weight,dose,P,LB,UB,r=30,epsilon=.001,grid=.01, N_dose=FALSE)
{
	change<-ifelse(N_dose==FALSE,0,1)

	#Required arguments
      n<-1
      p<-1
      it<-1
      e2<-epsilon
	nit<-r
	gr<-.01
      e1<-10^-7
	I5<-10^-10*diag(5)
	I4<-10^-10*diag(4)
      lb<-LB
      LB<-ifelse(change==0,log(LB),log(-UB))
      UB<-ifelse(change==0,log(UB),log(-lb))

#Setting parameter values
      P[3]<-ifelse(change==0,log(P[3]),log(-P[3]))
      T<-P

#Initial design points with weights
	X<-c(LB,LB+(UB-LB)/4,LB+2*(UB-LB)/4,LB+3*(UB-LB)/4,UB)
	wX<-length(X)
	W<-rep(1/wX,wX-1)

####Required functions to run the algorithm#######################

# Run V-algorithm to get initial design points
      M5<-upinfor(W,T,X,5)
	M4<-upinfor(W,T,X,4)

      while(n<nit){
            x<-seq(LB,UB,gr)
            ds<-rep(0,length(x))
            inv1<-Inv(M4,I4)
            inv<-Inv(M5,I5)
            for (i in 1:length(x))
                 ds[i]<-ds11(T,x[i],inv,inv1,5)
            newX<-x[which.max(ds)]
            newds<-max(ds)
            an<-1/(n+1)
            p<-abs(newds-1)
            newM5<-(1-an)*M5+an*f(T,newX,5)%*%t(f(T,newX,5))
		newM4<-(1-an)*M4+an*f(T,newX,4)%*%t(f(T,newX,4))
            M5<-newM5
		M4<-newM4
            X<-c(X,newX)
            n<-n+1
      }
      r<-length(X)
	X<-sort(unique(X[(r-wX):r]),decreasing=F)

#Searching optimal design using the initial design selected
      while(p>e2) {
            x<-seq(LB,UB,gr)
            ds<-rep(0,length(x))
            D<-S_weight(X,T,e1,DD_weight,I4,I5,5)
            X<-D[1,]
            W<-D[2,1:length(X)-1]
            M4<-upinfor(W,T,X,4)
            M5<-upinfor(W,T,X,5)
            inv1<-Inv(M4,I4)
            inv<-Inv(M5,I5)
            for (i in 1:length(x))
                  ds[i]<-ds11(T,x[i],inv,inv1,5)
            newX<-x[which.max(ds)]
            newds<-max(ds)
            X<-c(X,newX)
            X<-unique(sort(X,decreasing=F))
            newp<-abs(newds-1)
            if(abs(newp-p)<.0000001)   newp<-10^-20
            if(it>20)   newp<-10^-20
            p<-newp
            it<-it+1
      }
	Dspoints<-D[1,]
	Dsweight<-D[2,1:(length(Dspoints)-1)]
	if (change==0) Dose<-round(exp(D[1,]),2) else Dose<--round(exp(D[1,]),2)
	Weight<-round(D[2,],3)

#Ds-optimal design
	Dsopt<-rbind(Dose,Weight)

#Design given
      k<-length(dose)
      for (i in 1:k)
            if(dose[i]==0) dose[i]<-.0001
      if(change==1) dose<--dose
      dose<-log(dose)

	eff<-(det(upinfor(weight,T,dose,5))/det(upinfor(weight,T,dose,4)))/(det(upinfor(Dsweight,T,Dspoints,5))/det(upinfor(Dsweight,T,Dspoints,4)))
	size<-100*(1/eff-1)
	eff<-round(eff,2)
	size<-round(size,2)
	cat(format("Ds-optimal design", width=50),"\n")
      print(Dsopt)
      cat(format("Ds-efficiency", width=50),"\n")
      print(eff)
	cat(format("%More Samples Needed", width=50),"\n")
	print(size)
}



