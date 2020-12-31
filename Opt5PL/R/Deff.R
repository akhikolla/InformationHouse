Deff<-function(weight,dose,model,P,LB,UB,grid=.01,N_dose=FALSE)
{
#checking negative dose levels
      change<-ifelse(N_dose==FALSE,0,1)

#requested arguments
	np<-model
	n<-1
      p<-1
	it<-1
	e1<-10^-7
	e2<-.001
	nit<-30
	gr<-grid
      I5<-10^-10*diag(5)
      I4<-10^-10*diag(4)
      I3<-10^-10*diag(3)
      lb<-LB
      LB<-ifelse(change==0,log(LB),log(-UB))
      UB<-ifelse(change==0,log(UB),log(-lb))

#Setting parameter values
      P[3]<-ifelse(change==0,log(P[3]),log(-P[3]))
      nP<-P
      T5<-nP
      T4<-c(nP[1:4],1)
      T3<-c(nP[1:3],0,1)
	T<-rbind(T3,T4,T5)

#Distinguishing for 3, 4, 5PL models
	if (np==3)
      	{q=c(1,0,0)} else if (np==4) {q=c(0,1,0)} else {q=c(0,0,1)}
	q1<-q[1]
	q2<-q[2]
	q3<-1-(q1+q2)

#Initial design points with weights
      X<-c(LB,LB+(UB-LB)/4,LB+2*(UB-LB)/4,LB+3*(UB-LB)/4,UB)
      wX<-length(X)
      W<-rep(1/wX,wX-1)

####################################################################################

# Run V-algorithm to get initial design points
      while(n<nit){
            x<-seq(LB,UB,gr)
            ds<-rep(0,length(x))
	      M3<-upinfor(W,T[1,],X,3)
            M4<-upinfor(W,T[2,],X,4)
            M5<-upinfor(W,T[3,],X,5)
            inv3<-Inv(M3,I3)
            inv4<-Inv(M4,I4)
            inv5<-Inv(M5,I5)
 	      for (i in 1:length(x)) {
                 ds[i]<-q1*smallds1(T[1,],x[i],inv3,3)+q2*smallds1(T[2,],x[i],inv4,4)+q3*smallds1(T[3,],x[i],inv5,5)
            }
            newX<-x[which.max(ds)]
            newds<-max(ds)
            p<-abs(newds-1)
            X<-c(X,newX)
		W<-c(W,1-sum(W))
		newW<-(1-1/(n+1))*W
		W<-newW
            n<-n+1
      }
	X<-sort(unique(X[(length(X)-wX):length(X)]),decreasing=F)

#Searching optimal design using the initial design selected
      while(p>e2) {
            x<-seq(LB,UB,gr)
            ds<-rep(0,length(x))
            D<-S_weight(X,T,e1,D_weight,c(q1,q2,q3))
            X<-D[1,]
            W<-D[2,1:(length(X)-1)]
	      M3<-upinfor(W,T[1,],X,3)
            M4<-upinfor(W,T[2,],X,4)
            M5<-upinfor(W,T[3,],X,5)
            inv3<-Inv(M3,I3)
            inv4<-Inv(M4,I4)
            inv5<-Inv(M5,I5)
 	      for (i in 1:length(x))
                ds[i]<-q1*smallds1(T[1,],x[i],inv3,3)+q2*smallds1(T[2,],x[i],inv4,4)+q3*smallds1(T[3,],x[i],inv5,5)
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
	Dpoints<-D[1,]
	Dlength<-length(Dpoints)
	Dweight<-D[2,1:(Dlength-1)]
	if (change==0) {
            Dose<-round(exp(D[1,]),2)
      }
      else {
            Dose<--round(exp(D[1,]),2)
      }
	Weight<-round(D[2,],3)
#D-optimal design
	Dopt<-rbind(Dose,Weight)

#Design given
      for (i in 1:length(dose)){
          if (dose[i]==0) dose[i]<-.0001
      }
      if (change==1) dose<--dose
      dose<-log(dose)

      eff<-(det(upinfor(weight,T[(np-2),],dose,np))/det(upinfor(Dweight,T[(np-2),],Dpoints,np)))^(1/np)
	size<-100*(1/eff-1)
	eff<-round(eff,2)
	size<-round(size,2)
	cat(format("D-optimal design", width=50),"\n")
      print(Dopt)
      cat(format("D-efficiency", width=50),"\n")
      print(eff)
	cat(format("%More Samples Needed", width=50),"\n")
	print(size)
}

