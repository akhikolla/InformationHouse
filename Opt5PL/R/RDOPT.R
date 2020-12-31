RDOPT<-function(LB,UB,P3=0,P4=0,P5,q,grid=.01,r=30,epsilon=.001,N_dose=FALSE,log_scale=TRUE)
{
      change<-ifelse(N_dose==FALSE,0,1)
      change1<-ifelse(log_scale==FALSE,0,1)

#Required arguments
      n<-1
      p<-1
	it<-1
	e1<-10^-7
      I5<-10^-10*diag(5)
      I4<-10^-10*diag(4)
      I3<-10^-10*diag(3)
      lb<-LB
      LB<-ifelse(change==0,log(LB),log(-UB))
      UB<-ifelse(change==0,log(UB),log(-lb))

#Setting parameter values

      if(length(P3)==1) P3<-c(P5[1:3],0,1)
      if(length(P4)==1) P4<-c(P5[1:4],1)
      if(change==0) {
          P3<-c(P3[1:2],log(P3[3]),0,1)
          P4<-c(P4[1:2],log(P4[3]),P4[4],1)
          P5<-c(P5[1:2],log(P5[3]),P5[4:5])
      } else {
          P3<-c(P3[1:2],log(-P3[3]),0,1)
          P4<-c(P4[1:2],log(-P4[3]),P4[4],1)
          P5<-c(P5[1:2],log(-P5[3]),P5[4:5])
      }

      T5<-P5[1:5]
      T4<-c(P4[1:4],1)
      T3<-c(P3[1:3],0,1)
      T<-rbind(T3,T4,T5)

#Setting the weights for 3, 4, 5PL models
      q1<-q[1]
      q2<-q[2]
      q3<-1-(q1+q2)

#Initial design points with weights
      X<-c(LB,LB+(UB-LB)/4,LB+2*(UB-LB)/4,LB+3*(UB-LB)/4,UB)
      wX<-length(X)
      W<-rep(1/wX,wX-1)

####################################################################################

# Run V-algorithm to get initial design points
      while(n<r){
            x<-seq(LB,UB,grid)
            ds<-rep(0,length(x))
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
            p<-abs(newds-1)
            X<-c(X,newX)
		W<-c(W,1-sum(W))
		newW<-(1-1/(n+1))*W
		W<-newW
            n<-n+1
      }
	X<-sort(unique(X[(length(X)-wX):length(X)]),decreasing=FALSE)

#Searching optimal design using the initial design selected
      cat(format("Computing the difference between the sensitivity function and the upper bound", width=80),"\n")
      while(p>epsilon) {
            x<-seq(LB,UB,grid)
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
            X<-unique(sort(X,decreasing=FALSE))
            newp<-abs(newds-1)
            if(abs(newp-p)<.0000001) newp<-10^-20
            if(it>20) newp<-10^-20
            p<-newp
            it<-it+1
            cat(p,"\n")
      }

#Verification of the optimal design
	X<-D[1,]
      W<-D[2,1:(length(X)-1)]
      x<-seq(LB,UB,.001)
      ds<-rep(0,length(x))
	maxds<-rep(0,length(X))
	M3<-upinfor(W,T[1,],X,3)
      M4<-upinfor(W,T[2,],X,4)
      M5<-upinfor(W,T[3,],X,5)
      inv3<-Inv(M3,I3)
      inv4<-Inv(M4,I4)
      inv5<-Inv(M5,I5)
 	for (i in 1:length(x))
            ds[i]<-q1*smallds1(T[1,],x[i],inv3,3)+q2*smallds1(T[2,],x[i],inv4,4)+q3*smallds1(T[3,],x[i],inv5,5)
	for (i in 1:length(X))
            maxds[i]<-q1*smallds1(T[1,],X[i],inv3,3)+q2*smallds1(T[2,],X[i],inv4,4)+q3*smallds1(T[3,],X[i],inv5,5)
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
      names(L)<-"D-optimal design"
      return(L)
}


