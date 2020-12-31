weightedbinCItest <- function(x,y,S,suffStat) {
  ## p-value:
  CIbincoretest(x, y, S, suffStat)
}

weightedcatCItest <- function(x,y,S,suffStat) {
  ## p-value:
  CIcatcoretest(x, y, S, suffStat)
}

#user defined CI test
usrCItest <- function(x,y,S,suffStat) {
  ## p-value:
  #insert test here (I have inserted categorical data CI test)
  CIcatcoretest(x, y, S, suffStat)
}

#user defined core CI test
CIusrcoretest<-function(j,k,parentnodes,suffStat) {
  #insert core function here (I have inserted categorical data CI test)
  lp <- length(parentnodes) # number of parents
  
  Cj <- suffStat$Cvec[j] # number of levels of j
  Ck <- suffStat$Cvec[k] # number of levels of k
  
  switch(as.character(lp),
         "0"={# no parents
           Cp <- 1 # effectively 1 parent level
           summys <- rep(0, nrow(suffStat$data))
         },
         "1"={# one parent
           Cp <- suffStat$Cvec[parentnodes] # number of parent levels
           summys <- suffStat$data[,parentnodes]
         },     
         { # more parents
           Cp <- prod(suffStat$Cvec[parentnodes])
           # use mixed radix mapping to unique parent states
           summys<-colSums(cumprod(c(1,suffStat$Cvec[parentnodes[-lp]]))*t(suffStat$data[,parentnodes]))
         })
  
  datasummy <- colSums(c(1,Cj)*t(suffStat$data[,c(j,k)]))
  
  # tabulate the observed counts
  
  if(!is.null(suffStat$weightvector)){
    Ns <- collectCcatwt(summys, datasummy, suffStat$weightvector, Cp, Cj*Ck)
  } else{
    Ns <- collectCcat(summys, datasummy, Cp, Cj*Ck)
  }
  
  ### We could run over the data again for the marginals  
  # margj <- collectCcat(summys, suffStat$data[,j], Cp, Cj)
  # margk <- collectCcat(summys, suffStat$data[,k], Cp, Ck) 
  
  ### instead we combine the sums already  
  
  margj <- Ns[,1:Cj] + Ns[,Cj+1:Cj]
  if(Ck>2){
    for(ii in 3:Ck-1){
      margj <- margj + Ns[,ii*Cj+1:Cj]
    }
  }
  
  margk <- Ns[,Cj*(1:Ck-1)+1] + Ns[,Cj*(1:Ck-1)+2]
  if(Cj>2){
    for(ii in 3:Cj-1){
      margk <- margk + Ns[,Cj*(1:Ck-1)+ii+1]
    }
  }
  
  if(lp==0){# then we have vectors!
    Es <- c(margj*margk[1],margj*margk[2])
    if(Ck>2){
      for(ii in 3:Ck){
        Es <- c(Es,margj*margk[ii])
      }
    }
  } else {
    Es <- cbind(margj*margk[,1],margj*margk[,2])
    if(Ck>2){
      for(ii in 3:Ck){
        Es <- cbind(Es,margj*margk[,ii])
      }
    }
  }
  Es <- Es/rowSums(Ns) # normalise
  
  G2 <- 2*Ns*log(Ns/Es)
  G2[is.nan(G2)] <- 0
  
  Gsquared <- sum(G2)
  df <- Cp*(Cj-1)*(Ck-1)
  
  pvally<-pchisq(Gsquared, df, lower.tail = FALSE)
  
  return(pvally)
}

# this version uses the C code and runs through the data 4 times
CIbincoretest<-function(j,k,parentnodes,suffStat){
  
  lp<-length(parentnodes) # number of parents
  noparams<-2^lp # number of binary states of the parents
  
  switch(as.character(lp),
         "0"={# no parents # note the weighting is already in the first term
           P11<-sum(suffStat$d1[,j]*suffStat$data[,k],na.rm=TRUE)
           P10<-sum(suffStat$d1[,j]*(1-suffStat$data[,k]),na.rm=TRUE)
           P01<-sum(suffStat$d0[,j]*suffStat$data[,k],na.rm=TRUE)
           P00<-sum(suffStat$d0[,j]*(1-suffStat$data[,k]),na.rm=TRUE)
           PT<-P11+P10+P01+P00
           N1<-P11+P10
           N0<-P01+P00
           NT<-N0+N1
           M1<-P11+P01
           M0<-P10+P00
           MT<-M0+M1
           
           # calculate the statistic in the tedious way    
           
           if(P11>0){
             part1<-P11*(log(P11)+log(PT)-log(N1)-log(M1))
           } else{
             part1<-0
           }
           if(P10>0){
             part2<-P10*(log(P10)+log(PT)-log(N1)-log(M0))
           } else{
             part2<-0
           }
           if(P01>0){
             part3<-P01*(log(P01)+log(PT)-log(N0)-log(M1))
           } else{
             part3<-0
           }
           if(P00>0){
             part4<-P00*(log(P00)+log(PT)-log(N0)-log(M0))
           } else{
             part4<-0
           }
           
           Gsquared<-2*(part1+part2+part3+part4)
           
         },
         "1"={ # 1 parent
           Gsquared<-0  
           summys<-suffStat$data[,parentnodes]
           
           for(i in 1:noparams-1){
             totest<-which(summys==i) # note the weighting is already in the first term
             P11<-sum(suffStat$d1[totest,j]*suffStat$data[totest,k],na.rm=TRUE)
             P10<-sum(suffStat$d1[totest,j]*(1-suffStat$data[totest,k]),na.rm=TRUE)
             P01<-sum(suffStat$d0[totest,j]*suffStat$data[totest,k],na.rm=TRUE)
             P00<-sum(suffStat$d0[totest,j]*(1-suffStat$data[totest,k]),na.rm=TRUE)
             PT<-P11+P10+P01+P00
             N1<-P11+P10
             N0<-P01+P00
             NT<-N0+N1
             M1<-P11+P01
             M0<-P10+P00
             MT<-M0+M1
             
             # calculate the statistic in the tedious way    
             
             if(P11>0){
               part1<-P11*(log(P11)+log(PT)-log(N1)-log(M1))
             } else{
               part1<-0
             }
             if(P10>0){
               part2<-P10*(log(P10)+log(PT)-log(N1)-log(M0))
             } else{
               part2<-0
             }
             if(P01>0){
               part3<-P01*(log(P01)+log(PT)-log(N0)-log(M1))
             } else{
               part3<-0
             }
             if(P00>0){
               part4<-P00*(log(P00)+log(PT)-log(N0)-log(M0))
             } else{
               part4<-0
             }
             
             Gsquared<-Gsquared+2*(part1+part2+part3+part4)
           }
           
         },
         { # more parents
           
           summys<-colSums(2^(c(0:(lp-1)))*t(suffStat$data[,parentnodes]))
           
           P11vec<-suffStat$d1[,j]*suffStat$data[,k] # these include the weighting in the first term
           P01vec<-suffStat$d0[,j]*suffStat$data[,k]
           
           tokeep<-which(!is.na(summys+P11vec)) # remove NAs either in the parents or the children
           if(length(tokeep)<length(summys)){
             P11s<-collectC(summys[tokeep],P11vec[tokeep],noparams)
             P10s<-collectC(summys[tokeep],suffStat$d1[tokeep,j],noparams)-P11s
             P01s<-collectC(summys[tokeep],P01vec[tokeep],noparams)
             P00s<-collectC(summys[tokeep],suffStat$d0[tokeep,j],noparams)-P01s
           } else {
             P11s<-collectC(summys,P11vec,noparams)
             P10s<-collectC(summys,suffStat$d1[,j],noparams)-P11s
             P01s<-collectC(summys,P01vec,noparams)
             P00s<-collectC(summys,suffStat$d0[,j],noparams)-P01s
           }
           
           PTs<-P11s+P10s+P01s+P00s # collect all marginal counts
           N1s<-P11s+P10s
           N0s<-P01s+P00s
           NTs<-N0s+N1s
           M1s<-P11s+P01s
           M0s<-P10s+P00s
           MTs<-M0s+M1s
           
           # calculate the statistic in the tedious way    
           
           part1s<-P11s*(log(P11s)+log(PTs)-log(N1s)-log(M1s))
           part1s[which(P11s==0)]<-0 # put in limit by hand
           part1<-sum(part1s)
           
           part2s<-P10s*(log(P10s)+log(PTs)-log(N1s)-log(M0s))
           part2s[which(P10s==0)]<-0
           part2<-sum(part2s)
           
           part3s<-P01s*(log(P01s)+log(PTs)-log(N0s)-log(M1s))
           part3s[which(P01s==0)]<-0
           part3<-sum(part3s)
           
           part4s<-P00s*(log(P00s)+log(PTs)-log(N0s)-log(M0s))
           part4s[which(P00s==0)]<-0
           part4<-sum(part4s)
           
           Gsquared<-2*(part1+part2+part3+part4)
           
         })
  
  pvally<-pchisq(Gsquared, noparams, lower.tail = FALSE)
  
  return(pvally)
}



# this version uses the C code and runs through the data once

CIcatcoretest<-function(j,k,parentnodes,suffStat){
  
  lp <- length(parentnodes) # number of parents
  
  Cj <- suffStat$Cvec[j] # number of levels of j
  Ck <- suffStat$Cvec[k] # number of levels of k
  
  switch(as.character(lp),
         "0"={# no parents
           Cp <- 1 # effectively 1 parent level
           summys <- rep(0, nrow(suffStat$data))
         },
         "1"={# one parent
           Cp <- suffStat$Cvec[parentnodes] # number of parent levels
           summys <- suffStat$data[,parentnodes]
         },     
         { # more parents
           Cp <- prod(suffStat$Cvec[parentnodes])
           # use mixed radix mapping to unique parent states
           summys<-colSums(cumprod(c(1,suffStat$Cvec[parentnodes[-lp]]))*t(suffStat$data[,parentnodes]))
         })
  
  datasummy <- colSums(c(1,Cj)*t(suffStat$data[,c(j,k)]))
  
  # tabulate the observed counts
  
  if(!is.null(suffStat$weightvector)){
    Ns <- collectCcatwt(summys, datasummy, suffStat$weightvector, Cp, Cj*Ck)
  } else{
    Ns <- collectCcat(summys, datasummy, Cp, Cj*Ck)
  }
  
  ### We could run over the data again for the marginals  
  # margj <- collectCcat(summys, suffStat$data[,j], Cp, Cj)
  # margk <- collectCcat(summys, suffStat$data[,k], Cp, Ck) 
  
  ### instead we combine the sums already  
  
  margj <- Ns[,1:Cj] + Ns[,Cj+1:Cj]
  if(Ck>2){
    for(ii in 3:Ck-1){
      margj <- margj + Ns[,ii*Cj+1:Cj]
    }
  }
  
  margk <- Ns[,Cj*(1:Ck-1)+1] + Ns[,Cj*(1:Ck-1)+2]
  if(Cj>2){
    for(ii in 3:Cj-1){
      margk <- margk + Ns[,Cj*(1:Ck-1)+ii+1]
    }
  }
  
  if(lp==0){# then we have vectors!
    Es <- c(margj*margk[1],margj*margk[2])
    if(Ck>2){
      for(ii in 3:Ck){
        Es <- c(Es,margj*margk[ii])
      }
    }
  } else {
    Es <- cbind(margj*margk[,1],margj*margk[,2])
    if(Ck>2){
      for(ii in 3:Ck){
        Es <- cbind(Es,margj*margk[,ii])
      }
    }
  }
  Es <- Es/rowSums(Ns) # normalise
  
  G2 <- 2*Ns*log(Ns/Es)
  G2[is.nan(G2)] <- 0
  
  Gsquared <- sum(G2)
  df <- Cp*(Cj-1)*(Ck-1)
  
  pvally<-pchisq(Gsquared, df, lower.tail = FALSE)
  
  return(pvally)
}




