################################
#### HARMONIC FIT LIBRARY ######
################################

# CONTAINS
  # get.period.stats() - determine stats for each period (wet & dry or partioned wet or dry)
  # fit.harmonic.opts()- determine paramters of harmonic using procedure in Wilks, Statistical Methods in Atmospheric Sciences
  # harfuncEst() - function used in fit.harmonic() to calculate amp & phase angle pars
  # harmonicFunc() - applies harmonic pars to produce the resultant curve
  # calc.phase.ang -phase angle par using procedure in Wilks, Statistical Methods in Atmospheric Sciences
    # get.par.sets() - calculates the vectors of m & sd pars
    # get.residuals.wd.conditional() - (for temp or other pars)

#--------------------------------------------------------------------------------------
#for variables conditional on wet/dry status
#determine stats for each period (wet & dry or partioned wet or dry)
get.period.stats<-function(data=NULL,            #ts vector for data
                           nperiod=NULL,         #no. periods
                           i.hh=NULL,            #harmonic indices
                           sep.wd=TRUE,          #separate wet and dry?
                           omit.drypars=FALSE,   #choice to not calculate pars on dry days
                           rain=NULL,            #vector of rain
                           threshold=0           #set wet-dry threshold to be 0
                           ){
  out=list()
  
  if(sep.wd==TRUE){             #separate calculations for wet and dry days
    out$wet$m=rep(NA,nperiod);out$wet$sd=rep(NA,nperiod);out$wet$cv=rep(NA,nperiod)
    if(omit.drypars==FALSE){out$dry$m=rep(NA,nperiod);out$dry$sd=rep(NA,nperiod);out$dry$cv=rep(NA,nperiod)}
    for(h in 1:nperiod){
     w.ind=i.hh[[h]][which(rain[i.hh[[h]]]>threshold)]
     #calculate mean, sd and cv for wet series
     out$wet$m[h]=mean(data[w.ind],na.rm=TRUE)
     out$wet$sd[h]=sd(data[w.ind],na.rm=TRUE)
     #out$wet$cv[h]=out$wet$sd[h]/out$wet$m[h]
     
     if(omit.drypars==FALSE){
       #calculate mean, sd and cv for dry series
       d.ind=i.hh[[h]][which(rain[i.hh[[h]]]<=threshold)]
       out$dry$m[h]=mean(data[d.ind],na.rm=TRUE)
       out$dry$sd[h]=sd(data[d.ind],na.rm=TRUE)
       #out$dry$cv[h]=out$dry$sd[h]/out$dry$m[h]
     }
       
    }
  }else{                        #calculate for all days together
    out$all$m=rep(NA,nperiod);out$all$sd=rep(NA,nperiod);out$all$cv=rep(NA,nperiod)
    for(h in 1:nperiod){
      out$all$m[h]=mean(data[i.hh[[h]]],na.rm=TRUE)
      out$all$sd[h]=sd(data[i.hh[[h]]],na.rm=TRUE)
      #out$all$cv[h]=out$all$sd[h]/out$all$m[h]
    }
  }
  return(out)  
}
#out=get.period.stats(data=obs,nperiod=26,i.hh=i.hh,sep.wd=TRUE,omit.drypars=TRUE,rain=obs)

#TAKEN FROM WGEN MANUAL/statistical methods fro atmospheric sciences, D.S. Wilks 2005 pp374-382
harfuncEst<-function(f=cos,            #function (cos or sine)     
                  nperiod=NULL,        #no. of periods
                  period=NULL,         #current periods
                  dat=NULL,            #value for the periods
                  est.mean=NULL,       #estimated mean across all periods
                  k=1                  #harmonic number
                  ){
  
  x1=(dat-est.mean)*f(2*pi*period*k/nperiod)
  return(x1)
  
}

#TAKEN FROM WGEN MANUAL/statistical methods fro atmospheric sciences, D.S. Wilks 2005 pp374-382
fit.harmonic.opts<-function(nperiod=NULL,    # no. of periods 
                            v.stat=NULL,     # vector of stat for all periods, length=nperiod
                            k=1              # no. of harmonics to apply
){
  
  #Determine means
  m_bar=mean(v.stat,na.rm=TRUE)
  
  #Make space to store amplitude and phase angle
  c=rep(NA,k)
  ang=rep(NA,k)
  
  #LOOP OVER DIFFERENT HARMONICS 1 to k
  for(i in 1:k){
    #Calc estimators using cos & sin for first harmonic
    m_a=sum(harfuncEst(f=cos,nperiod=nperiod,period=seq(1,nperiod),dat=v.stat,est.mean=m_bar,k=i)) 
    m_b=sum(harfuncEst(f=sin,nperiod=nperiod,period=seq(1,nperiod),dat=v.stat,est.mean=m_bar,k=i))
    aa=m_a*2.0/nperiod
    bb=m_b*2.0/nperiod
    
    #get phase angle & amplitue for first harmonic
    c[i]=sqrt(aa^2.0+bb^2.0) 
    ang[i]=calc.phase.ang(aa=aa,bb=bb)
    rm(m_a,m_b,aa,bb)
  }
  
  out=list(mean=m_bar,amp=c,phase.ang=ang)
  return(out)
}
#phase angle calculator
#statistical methods for atmospheric sciences, D.S. Wilks 2005 pp374-382 (note p379)
calc.phase.ang<-function(aa,     
                         bb
){
  if(aa>0){
    ang=atan(bb/aa)
    if(ang<0){ang=ang+pi}    #make positive if less than zero by adding pi
  }
  if(aa<0){
    angA=atan(bb/aa)+pi; angB=atan(bb/aa)-pi
    tmp=c(angA,angB)
    ind=which((tmp>0) & (tmp<(2*pi))) #choose the one between 0 & 2pi
    ang=tmp[ind]
  }
  if(aa==0){ang=pi/2.0}
  return(ang)
}

#Calulates points on harmonic curve - Allows higher harmonics (up to 3)
harmonicFunc<- function(x,           # series of times
                        mean,        # mean
                        amp,         # vector of amplitudes
                        phase.ang,   # vector of phase angles
                        k,            # no. of harmonics (max 3)
                        nperiod=26   #No. of periods 
                        ) {
  # First line of switch modified by SW 31/12/2019 for numerical speedup. To be extended to additional harmonics, and also better handling of leap years.
  switch(k,
         y <- rep(mean+amp[1]*cos(2*pi*x[1:nperiod]/nperiod-phase.ang[1]), length.out=length(x)),
         y <- mean+amp[1]*cos(2*pi*x/nperiod-phase.ang[1])+amp[2]*cos(2*pi*x/nperiod-phase.ang[2]),
         y <- mean+amp[1]*cos(2*pi*x/nperiod-phase.ang[1])+amp[2]*cos(2*pi*x/nperiod-phase.ang[2])+amp[3]*cos(2*pi*x/nperiod-phase.ang[3]),
         -999           # default
         )
  return(y)
}



#EXAMPLE 
# pp=out$wet$cv
# plot(x=seq(1,nperiod),y=pp,type="p",col="blue",pch=1)
# kk=1
# test=fit.harmonic.opts(nperiod=26, v.stat=pp, k=kk)
# test.har=harmonicFunc(x=seq(1,nperiod),mean=test$mean,amp=test$amp,phase.ang=test$amp,k=kk)
# lines(x=seq(1,nperiod),y=test.har,col="red")
# kk=2
# test=fit.harmonic.opts(nperiod=26, v.stat=pp, k=kk)
# test.har=harmonicFunc(x=seq(1,nperiod),mean=test$mean,amp=test$amp,phase.ang=test$amp,k=kk)
# lines(x=seq(1,nperiod),y=test.har,col="purple")
# kk=3
# test=fit.harmonic.opts(nperiod=26, v.stat=pp, k=kk)
# test.har=harmonicFunc(x=seq(1,nperiod),mean=test$mean,amp=test$amp,phase.ang=test$amp,k=kk)
# lines(x=seq(1,nperiod),y=test.har,col="green")



#calculation of residual elements
 #for temp - wet and dry days separated
  #equation (4) Richardson 1981
#THIS RELATES TO TEMP & OTHER SERIES
get.par.sets<-function(nperiod=NULL,        #No. of periods in harmonic
                       Hpar_m=NULL,         #harmonic pars for means (mean, amp, phase.angle)
                       Hpar_sd=NULL,        #harmonic pars for means (mean, amp, phase.angle)
                       k=NULL               # no. of harmonics (max 3)
                       ){
  if(k>3){
    logfile("ERROR: No. of cycles > 3. Invalid selection.")
    stop("No. of chosen cycles too large")
  }
  
  #get period means and sd
  mPar=harmonicFunc(x=seq(1,nperiod),mean=Hpar_m$mean,amp=Hpar_m$amp,phase.ang=Hpar_m$phase.ang,k)
  sdPar=harmonicFunc(x=seq(1,nperiod),mean=Hpar_sd$mean,amp=Hpar_sd$amp,phase.ang=Hpar_sd$phase.ang,k)
  
  out=list(mPar=mPar,sdPar=sdPar)
  return(out)
}


#FINISH FOR TEMPERATURE & OTHER SERIES
# get.residuals.wd.conditional<-function(data=NULL,            #ts vector for data
#                                        nperiod=NULL,                        #no. periods
#                                        i.hh=NULL                           #harmonic indices
# 
# 
# ){
#   
#   
#   est.stats=get.period.stats(data=data,nperiod=nperiod,i.hh=i.hh,rain=NULL,sep.wd=TRUE,omit.drypars=FALSE,threshold=0)
# 
#   fit.harmonic.opts()  #do twice for each wet and dry
#   
#   get.pars()     #do twice   - store these logically for later use (save time)
# #need good naming convention for these
#   
#   res=rep(NA,length(data))
#   #do as vector operations
#   #if wet day  res=(x - xbar_w[period])/sd_w[period]
#   #if dry day  res=(x - xbar_d[period])/sd_d[period]
#   
#   for(h in 1:nperiod){
#     w.ind=i.hh[[h]][which(rain[i.hh[[h]]]>threshold)]   #this has be calculated before hmmm
#     d.ind=i.hh[[h]][which(rain[i.hh[[h]]]<=threshold)]
#     res[w.ind]=(data[w.ind]-)/
#     
# 
#   }
# 
# }
#correlations between residual elements (how many climate variables?)
#this will change model choice

#-----------------------------------------------------------------

