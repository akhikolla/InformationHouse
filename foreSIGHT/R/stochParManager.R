 ################################
#### STOCHASTIC PAR MANAGER ####
################################

#CONTAINS
  #get.multi.model.info() - get modle info for multiple modelTags
  #par.manager() - based on model tag converts input vector of pars into pdd, pwd, alpha & beta vectors of nperiod length
  #init.calib() - initial calibration function. Uses obs data to fit model to baseline climate.
    #uses model tag to determine what should be fitted (e.g. how many periods, using harmonic, etc.)
    #gammaParsMLE2() - args(dat, wetThresh)
    #pdd.pwd.estimator() -args (dat,ind,threshold)
  #simHarTS.parmanager() - based on model tag converts input vector of pars into cor0, cor1, Hmean, Hsd vectors of period length
  #whichPars
  #update.simPriority() -update modelTag via simPriority
#-------------------------------------------------------------------------------------------------------------
#Get info for multiple models
 get.multi.model.info<-function(modelTag=NULL){
   nMod=length(modelTag)   
   if(nMod==1){
     modelInfo=list()
     modelInfo[[modelTag[1]]]=get.model.info(modelTag[1])                     #even if 1 model still stored in list format
   }else{
     modelInfo=sapply(X = modelTag,FUN=get.model.info,USE.NAMES=TRUE)
   }
   return(modelInfo)
 }
 

#UPDATE MODEL INFO IF FIXED PARAMETERS
update.model.info<-function(modelTag=NULL, modelInfo=NULL,fixedPars=NULL,minUserBound=NULL,maxUserBound=NULL,file=NULL){
  #check if model info correct
  if((modelTag =="P-har12-wgen-FS")){  #only for FS currently
    #Check if the correct number of fixed parameters is supplied
    if(length(fixedPars)!=4){
      logfile("Error: Incorrect number of fixed parameters values provided in fixedPar (4 needed)",file)
      logfile("Program terminated",file)
      stop("Incorrect number of fixed parameters values provided in fixedPar (4 needed)")
    }  
    
    temp_modelInfo=get.model.info(modelTag="P-har12-wgen")
    modelInfo$parNam=temp_modelInfo$parNam
    modelInfo$npars=temp_modelInfo$npars
    modelInfo$fixedPars=NA   #update to no fixed pars
    #Insert fixed pars into vector quick way (only 1 option uses this so far)
    make_minBound=c(modelInfo$minBound[1:2],fixedPars[1],
                    modelInfo$minBound[3:4],fixedPars[2],
                    modelInfo$minBound[5:6],fixedPars[3],
                    modelInfo$minBound[7:8],fixedPars[4])

    make_maxBound=c(modelInfo$maxBound[1:2],fixedPars[1],
                    modelInfo$maxBound[3:4],fixedPars[2],
                    modelInfo$maxBound[5:6],fixedPars[3],
                    modelInfo$maxBound[7:8],fixedPars[4])
    modelInfo$minBound=make_minBound
    modelInfo$maxBound=make_maxBound
  }
  #USER SPECIFIED BOUNDS CASE
  if((!is.null(minUserBound)) & (!is.null(maxUserBound))){  #
    #Check if the correct number of fixed parameters is supplied
    if((length(minUserBound)!=modelInfo$npars)|(length(maxUserBound)!=modelInfo$npars)){
      dummy=paste0("Error: Incorrect length of supplied bounds ", modelInfo$npars, " needed")
      logfile(dummy,file)
      logfile("Program terminated",file)
      stop(dummy)
    } else{
      modelInfo$minBound=minUserBound
      modelInfo$maxBound=maxUserBound
    }
    
  }
  
  #Otherwise do nothing
  return(modelInfo)
}

# Anjana : function for ar1 parameter calculation (moved up from modelSequencer.R)
#-----------------------------------------------------------------------

#Culley 2019 new loop to add ar(1)
add_ar1Param <- function(modelTag, modelInfo, datInd) {

  for(mod in 1:length(modelTag)){
    if(modelTag[mod]=="P-har-wgen"){
      # hard coded parameters
      ar1ParMult=0.001 # correlation between MULTIPLIER of monthly toals (not same as correlation between monthly totals) was 0.97
      multRange=0.8 # i.e. 0.1 is +/10%, so multiplier 95% limit is 0.9 to 1.1
      # translated param values needed for AR1
      multiplierMean=1
      sdJumpDistr=multRange/1.96*sqrt(1-ar1ParMult^2)
      # simulation
      X=arima.sim(n = datInd[[modelTag[mod]]]$nyr*12, list(ar =ar1ParMult),sd = sdJumpDistr)
      X=X@.Data+multiplierMean # extract data and add on mean
      multSim=rep(NA,datInd[[modelTag[mod]]]$ndays)
      for(iy in 1:datInd[[modelTag[mod]]]$nyr){
        for(im in 1:12){
          ind=datInd[[modelTag[mod]]]$i.yy[[iy]][which(datInd[[modelTag[mod]]]$i.yy[[iy]]%in%datInd[[modelTag[mod]]]$i.mm[[im]])]  # to save time this step can be pre-processed into datInd$i.yymm, a list of length 12*nyr
          multSim[ind]=X[(iy-1)*12+im]
        }
      }
      multSim<-pmax(multSim,0)
      modelInfo[[modelTag[mod]]]$ar1=multSim
    } else {
      if (modelTag[mod]!="Simple-ann") {
        modelInfo[[modelTag[mod]]]$ar1=NULL
      }
    }
  }
  return(modelInfo)
}
#------------------------------------------------------------------------
# #EXAMPLE 1
# modelTag="P-har12-wgen-FS"
# modelInfo=get.model.info(modelTag)
# update.model.info(modelTag=modelTag, modelInfo=modelInfo,fixedPars=c(1,2,3,4),minUserBound=NULL,maxUserBound=NULL)
# #EXAMPLE 2
# modelTag="P-ann-wgen"
# modelInfo=get.model.info(modelTag)
# update.model.info(modelTag=modelTag, modelInfo=modelInfo,fixedPars=NULL,minUserBound=c(0,0,3,4),maxUserBound=c(1,1,5,5))
#EXAMPLE 3
# modelTag="P-ann-wgen"
# modelInfo=get.model.info(modelTag)
# update.model.info(modelTag=modelTag, modelInfo=modelInfo,fixedPars=NULL,minUserBound=NULL,maxUserBound=NULL)

#RETURN VARIOUS PARS
return.simPriority<-function(modelInfo=NULL){
  return(modelInfo$simPriority)
}
return.simVar<-function(modelInfo=NULL){
  return(modelInfo$simVar)
}


parManager <- function(parS, modelEnv) {
  UseMethod("parManager", parS)
}

#PARAMETER MANAGER FOR WGEN STYLE RAIN SIMULATOR
  #NPERIOD, I.PP , DATIND INFO
  #IF NEEDED (E.G. HARMONIC) 
parManager.wgen <- function(parS = NULL,        # pars to split
                            modelEnv = NULL    # modelEnv that stores modelInfo, modelTag & datInd
  ){
  
  modelInfo <- modelEnv$modelInfo
  datInd <- modelEnv$datInd
  
  if (length(parS) != modelInfo$npars) {
    stop("Error: The number of parameters passed to the par manager does not match the number of parameters of the selected model")
  }
  
  #IF NO HARMONIC OR PARAMETER FIXING APPLIED IN MODEL VERSION
  if(is.na(modelInfo$ncycle) & is.na(modelInfo$fixedPars)){
    
    if(modelInfo$nperiod == 1) {
      #stop("Error: nperiod is not equal to 1, rep to create parameters of length `ndays` in the par manager will not work")
    #check length(pars == modelInfo$npars)   # if it fails put in a warning
    pdd=rep_len(parS[1:modelInfo$nperiod], datInd$ndays)                                 # extract first set of pars (pars evenly split across vector)
    pwd=rep_len(parS[(modelInfo$nperiod+1):(2*modelInfo$nperiod)], datInd$ndays)
    alpha=rep_len(parS[(2*modelInfo$nperiod+1):(3*modelInfo$nperiod)], datInd$ndays)
    beta=rep_len(parS[(3*modelInfo$nperiod+1):(4*modelInfo$nperiod)], datInd$ndays)
    } else if(modelInfo$nperiod == 4) {
      if(datInd$ndays != (length(unlist(datInd[["i.ss"]])))) stop("seasonal parameter assignment doesn't work properly.")
      pdd <- assignSeasPars(parS[1], parS[2], parS[3], parS[4], datInd[["i.ss"]])
      pwd <- assignSeasPars(parS[5], parS[6], parS[7], parS[8], datInd[["i.ss"]])
      alpha <- assignSeasPars(parS[9], parS[10], parS[11], parS[12], datInd[["i.ss"]])
      beta <- assignSeasPars(parS[13], parS[14], parS[15], parS[16], datInd[["i.ss"]])
    }
  }
  
  #if harmonics are required and no pars fixed - fit them
   if(!is.na(modelInfo$ncycle) & is.na(modelInfo$fixedPars)){
     
     
     pdd<-vector(mode="numeric",datInd$ndays) #Initialise vectors
     pwd<-vector(mode="numeric",datInd$ndays)
     alpha<-vector(mode="numeric",datInd$ndays)
     beta<-vector(mode="numeric",datInd$ndays)
      
     #Culley 2019 the below code includes leap years
     # for (i in 1:datInd$nyr){                 #There are probably more speedups here, but this is my leap year solution. It moves through each year, and either takes 365 or 366 points from a harmonic.
     #   pdd[datInd$i.yy[[i]]]<-harmonicFunc(x=seq(1:length(datInd$i.yy[[i]])),mean=parS[1],amp=parS[2],phase.ang = parS[3],k=1,nperiod=length(datInd$i.yy[[i]]))
     #   pwd[datInd$i.yy[[i]]]<-harmonicFunc(x=seq(1:length(datInd$i.yy[[i]])),mean=parS[4],amp=parS[5],phase.ang = parS[6],k=1,nperiod=length(datInd$i.yy[[i]]))
     #   alpha[datInd$i.yy[[i]]]<-harmonicFunc(x=seq(1:length(datInd$i.yy[[i]])),mean=parS[7],amp=parS[8],phase.ang = parS[9],k=1,nperiod=length(datInd$i.yy[[i]]))
     #   beta[datInd$i.yy[[i]]]<-harmonicFunc(x=seq(1:length(datInd$i.yy[[i]])),mean=parS[10],amp=parS[11],phase.ang = parS[12],k=1,nperiod=length(datInd$i.yy[[i]]))
     # }
     
     #Culley 2019 these parameter generators ignore leap years, so for long time series will become out of sync.
       pdd<-harmonicFunc(x=seq(1:datInd$ndays),mean=parS[1],amp=parS[2],phase.ang = parS[3],k=1,nperiod=365)
       pwd<-harmonicFunc(x=seq(1:datInd$ndays),mean=parS[4],amp=parS[5],phase.ang = parS[6],k=1,nperiod=365)
       alpha<-harmonicFunc(x=seq(1:datInd$ndays),mean=parS[7],amp=parS[8],phase.ang = parS[9],k=1,nperiod=365)
       beta<-harmonicFunc(x=seq(1:datInd$ndays),mean=parS[10],amp=parS[11],phase.ang = parS[12],k=1,nperiod=365)
     

     #Culley 2019 setting 0-1 limits for pdd,pwd.
     pdd<-pmin(pdd,1)
     pdd<-pmax(pdd,0)

     pwd<-pmin(pwd,1)
     pwd<-pmax(pwd,0)
     
     #Culley 2019 non negative limits for alpha,beta
     #alpha<-pmax(alpha,0)
     #beta<-pmax(beta,0)
     alpha<-pmax(alpha,.Machine$double.xmin) # Modified by SW 1/1/2020 as gamma distribution can't handle zero values for shape and scale
     beta<-pmax(beta,.Machine$double.xmin)       
    
   }
  

   #out is to - CALCULATE PAR VECTORS (PDD,PWD,ALPA,BETA)
  out=list(pdd=pdd,
           pwd=pwd,
           alpha=alpha,
           beta=beta)
   return(out)
}



assignSeasPars <- function(par1, par2, par3, par4, seasInd) {
  ndays <- length(unlist(seasInd))
  parAllDays <- rep_len(NA, length.out = ndays)
  parAllDays[seasInd[[1]]] <- par1
  parAllDays[seasInd[[2]]] <- par2
  parAllDays[seasInd[[3]]] <- par3
  parAllDays[seasInd[[4]]] <- par4
  if (any(is.na(parAllDays))) stop("NA values in assigned seasonal parameters.")
  return(parAllDays)
}


# Parameter manager for latent model WGEN
#-------------------------------------------------------------------------------------------------------------

parManager.latent <- function(parS = NULL,          # pars to split
                              modelEnv = NULL       # modelEnv that stores modelInfo, modelTag & datInd
) {
  
  modelInfo <- modelEnv$modelInfo
  datInd <- modelEnv$datInd
  
  if (length(parS) != modelInfo$npars) {
    stop("Error: The number of parameters passed to the par manager does not match the number of parameters of the selected model")
  }
  
  # IF NO HARMONIC OR PARAMETER FIXING APPLIED IN MODEL VERSION
  if(is.na(modelInfo$ncycle) & is.na(modelInfo$fixedPars)){
    
    if(modelInfo$nperiod != 1) {
      stop("Error: nperiod is not equal to 1, rep to create parameters of length `ndays` in the par manager will not work")
    }
    
    alpha <- rep_len(parS[1:modelInfo$nperiod], datInd$ndays)                                 # extract first set of pars (pars evenly split across vector)
    sigma <- rep_len(parS[(modelInfo$nperiod+1):(2*modelInfo$nperiod)], datInd$ndays)
    mu <- rep_len(parS[(2*modelInfo$nperiod+1):(3*modelInfo$nperiod)], datInd$ndays)
    lambda <- rep_len(parS[(3*modelInfo$nperiod+1):(4*modelInfo$nperiod)], datInd$ndays)
  }
  
  # if harmonics are required and no pars fixed - fit them
  if(!is.na(modelInfo$ncycle) & is.na(modelInfo$fixedPars)){
    
    alpha <- vector(mode = "numeric", datInd$ndays) # Initialise vectors - required? check if it makes a difference for computational time
    sigma <- vector(mode = "numeric", datInd$ndays)
    mu <- vector(mode = "numeric", datInd$ndays)
    lambda <- vector(mode = "numeric", datInd$ndays)
    
    # Culley 2019 these parameter generators ignore leap years, so for long time series will become out of sync.
    alpha <- harmonicFunc(x = seq(1:datInd$ndays), mean=parS[1], amp=parS[2], phase.ang = parS[3], k = 1, nperiod = 365)
    sigma <- harmonicFunc(x = seq(1:datInd$ndays), mean=parS[4], amp=parS[5], phase.ang = parS[6], k = 1, nperiod = 365)
    mu <- harmonicFunc(x = seq(1:datInd$ndays), mean=parS[7], amp=parS[8], phase.ang = parS[9], k = 1, nperiod = 365)
    lambda <- harmonicFunc(x = seq(1:datInd$ndays), mean=parS[10], amp=parS[11], phase.ang = parS[12], k = 1, nperiod = 365)
    
    # set limits on parameters after get point values from the harmonic function
    alpha <- pmax(alpha, -1)
    alpha <- pmin(alpha, 1)
    sigma <- pmax(sigma, 0)
    
  }
  
  out=list(alpha = alpha,
           sigma = sigma,
           mu = mu,
           lambda = lambda)
  return(out)
}

#-------------------------------------------------------------------------------------------------------------


# init.calib<- function(modelTag=NULL,  #model identifier
#                       modelInfo=NULL, #model information based on model identifier
#                       data=NULL,      #observed data (data frame) to be used in calibration
#                       datInd=NULL     #date index information
#   ){
# 
#   #FOR ALL MODELS
#     pdd=rep(0,modelInfo$nperiod); alpha=beta=pwd=pdd    #make space to store fitted pars
#     for(p in 1:modelInfo$nperiod){
#       #fit pwd and pdd
#       tmp=pdd.pwd.estimator(dat=data,ind=datInd$i.pp[[p]],threshold=0.00)
#       pdd[p]=tmp$pdd; pwd[p]=tmp$pwd
# 
#       #fit gamma pars (alpha & beta - using wgen manual labelling convention)
#       tmp=gammaParsMLE2(dat=data[datInd$i.pp[[p]]],wetThresh=0.00)
#       alpha[p]=tmp$shape; beta[p]=tmp$scale
#     }
#   
#     
#   
#   # MODELS WHERE A HARMONIC IS USED
#     #EACH PAR (PDD, PWD, ALPHA, BETA) IS FIT TO A HARMONIC (SOME EXCEPTIONS/SPECIAL CASES)
#     if(!is.na(modelInfo$ncycle)){
#       #fit pwd, pdd, alpha, beta
#       pdd.fit=fit.harmonic.opts(v.stat=pdd,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)  #from harmonicFit.R
#       pwd.fit=fit.harmonic.opts(v.stat=pwd,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
#       alpha.fit=fit.harmonic.opts(v.stat=alpha,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
#       beta.fit=fit.harmonic.opts(v.stat=beta,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
#       #PARS MUST BE ARANGED IN ORDER (MEAN, AMP, PHASE ANGLE) FOR EACH FITTED PAR
#       initCalibPars=c(unlist(pdd.fit,use.names = FALSE),
#                       unlist(pwd.fit,use.names = FALSE),
#                       unlist(alpha.fit,use.names = FALSE),
#                       unlist(beta.fit,use.names = FALSE))
#     }else{
#       #MAKE PAR VECTOR C(PAR1 X NPERIOD),(PAR2 X NPERIOD),...)
#       initCalibPars=c(pdd,pwd,alpha,beta)
#     }
#     
#     #PLOT UP
#     # windows();par(mfrow=c(2,2))
#     # plot(x=seq(1,modelInfo$nperiod),pdd,pch=16,xlab="period")
#     # cycle.pdd=harmonicFunc(x=seq(1,modelInfo$nperiod),mean=pdd.fit$mean,amp=pdd.fit$amp,phase.ang=pdd.fit$phase.ang,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
#     # lines(x=seq(1,modelInfo$nperiod),cycle.pdd,col="red")
#     # title("pdd")
#     # plot(x=seq(1,modelInfo$nperiod),pwd,pch=16,xlab="period")
#     # cycle.pwd=harmonicFunc(x=seq(1,modelInfo$nperiod),mean=pwd.fit$mean,amp=pwd.fit$amp,phase.ang=pwd.fit$phase.ang,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
#     # lines(x=seq(1,modelInfo$nperiod),cycle.pwd,col="red")
#     # title("pwd")
#     # plot(x=seq(1,modelInfo$nperiod),alpha,pch=16,xlab="period")
#     # alpha.cycle=harmonicFunc(x=seq(1,modelInfo$nperiod),mean=alpha.fit$mean,amp=alpha.fit$amp,phase.ang=alpha.fit$phase.ang,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
#     # lines(x=seq(1,modelInfo$nperiod),alpha.cycle,col="red")
#     # title("alpha")
#     # plot(x=seq(1,modelInfo$nperiod),beta,pch=16,xlab="period")
#     # beta.cycle=harmonicFunc(x=seq(1,modelInfo$nperiod),mean=beta.fit$mean,amp=beta.fit$amp,phase.ang=beta.fit$phase.ang,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
#     # lines(x=seq(1,modelInfo$nperiod),beta.cycle,col="red")
#     # title("beta")
#     
#     
#     #return pars from initial calibration
#     return(initCalibPars)
# }
# 
# #------------------------------------------------------------------------------------------
# gammaParsMLE2<- function(dat=NULL, #vector timeseries of rainfall amounts
#                          wetThresh=0.01, #threshold at which day is deemed wet
#                          ...){
#   
#   x=dat[which(dat>wetThresh)]    #get wet day amounts
#   nw=length(x)                   #no. wet days
#   x.mean=mean(x)                 #arithmetic mean of wet days
#  
#   s=log(x.mean)-sum(log(x))/nw
#   #est.shape - note in wgen manual shape is denoted using alpha
#   # ML method to estimate the 2 parameters of the gamma distribution from
#   # Wiki - https://en.wikipedia.org/wiki/Gamma_distribution#Characterization_using_shape_.CE.B1_and_rate_.CE.B2
#  # k.est=(3-s+sqrt((s-3)^2+24*s))/(12*s)
#   #estimator from wgen manual
#   Anum=8.898919+9.05995*s+0.9775373*s^2.0
#   Adom=s*(17.79728+11.968477*s+s^2.0)
#   k.est=Anum/Adom
#   if(k.est >=1.0){k.est=0.998}
#   
#   #est.scale - Note in wgen manual scale is denoted by beta
#   theta.est=x.mean/k.est
#   
#   out=list(scale=theta.est,shape=k.est)
#   return(out)
#   
# }
# 
# #---------------------------------------------------------------------------
# #Estimate pdd and pwd
# pdd.pwd.estimator<-function(dat=NULL,       # vector of rainfall values
#                             ind=NULL,       # indexes of days to assess
#                             threshold=0.01  # wet threshold
# ){
#   n=length(ind)
#   nw=length(which(dat[ind]>threshold))
#   nd=n-nw
#   #pDry=nd/n
#   
#   ind.prior=ind[-n]; ind.now=ind[-n]+1  #for clarity spell out which is the prior day series and current day series (DROP LAST VALUE TO AVOID ARRAY OVERFLOW)
#   
#   ind.wd.p=ind[which((dat[ind.prior])>threshold & (dat[ind.now])<=threshold)]  # GET INDICES OF WET(i-1) - DRY(i) PAIRS
#   n.wd=length(ind.wd.p)
#   
#   ind.dw.p=ind[which((dat[ind.prior])<=threshold & (dat[ind.now])>threshold)]  # GET INDICES OF DRY(i-1) - WET(i) PAIRS
#   n.dw=length(ind.dw.p)
#   
#   probs=list(pwd=(n.wd/nw),pdd=(1-(n.dw/nd)))
#   return(probs)
# }


#parMangement for TS generation
simHarTS.parmanager<-function(parS=NULL,   #par vector to be divided up (cors, wet (mean par, sd pars), dry (mean pars,sd pars))
                              modelTag=NULL,
                              modelInfo=NULL,
                              initCalibPars=NULL
){
  
  #position calculator function for harmonic pars
  posCalc<-function(npos,nAssocSeries,ncycle){st.pos=(npos-1+nAssocSeries)*(2*ncycle+1)+2}
  
  Hmean=list(); Hsd=list()
  
  #if harmonics are required, no pars fixed & conditional on wd status - fit them
  if(!is.na(modelInfo$ncycle) & is.na(modelInfo$fixedPars)){
    if(modelInfo$nAssocSeries ==0){  #If series residual generated alone
      cor1=parS[1]
      cor0=1
    }else{
      #INSERT WARNING
      #put on to-do list
      
      stop("Missing capcity to generate correl residuals")
    }
    
    if(modelInfo$WDcondition == TRUE){
      if(modelInfo$wdCycle == "sCycle"){
        #MAKE H MEANS WET AND DRY THE SAME
        st.parset=posCalc(npos=1,nAssocSeries=modelInfo$nAssocSeries,ncycle=modelInfo$ncycle)     # start position in vector for this parameter set
        Hmean$W=harmonicFunc(x=seq(1,modelInfo$nperiod),
                             mean=parS[st.parset],
                             amp=parS[(st.parset+1):(st.parset+modelInfo$ncycle)],
                             phase.ang=parS[(st.parset+modelInfo$ncycle+1):((st.parset+2*modelInfo$ncycle))],
                             k=modelInfo$ncycle,
                             nperiod=modelInfo$nperiod)
        
        Hmean$D=harmonicFunc(x=seq(1,modelInfo$nperiod),
                             mean=parS[st.parset],
                             amp=parS[(st.parset+1):(st.parset+modelInfo$ncycle)],
                             phase.ang=parS[(st.parset+modelInfo$ncycle+1):((st.parset+2*modelInfo$ncycle))],
                             k=modelInfo$ncycle,
                             nperiod=modelInfo$nperiod)
        
        #MAKE H SD WET AND DRY CONDITIONAL
        st.parset=posCalc(npos=2,nAssocSeries=modelInfo$nAssocSeries,ncycle=modelInfo$ncycle)         # start position in vector for this parameter set
        Hsd$W=harmonicFunc(x=seq(1,modelInfo$nperiod),
                           mean=parS[st.parset],
                           amp=parS[(st.parset+1):(st.parset+modelInfo$ncycle)],
                           phase.ang=parS[(st.parset+modelInfo$ncycle+1):((st.parset+2*modelInfo$ncycle))],
                           k=modelInfo$ncycle,
                           nperiod=modelInfo$nperiod)
        
        st.parset=posCalc(npos=3,nAssocSeries=modelInfo$nAssocSeries,ncycle=modelInfo$ncycle)          # start position in vector for this parameter set
        Hsd$D=harmonicFunc(x=seq(1,modelInfo$nperiod),
                           mean=parS[st.parset],
                           amp=parS[(st.parset+1):(st.parset+modelInfo$ncycle)],
                           phase.ang=parS[(st.parset+modelInfo$ncycle+1):((st.parset+2*modelInfo$ncycle))],
                           k=modelInfo$ncycle,
                           nperiod=modelInfo$nperiod)
        
      }else{
        st.parset=posCalc(npos=1,nAssocSeries=modelInfo$nAssocSeries,ncycle=modelInfo$ncycle)     # start position in vector for this parameter set
        Hmean$W=harmonicFunc(x=seq(1,modelInfo$nperiod),
                             mean=parS[st.parset],
                             amp=parS[(st.parset+1):(st.parset+modelInfo$ncycle)],
                             phase.ang=parS[(st.parset+modelInfo$ncycle+1):((st.parset+2*modelInfo$ncycle))],
                             k=modelInfo$ncycle,
                             nperiod=modelInfo$nperiod)
        
        st.parset=posCalc(npos=2,nAssocSeries=modelInfo$nAssocSeries,ncycle=modelInfo$ncycle)         # start position in vector for this parameter set
        Hsd$W=harmonicFunc(x=seq(1,modelInfo$nperiod),
                           mean=parS[st.parset],
                           amp=parS[(st.parset+1):(st.parset+modelInfo$ncycle)],
                           phase.ang=parS[(st.parset+modelInfo$ncycle+1):((st.parset+2*modelInfo$ncycle))],
                           k=modelInfo$ncycle,
                           nperiod=modelInfo$nperiod)
        
        st.parset=posCalc(npos=3,nAssocSeries=modelInfo$nAssocSeries,ncycle=modelInfo$ncycle)         # start position in vector for this parameter set
        Hmean$D=harmonicFunc(x=seq(1,modelInfo$nperiod),
                             mean=parS[st.parset],
                             amp=parS[(st.parset+1):(st.parset+modelInfo$ncycle)],
                             phase.ang=parS[(st.parset+modelInfo$ncycle+1):((st.parset+2*modelInfo$ncycle))],
                             k=modelInfo$ncycle,
                             nperiod=modelInfo$nperiod)
        
        st.parset=posCalc(npos=4,nAssocSeries=modelInfo$nAssocSeries,ncycle=modelInfo$ncycle)          # start position in vector for this parameter set
        Hsd$D=harmonicFunc(x=seq(1,modelInfo$nperiod),
                           mean=parS[st.parset],
                           amp=parS[(st.parset+1):(st.parset+modelInfo$ncycle)],
                           phase.ang=parS[(st.parset+modelInfo$ncycle+1):((st.parset+2*modelInfo$ncycle))],
                           k=modelInfo$ncycle,
                           nperiod=modelInfo$nperiod)
      }
      
    }else{
      #NOT CONDITIONAL ON WET-DRY
      st.parset=posCalc(npos=1,nAssocSeries=modelInfo$nAssocSeries,ncycle=modelInfo$ncycle)     # start position in vector for this parameter set
      Hmean$WD=harmonicFunc(x=seq(1,modelInfo$nperiod),
                            mean=parS[st.parset],
                            amp=parS[(st.parset+1):(st.parset+modelInfo$ncycle)],
                            phase.ang=parS[(st.parset+modelInfo$ncycle+1):((st.parset+2*modelInfo$ncycle))],
                            k=modelInfo$ncycle,
                            nperiod=modelInfo$nperiod)
      
      st.parset=posCalc(npos=2,nAssocSeries=modelInfo$nAssocSeries,ncycle=modelInfo$ncycle)         # start position in vector for this parameter set
      Hsd$WD=harmonicFunc(x=seq(1,modelInfo$nperiod),
                          mean=parS[st.parset],
                          amp=parS[(st.parset+1):(st.parset+modelInfo$ncycle)],
                          phase.ang=parS[(st.parset+modelInfo$ncycle+1):((st.parset+2*modelInfo$ncycle))],
                          k=modelInfo$ncycle,
                          nperiod=modelInfo$nperiod)
    }
  }else{
    stop("Fixed parameter functionality yet to come")
  }
  
  
  #CONDITIONS WHERE SOME PARS ARE FIXED
  #FIXED SEASONALITY OR OCCURANCE OR ...
  #INITICALIBPARS
  #need sytem for each type of fixing
  # if modelInfo$fixedPars="phase.angle"...
  #if(!is.na(modelInfo$ncycle) & (modelInfo$fixedPars == "phase.angle")){}
  # if modelInfo$fixedPars=
  
  #out is to - CALCULATE PAR VECTORS (corr_1,corr_0,Hmean,Hsd)
  out=list(cor1=cor1,
           cor0=cor0,
           Hmean=Hmean,
           Hsd=Hsd
  )
  return(out)
  
}

#tester
#simHarTS.parmanager(parS=seq(1,13),modelTag="Temp-har26-wgen",modelInfo=modelInfo,initCalibPars=NULL)

#WHICHPARS
whichPars<-function(simVar=NULL,
                    modelInfo=NULL
){
  parLoc=list()
  pos.start=1
  for(i in 1:length(simVar)){
    dummyA=modelInfo[[i]]$npars
    pos.end=(dummyA-1) + pos.start
    tmp=c(pos.start,pos.end)
    parLoc[[i]]=tmp         # store in list
    pos.start=pos.end + 1     # update pos.start ready for next model
  }
  return(parLoc)
}


#Update modelTag order
update.simPriority<-function(modelInfo=NULL){
  simPriority=sort(sapply(X=modelInfo,FUN=return.simPriority,USE.NAMES=TRUE)) #get simulation priority of each model
  modelTag=names(simPriority)                                                 # Force simVar="P" to come first via sorting by $simPrority
  return(modelTag)
}



