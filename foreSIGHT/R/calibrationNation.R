############################
#### CALIBRATION NATION ####
############################

#CONTAINS
  #modCalibrator() - visible function - calibrates models based on modelTags
    #init.calib.master() - master calibration function
    #init.calib.Pwgen() - initial calibration function (rainfall). Uses obs data to fit model to baseline climate.
      #uses model tag to determine what should be fitted (e.g. how many periods, using harmonic, etc.)
      #gammaParsMLE2() - args(dat, wetThresh)
      #pdd.pwd.estimator() -args (dat,ind,threshold)
    # init.calib.harTS()  - fitter for harmonic (non rainfall models)
    # calcDayResidFunc()
    # calcDayResidFunc_wdSep()


#------------------------------------------------
modCalibrator<-function(obs=NULL,
                        modelTag=NULL,
                        window=NULL     #sets moving average window to calibrate daily gamma parameters for "P-har-WGEN"
){
  
  #check obs data, modelTag (abbridged version from control)
  argCheck=argument_check_calibrator(names=names(obs),  #Need an invalid modelTage Choice flag if FS type
                                     obs=obs,
                                     modelTag=modelTag)
  #Update inputs (e.g truncate incomplete years)
  obs=input_check_calibrator(obs=obs)
  
  #GET ADDITIONAL MODEL INFO, SIMVARS etc
  modelInfo=get.multi.model.info(modelTag=modelTag)
  modelTag=update.simPriority(modelInfo=modelInfo)
  simVar=sapply(X=modelInfo[modelTag],FUN=return.simVar,USE.NAMES=TRUE)       #?CREATE MODEL MASTER INFO - HIGHER LEVEL?
  
  #Get date information
  datInd=mod.get.date.ind(obs=obs[,c("year","month","day")],modelTag=modelTag,modelInfo=modelInfo) #Get datInd for all modelTags
  
  #calibrate each model
  parDat=list()
  for(i in 1:length(modelTag)){
    if(is.null(obs$P)){rain=NULL}else{rain=obs$P}  #little temp switch here
    
    tempVar=modelInfo[[modelTag[i]]]$simVar
    
    parOut=init.calib.master(modelTag=modelTag[i],
                             modelInfo=modelInfo[[modelTag[i]]],
                             data=obs[[tempVar]],
                             datInd=datInd[[modelTag[i]]],
                             rain=rain,  
                             threshold=0,
                             window=window)
    parDat[[modelTag[i]]]=parOut
  }

  #return parameters in list
  return(parDat)
}

#obs<-read.csv("C:\\Users\\Phil\\Dropbox\\Goyder\\test\\AliceSpringsAWAP.csv")
#modCalibrator(obs=obs,modelTag="P-har6-wgen")

init.calib.master<-function(modelTag=NULL,        #these are set to match the arguments for init.calib
                            modelInfo=NULL,
                            data=NULL,            #takes a vector of temperature
                            datInd=NULL,
                            rain=NULL,
                            threshold=0,
                            window=NULL
){
  
  #SWITCH BASED ON WHAT simVar/modelType
  
  switch(modelInfo$simVar,
         "All" = { #If simple scaling selected
                  stop("Invalid request - simple scaling model does not require calibration")
                  },   
         "P" = { #If a rainfall model is selected
                parSet=init.calib.Pwgen(modelTag=modelTag,modelInfo=modelInfo,data=data,datInd=datInd)
                  },     
               { #If non-rainfall harmonic model is selected
                parSet=init.calib.harTS(modelTag=modelTag,
                                        modelInfo=modelInfo,
                                        data=data,
                                        datInd=datInd,
                                        rain=rain,
                                        threshold=threshold)  
                }            
         )
  
  return(parSet)
}

init.calib.Pwgen<- function(modelTag=NULL,  #model identifier
                            modelInfo=NULL, #model information based on model identifier
                            data=NULL,      #observed data (data frame) to be used in calibration
                            datInd=NULL,     #date index information
                            window=15 #Culley 27/11/19 added argument for width of window to collect data for daily parameters
){
  
  #RECLASSIFY IF "P-har12-wgen-FS"
  if(modelTag=="P-har12-wgen-FS"){modelTag="P-har12-wgen"} #RECLASSIFY IF "P-har12-wgen-FS"
  
  #FOR ALL RAINFALL MODELS
  pdd=rep(0,modelInfo$nperiod); alpha=beta=pwd=pdd    #make space to store fitted pars
  
  for(p in 1:modelInfo$nperiod){
    
    #Culley 27/11/19 - Changing modCalibrate for "P-har-wgen"
    if(modelTag=="P-har-wgen"){ #This will still fit harmonics daily, (i.e. modelInfor$nperiod=365), but the indexing needs to include a window of data either side of the single day
      
      i.ww<-rep(NA,1)
      for(j in 1:length(datInd$i.pp[[p]])){       #extend each index (one point each year) to include window days either side
        datarun<-seq(datInd$i.pp[[p]][j]-window,datInd$i.pp[[p]][j]+window)
        if(j==1){datarun<-datarun[datarun>0]} #in first and final years, trim data, as sometimes cannot grab from left or right, respectively
        if(j==length(datInd$i.pp[[p]])){datarun<-datarun[datarun<=tail(datInd$i.pp[[modelInfo$nperiod]], n=1)]}
        i.ww<-append(i.ww,datarun)
      }
      i.ww<-i.ww[-(1)]
      
      tmp=pdd.pwd.estimator(dat=data,ind=i.ww,threshold=0.00)
      pdd[p]=tmp$pdd; pwd[p]=tmp$pwd
      
      #fit gamma pars (alpha & beta - using wgen manual labelling convention)
      tmp=gammaParsMLE2(dat=data[i.ww],wetThresh=0.00)
      alpha[p]=tmp$shape; beta[p]=tmp$scale
      
    } else {
    #fit pwd and pdd
      tmp=pdd.pwd.estimator(dat=data,ind=datInd$i.pp[[p]],threshold=0.00)
      pdd[p]=tmp$pdd; pwd[p]=tmp$pwd
      
      #fit gamma pars (alpha & beta - using wgen manual labelling convention)
      tmp=gammaParsMLE2(dat=data[datInd$i.pp[[p]]],wetThresh=0.00)
      alpha[p]=tmp$shape; beta[p]=tmp$scale
    }
  }
  
  # MODELS WHERE A HARMONIC IS USED
  #EACH PAR (PDD, PWD, ALPHA, BETA) IS FIT TO A HARMONIC (SOME EXCEPTIONS/SPECIAL CASES)
  if(!is.na(modelInfo$ncycle)){
    #fit pwd, pdd, alpha, beta
    pdd.fit=fit.harmonic.opts(v.stat=pdd,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)  #from harmonicFit.R
    pwd.fit=fit.harmonic.opts(v.stat=pwd,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
    alpha.fit=fit.harmonic.opts(v.stat=alpha,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
    beta.fit=fit.harmonic.opts(v.stat=beta,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
    #PARS MUST BE ARANGED IN ORDER (MEAN, AMP, PHASE ANGLE) FOR EACH FITTED PAR
    initCalibPars=c(unlist(pdd.fit,use.names = FALSE),
                    unlist(pwd.fit,use.names = FALSE),
                    unlist(alpha.fit,use.names = FALSE),
                    unlist(beta.fit,use.names = FALSE))
  }else{
    #MAKE PAR VECTOR C(PAR1 X NPERIOD),(PAR2 X NPERIOD),...)
    initCalibPars=c(pdd,pwd,alpha,beta)
  }
  
  #PLOT UP
  # windows();par(mfrow=c(2,2))
  # plot(x=seq(1,modelInfo$nperiod),pdd,pch=16,xlab="period")
  # cycle.pdd=harmonicFunc(x=seq(1,modelInfo$nperiod),mean=pdd.fit$mean,amp=pdd.fit$amp,phase.ang=pdd.fit$phase.ang,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
  # lines(x=seq(1,modelInfo$nperiod),cycle.pdd,col="red")
  # title("pdd")
  # plot(x=seq(1,modelInfo$nperiod),pwd,pch=16,xlab="period")
  # cycle.pwd=harmonicFunc(x=seq(1,modelInfo$nperiod),mean=pwd.fit$mean,amp=pwd.fit$amp,phase.ang=pwd.fit$phase.ang,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
  # lines(x=seq(1,modelInfo$nperiod),cycle.pwd,col="red")
  # title("pwd")
  # plot(x=seq(1,modelInfo$nperiod),alpha,pch=16,xlab="period")
  # alpha.cycle=harmonicFunc(x=seq(1,modelInfo$nperiod),mean=alpha.fit$mean,amp=alpha.fit$amp,phase.ang=alpha.fit$phase.ang,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
  # lines(x=seq(1,modelInfo$nperiod),alpha.cycle,col="red")
  # title("alpha")
  # plot(x=seq(1,modelInfo$nperiod),beta,pch=16,xlab="period")
  # beta.cycle=harmonicFunc(x=seq(1,modelInfo$nperiod),mean=beta.fit$mean,amp=beta.fit$amp,phase.ang=beta.fit$phase.ang,k=modelInfo$ncycle,nperiod=modelInfo$nperiod)
  # lines(x=seq(1,modelInfo$nperiod),beta.cycle,col="red")
  # title("beta")

  #return pars from initial calibration
  return(initCalibPars)
}

#------------------------------------------------------------------------------------------
gammaParsMLE2<- function(dat=NULL, #vector timeseries of rainfall amounts
                         wetThresh=0.01, #threshold at which day is deemed wet
                         ...){
  
  x=dat[which(dat>wetThresh)]    #get wet day amounts
  nw=length(x)                   #no. wet days
  x.mean=mean(x)                 #arithmetic mean of wet days
  
  s=log(x.mean)-sum(log(x))/nw
  #est.shape - note in wgen manual shape is denoted using alpha
  # ML method to estimate the 2 parameters of the gamma distribution from
  # Wiki - https://en.wikipedia.org/wiki/Gamma_distribution#Characterization_using_shape_.CE.B1_and_rate_.CE.B2
  # k.est=(3-s+sqrt((s-3)^2+24*s))/(12*s)
  #estimator from wgen manual
  Anum=8.898919+9.05995*s+0.9775373*s^2.0
  Adom=s*(17.79728+11.968477*s+s^2.0)
  k.est=Anum/Adom
  if(k.est >=1.0){k.est=0.998}
  
  #est.scale - Note in wgen manual scale is denoted by beta
  theta.est=x.mean/k.est
  
  out=list(scale=theta.est,shape=k.est)
  return(out)
  
}

#---------------------------------------------------------------------------
#Estimate pdd and pwd
pdd.pwd.estimator<-function(dat=NULL,       # vector of rainfall values
                            ind=NULL,       # indexes of days to assess
                            threshold=0.01  # wet threshold
){
  n=length(ind)
  nw=length(which(dat[ind]>threshold))
  nd=n-nw
  #pDry=nd/n
  
  ind.prior=ind[-n]; ind.now=ind[-n]+1  #for clarity spell out which is the prior day series and current day series (DROP LAST VALUE TO AVOID ARRAY OVERFLOW)
  
  ind.wd.p=ind[which((dat[ind.prior])>threshold & (dat[ind.now])<=threshold)]  # GET INDICES OF WET(i-1) - DRY(i) PAIRS
  n.wd=length(ind.wd.p)
  
  ind.dw.p=ind[which((dat[ind.prior])<=threshold & (dat[ind.now])>threshold)]  # GET INDICES OF DRY(i-1) - WET(i) PAIRS
  n.dw=length(ind.dw.p)
  
  probs=list(pwd=(n.wd/nw),pdd=(1-(n.dw/nd)))
  return(probs)
}

### !!!!!!!NOTE THIS IS SET UP FOR 1 CYCLE ONLY AT THE MOMENT!
init.calib.harTS<-function(modelTag=NULL,        #these are set to match the arguments for init.calib
                           modelInfo=NULL,
                           data=NULL,            #takes a vector of temperature
                           datInd=NULL,
                           rain=NULL,
                           threshold=0
){
  
    #First get harmonics for mean and stdev
    periodStats<-get.period.stats(data=data,
                                  nperiod=modelInfo$nperiod,
                                  i.hh=datInd$i.pp,
                                  sep.wd=modelInfo$WDcondition,
                                  # omit.drypars=FALSE,   #choice to not calculate pars on dry days
                                  rain=rain,            #vector of rain
                                  threshold=threshold)
    
    
    if(modelInfo$WDcondition==FALSE){ #IF NO WET-DRY SEPARATION
      #CALCULATE HARMONIC PARAMETERS FOR THE MEANS AND STANDARD DEVIATIONS
      par_meanHar=fit.harmonic.opts(nperiod = modelInfo$nperiod,v.stat=periodStats$all$m,k=modelInfo$ncycle)
      par_sdHar=fit.harmonic.opts(nperiod = modelInfo$nperiod,v.stat=periodStats$all$sd,k=modelInfo$ncycle)
      
      #CALCULATE PARAMETERS FOR EACH PERIOD BLOCK USING FITTED HARMONICS
      Hpar_m=harmonicFunc(x=seq(1,modelInfo$nperiod),
                          mean=par_meanHar$mean,
                          amp=par_meanHar$amp,
                          phase.ang=par_meanHar$phase.ang,
                          k=modelInfo$ncycle,
                          nperiod=modelInfo$nperiod)
      
      Hpar_sd=harmonicFunc(x=seq(1,modelInfo$nperiod),
                           mean=par_sdHar$mean,
                           amp=par_sdHar$amp,
                           phase.ang=par_sdHar$phase.ang,
                           k=modelInfo$ncycle,
                           nperiod=modelInfo$nperiod)
      
      #CALCULATE VALUE BASED ON PERIOD PARS & STUFF BACK INTO TS AT CORRECT POINT
      genRes=rep(NA,datInd$ndays)                                    #MAKE BLANK VECTOR FOR RESIDUAL TS
      for(p in 1:modelInfo$nperiod){
        genRes[datInd$i.pp[[p]]]=calcDayResidFunc(data=data[datInd$i.pp[[p]]],mean=Hpar_m[p],sd=Hpar_sd[p])
      }
      
      # CALCULATE LAG-1 CORRELATION IN RESIDUALS
      a=acf(genRes,na.action=na.pass,plot=F) # GET OBSERVED AT-SITE TEMPORAL AUTO-CORRELATION OF RESIDUALS
      r.res=a$acf[2]                         # STORE LAG-1 TEMPORAL AUTO-CORRELATION
      
      #COLLATE PARAMETER SET IN REQUIRED ORDER
      parSet<-as.numeric(c(r.res,unlist(par_meanHar),unlist(par_sdHar))) #cor1, mean.m, mean.amp,mean.phase, sd.m, sd.amp, sd.phase
      
      
    }else{   #IF WET-DRY SEPARATION
      
      #CALCULATE HARMONIC PARAMETERS FOR THE MEANS AND STANDARD DEVIATIONS
      #DRY
      par_meanHar_dry=fit.harmonic.opts(nperiod = modelInfo$nperiod,v.stat=periodStats$dry$m,k=modelInfo$ncycle)
      par_sdHar_dry=fit.harmonic.opts(nperiod = modelInfo$nperiod,v.stat=periodStats$dry$sd,k=modelInfo$ncycle)
      #WET
      par_meanHar_wet=fit.harmonic.opts(nperiod = modelInfo$nperiod,v.stat=periodStats$wet$m,k=modelInfo$ncycle)
      par_sdHar_wet=fit.harmonic.opts(nperiod = modelInfo$nperiod,v.stat=periodStats$wet$sd,k=modelInfo$ncycle)
      
      #CALCULATE PARAMETERS FOR EACH PERIOD BLOCK USING FITTED HARMONICS
      #DRY
      Hpar_m_dry=harmonicFunc(x=seq(1,modelInfo$nperiod),
                              mean=par_meanHar_dry$mean,
                              amp=par_meanHar_dry$amp,
                              phase.ang=par_meanHar_dry$phase.ang,
                              k=modelInfo$ncycle,
                              nperiod=modelInfo$nperiod)
      
      Hpar_sd_dry=harmonicFunc(x=seq(1,modelInfo$nperiod),
                               mean=par_sdHar_dry$mean,
                               amp=par_sdHar_dry$amp,
                               phase.ang=par_sdHar_dry$phase.ang,
                               k=modelInfo$ncycle,
                               nperiod=modelInfo$nperiod)
      #WET
      Hpar_m_wet=harmonicFunc(x=seq(1,modelInfo$nperiod),
                              mean=par_meanHar_wet$mean,
                              amp=par_meanHar_wet$amp,
                              phase.ang=par_meanHar_wet$phase.ang,
                              k=modelInfo$ncycle,
                              nperiod=modelInfo$nperiod)
      
      Hpar_sd_wet=harmonicFunc(x=seq(1,modelInfo$nperiod),
                               mean=par_sdHar_wet$mean,
                               amp=par_sdHar_wet$amp,
                               phase.ang=par_sdHar_wet$phase.ang,
                               k=modelInfo$ncycle,
                               nperiod=modelInfo$nperiod)
      
      #AND AN ODD CASE TOO FOR WHEN YOU ONLY VARY SD ON THE WET DRY SPLIT
      if(modelInfo$wdCycle=="sCycle"){
        periodStats_mCycle<-get.period.stats(data=data,nperiod=modelInfo$nperiod,i.hh=datInd$i.pp,sep.wd=FALSE)
        par_meanHar_all=fit.harmonic.opts(nperiod = modelInfo$nperiod,v.stat=periodStats_mCycle$all$m,k=modelInfo$ncycle)
        #CALCULATE PARAMETERS FOR EACH PERIOD BLOCK USING FITTED HARMONICS
        Hpar_m_all=harmonicFunc(x=seq(1,modelInfo$nperiod),
                                mean=par_meanHar_all$mean,
                                amp=par_meanHar_all$amp,
                                phase.ang=par_meanHar_all$phase.ang,
                                k=modelInfo$ncycle,
                                nperiod=modelInfo$nperiod)
      }
      
      
      #WHICH OF THE CYCLES ARE USED?
      switch(modelInfo$wdCycle,
             "All"={
               #CALCULATE VALUE BASED ON PERIOD PARS & WET DRY STATUS & STUFF BACK INTO TS AT CORRECT POINT
               genRes=rep(NA,datInd$ndays)                                    #MAKE BLANK VECTOR FOR RESIDUAL TS
               for(p in 1:modelInfo$nperiod){
                 genRes[datInd$i.pp[[p]]]=calcDayResidFunc_wdSep(data=data[datInd$i.pp[[p]]],
                                                                 mean_dry=Hpar_m_dry[p],
                                                                 sd_dry=Hpar_sd_dry[p],
                                                                 mean_wet=Hpar_m_wet[p],
                                                                 sd_wet=Hpar_sd_wet[p],
                                                                 rain=rain[datInd$i.pp[[p]]],
                                                                 threshold=threshold)
               }
               
               # CALCULATE LAG-1 CORRELATION IN RESIDUALS
               a=acf(genRes,na.action=na.pass,plot=F) # GET OBSERVED AT-SITE TEMPORAL AUTO-CORRELATION OF RESIDUALS
               r.res=a$acf[2]                         # STORE LAG-1 TEMPORAL AUTO-CORRELATION
               
               #COLLATE PARAMETER SET IN REQUIRED ORDER
               parSet<-as.numeric(c(r.res,
                                    unlist(par_meanHar_wet),
                                    unlist(par_sdHar_wet),
                                    unlist(par_meanHar_dry),
                                    unlist(par_sdHar_dry)
               )) 
               #("cor0","W-mCycle-m","W-mCycle-amp","W-mCycle-ang", "W-sCycle-m","W-sCycle-amp","W-sCycle-ang", "D-mCycle-m","D-mCycle-amp","D-mCycle-ang","D-sCycle-m","D-sCycle-amp","D-sCycle-ang")
               
               
             },
             "sCycle"={
               #THIS DOES NOT REALLY MAKE SENSE IN TERMS OF THIS FORWARD CALIBRATION
               #CALCULATE VALUE BASED ON PERIOD PARS & WET DRY STATUS & STUFF BACK INTO TS AT CORRECT POINT
               genRes=rep(NA,datInd$ndays)                                    #MAKE BLANK VECTOR FOR RESIDUAL TS
               for(p in 1:modelInfo$nperiod){
                 genRes[datInd$i.pp[[p]]]=calcDayResidFunc_wdSep(data=data[datInd$i.pp[[p]]],
                                                                 mean_dry=Hpar_m_all[p],
                                                                 sd_dry=Hpar_sd_dry[p],
                                                                 mean_wet=Hpar_m_all[p],
                                                                 sd_wet=Hpar_sd_wet[p],
                                                                 rain=rain[datInd$i.pp[[p]]],
                                                                 threshold=threshold)
               }
               # CALCULATE LAG-1 CORRELATION IN RESIDUALS
               a=acf(genRes,na.action=na.pass,plot=F) # GET OBSERVED AT-SITE TEMPORAL AUTO-CORRELATION OF RESIDUALS
               r.res=a$acf[2]                         # STORE LAG-1 TEMPORAL AUTO-CORRELATION
               
               #COLLATE PARAMETER SET IN REQUIRED ORDER
               parSet<-as.numeric(c(r.res,
                                    unlist(par_meanHar_all),
                                    unlist(par_sdHar_wet),
                                    unlist(par_sdHar_dry)
               )) 
             },
             {print("error")}
      ) #end switch
      
    } #end wet-dry separation
    

  
  return(parSet)
}

#CALCUATE RESIDUALS OF TIME SERIES
calcDayResidFunc<-function(data=NULL,
                           mean=NULL,
                           sd=NULL){
  err=(data-mean)/sd
}

#CALCULATE RESIDUALS BASED ON WET DRY STATUS
calcDayResidFunc_wdSep<-function(data=NULL,
                                 mean_dry=NULL,
                                 sd_dry=NULL,
                                 mean_wet=NULL,
                                 sd_wet=NULL,
                                 rain=NULL,
                                 threshold=NULL
){
  
  nLen=length(rain); err=rep(NA,nLen) #make err vector
  
  #CALCULATE WET ERRORS
  w.ind=which(rain>threshold)
  err[w.ind]=(data[w.ind]-mean_wet)/sd_wet
  #CALCULATE DRY ERROR
  d.ind=which(rain<=threshold)
  err[d.ind]=(data[d.ind]-mean_dry)/sd_dry
  
  #RETURN ERRORS
  err
}
