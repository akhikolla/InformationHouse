#######################################
##     WGEN FUNCTION LIBRARY         ##
#######################################

# CONTAINS

  # switch_simulator()
  #----------------------------
  # P_WGEN_master()
    # P_WGEN() - simulates rainfall using Richardson type model
        #Pstatus_WGEN()
        #Pamount_WGEN()
  #----------------------------
  # TS_WGEN_master()
    # TS_WGEN() -
      # calcDaySeries()
      # calcDayFunc()
      # residualGenerator()
      #****NOT DONE YET Also other correlated series models (via residuals etc) ---

#FUNCTIONS
#-----------------------------------------------------------------------------------------------------------

switch_simulator<-function(type=NULL,          # what vartype is being simulated
                           parS=NULL,
                           modelEnv = NULL,
                           randomVector = NULL,
                           wdSeries=NULL,      
                           resid_ts=NULL,
                           seed=NULL
                           ){

  switch(type,
         "P" = { 
           modelTag = modelEnv$P_modelEnv$modelTag
           switch(strsplit(modelTag, split="-")[[1]][3],
                        "latent" = { P_latent_master(parS=parS,
                                                     modelEnv = modelEnv$P_modelEnv,
                                                     randomVector = randomVector)
                          },
                        {
                          P_WGEN_master(parS=parS,              #RAIN SELECTED
                               modelEnv = modelEnv$P_modelEnv, 
                               randomVector = randomVector,
                               N=seed)
                        } )
               },
         "Temp" = { TS_WGEN_master(parS=parS,         
                                   modelTag=modelEnv$Temp_modelEnv$modelTag,      
                                   modelInfo=modelEnv$Temp_modelEnv$modelInfo,
                                   datInd=modelEnv$Temp_modelEnv$datInd,  
                                   randomVector = randomVector,
                                   initCalibPars=NULL, 
                                   wdSeries=wdSeries,      
                                   resid_ts=resid_ts,
                                   seed=seed)
           
                   },
         "PET" = { TS_WGEN_master(parS=parS,         
                                  modelTag=modelEnv$PET_modelEnv$modelTag,      
                                  modelInfo=modelEnv$PET_modelEnv$modelInfo,
                                  datInd=modelEnv$PET_modelEnv$datInd,   
                                  randomVector = randomVector,
                                  initCalibPars=NULL, 
                                  wdSeries=wdSeries,      
                                  resid_ts=resid_ts,
                                  seed=seed,
                                  trunc=0)
           
                 },
         "Radn" = { TS_WGEN_master(parS=parS,         
                                  modelTag=modelEnv$Radn_modelEnv$modelTag,      
                                  modelInfo=modelEnv$Radn_modelEnv$modelInfo,
                                  datInd=modelEnv$Radn_modelEnv$datInd,    
                                  randomVector = randomVector,
                                  initCalibPars=NULL, 
                                  wdSeries=wdSeries,      
                                  resid_ts=resid_ts,
                                  seed=seed,
                                  trunc=0)
           
         },
     

         -99.00
  )
}

P_WGEN_master<- function(parS,               # vector of pars (will change in optim)
                         modelEnv,
                         randomVector,
                         N=NULL             # seeds
){
  # Converts supplied pars into required format (e.g. if harmonic applied)
  class(parS) <- "wgen"
  parTS <- parManager(parS = parS, modelEnv = modelEnv)
  
  # Simulate rainfall timeseries
  sim <- P_WGEN(parTS = parTS,
                N = N,
                modelEnv = modelEnv,
                randomVector = randomVector)
  
  return(sim)  #return simulated rainfall
}

P_WGEN<-function(parTS,    
                 N=NULL,    #seed
                 modelEnv,
                 randomVector = NULL
){
  # unpack WGEN parameters
  parPwd <- parTS$pwd   # vector of pars (length = ndays)
  parPdd <- parTS$pdd
  parAlpha <- parTS$alpha
  parBeta <- parTS$beta
  ar1 <- modelEnv$modelInfo$ar1  # Culley 2019 add ar(1) model multipliers

  # sim length  
  ndays <- length(randomVector)
  
  # sim occurrence
  simS=Pstatus_WGEN(parPwd=parPwd,    # vector of pars for pwd (length = ndays)
                    parPdd=parPdd,    # vector of pars for pdd (length = ndays)
                    ndays=ndays,      
                    randomVector = randomVector 
  )
  
  # sim amounts
  simP=Pamount_WGEN(parAlpha=parAlpha,         # vector of pars for alpha (length = ndays)
                    parBeta=parBeta,           # vector of pars for beta (length = ndays)
                    status_ts=simS,            # TS vector of wet/dry statuses-obtained from the output of 'wvar_gen_Pstatus'
                    N=N,                       # random seed
                    ndays=ndays                
  )
  
  
  # multiply to inflate standard deviation at monthly scale and introduce monthly correlation
  if(!is.null(ar1)){simP$sim=simP$sim*ar1}
  
  syntP <- list(sim=simP$sim,
                seed=N)
  return(syntP)
}

Pstatus_WGEN <- function(parPwd,    # vector of pars for pwd (length = nperiod) - The modified cpp code expects vector of length ndays (not nperiod)
                         parPdd,    # vector of pars for pdd (length = nperiod) - The modified cpp code expects vector of length ndays (not nperiod)
                         ndays,
                         randomVector
                                                              
  ){
  
  drywet_TS <- Pstatus_WGEN_cpp(parPwd, parPdd, randomVector, ndays)
  return(drywet_TS)
}

Pamount_WGEN <- function(parAlpha=NULL,         # vector of pars for alpha (length = nperiod) - alpha has to be of length ndays (as per the code)
                         parBeta=NULL,          # vector of pars for beta (length = nperiod)
                         status_ts=NULL,             # TS vector of wet/dry statuses-obtained from the output of 'wvar_gen_Pstatus'
                         N=NULL,                # random seeds
                         ndays=NULL
){ 
  
  set.seed(N)     # seed seed to fix input needed for rgamma  ---  this creates a challenge for passing in a vector of random numbers
  rain <- vector(mode="numeric", ndays)  #allocate time series
  wet.days<-which(status_ts==1)         #index wet occurance
  rain[wet.days] <- rgamma(parAlpha[wet.days],shape=parAlpha[wet.days],scale=parBeta[wet.days])   #sample rain amount
  syntP <- list(sim=rain,  #
                seed=N)       
  return(syntP)
}

# Latent variable rainfall model
#-----------------------------------------------------------------------------------------------------

P_latent_master <- function(parS,                 # vector of pars (will change in optim)
                            modelEnv,
                            randomVector = NULL
) {
  
  # Converts supplied pars into required format (e.g. if harmonic applied)
  class(parS) <- "latent"
  parTS <- parManager(parS = parS, modelEnv = modelEnv)
  
  # Simulate rainfall timeseries
  sim <- P_latent(parTS = parTS,                 # wgen parameters
                  randomVector = randomVector)   # random vector of length ndays
  
  return(sim)  #return simulated rainfall
  
}


P_latent <- function(parTS,                  # list of parameters
                     randomVector = NULL
) {
  
  # Unpack WGEN parameters
  parAlpha <- parTS$alpha
  parSigma <- parTS$sigma
  parMu <- parTS$mu
  parLambda <- parTS$lambda
  ndays <- length(randomVector)
  
  # Calculate latent variable - latentX
  epsilonT <- qnorm(randomVector, mean = 0, sd = parSigma)
  latentX <- latentX_calc_cpp(parAlpha, epsilonT, ndays)
  latentX <- latentX + parMu
  
  # Transform latentX to rainfall
  rain <- rep_len(0, ndays)  
  latentX_pos_ind <- which(latentX > 0)
  rain[latentX_pos_ind] <- latentX[latentX_pos_ind] ^ parLambda[latentX_pos_ind]
  
  syntP <- list(sim = rain)       
  return(syntP)
  
}

#-----------------------------------------------------------------------------------------------------

TS_WGEN_master<- function(parS=NULL,         # vector of pars (will change in optim)
                          modelTag=NULL,      # tag to link/identify model 
                          modelInfo=NULL,
                          datInd=NULL,        # dat ind
                          randomVector = NULL,
                          initCalibPars=NULL, # vector of pars from initial baseline calibration
                          wdSeries=NULL,      # rain  series
                          resid_ts=NULL,
                          seed=NULL,
                          trunc=NULL
){
  #Converts supplied pars into required format (e.g. if harmonic applied)
  par=simHarTS.parmanager(parS=parS,modelTag=modelTag,modelInfo=modelInfo,initCalibPars=initCalibPars)
  
  #SIMULATE REQUIRED TIMESERIES
  sim=TS_WGEN(parCor0=par$cor0,         # correl pars
              parCor1=par$cor1, 
              parHmean=par$Hmean,       # mean harmonic pars (for $WD or $W & $D)
              parHsd=par$Hsd,           # sd harmonic pars (for $WD or $W & $D)
              k=modelInfo$ncycle,
              nperiod=modelInfo$nperiod,
              ndays=datInd$ndays,
              nAssocSeries=modelInfo$nAssocSeries,
              i.pp=datInd$i.pp,
              initCalibPars=initCalibPars,  # initial model calibration pars
              WDcondition=modelInfo$WDcondition,   # generate ts conditional on wet-dry series
              wdSeries=wdSeries,        # wet-dry series
              resid_ts=resid_ts,        # leave capacity to generate or receive residuals
              seed=seed,                 # seed for residuals generation
              randomVector=randomVector
  )
  #Make data below threshold zero (i.e. no negative radiation)
  if(!is.null(trunc)){
   badIndx=which(sim$sim<trunc)
   if(length(badIndx)>0){sim$sim[badIndx]=0}
  }
  
  return(sim)
}

TS_WGEN<-function(parCor0=NULL,         # correl pars
                  parCor1=NULL,
                  parHmean=NULL,      # mean harmonic pars
                  parHsd=NULL,        # sd harmonic pars
                  k=NULL,
                  nperiod=NULL,
                  ndays=NULL,
                  nAssocSeries=NULL,
                  i.pp=NULL,
                  initCalibPars=NULL,  # initial model calibration pars
                  WDcondition=FALSE,   # generate ts conditional on wet-dry series
                  wdSeries=NULL,       # wet-dry series
                  resid_ts=NULL,       # leave capacity to generate or receive residuals
                  seed=NULL,            # seed for residuals generation
                  randomVector=NULL
){
  
  #generate residual series if none supplied
  if(is.null(resid_ts)){
    resid_ts=residualGenerator(parCor0=parCor0,
                               parCor1=parCor1,
                               ndays=ndays,
                               nAssocSeries=nAssocSeries,
                               randomVector=randomVector
    )
  }
  
  #divy up parameters and simulate
  sim=calcDaySeries(Hpar_m=parHmean,       #
                    Hpar_sd=parHsd,     
                    k=k,           
                    nperiod=nperiod,     
                    ndays=ndays,  
                    resid_ts=resid_ts,
                    i.pp=i.pp,         
                    WDcondition=WDcondition, 
                    wdSeries=wdSeries     
  )
  out=list(sim=sim,seed=seed)
  return(out)
}


#Generate the daily series
calcDaySeries<-function(Hpar_m=NULL,       #harmonic pars for means of length nperiod 
                        Hpar_sd=NULL,      #harmonic pars for std dev 
                        k=NULL,            #No. cycles in model specified
                        nperiod=NULL,      #No. period require by model specified
                        ndays=NULL,        # No. of days simulated
                        resid_ts=NULL,     #residual error timeseries
                        i.pp=NULL,         #period  indices
                        WDcondition=FALSE, #generate ts conditional on wet-dry series
                        wdSeries=NULL      #wetdry series
){
  
  genTS=rep(NA,ndays)           #MAKE BLANK VECTOR FOR TS
  if(WDcondition==FALSE){       #IF NO WET-DRY CONDITION
    #CALCULATE VALUE BASED ON PERIOD PARS & STUFF BACK INTO TS AT CORRECT POINT
    for(p in 1:nperiod){
      genTS[i.pp[[p]]]=calcDayFunc(mean=Hpar_m$WD[p],sd=Hpar_sd$WD[p],err=resid_ts[i.pp[[p]]])
    }
  }else{                        #IF WET-DRY CONDITIONAL
    #DETERMINE WET DAY INDICIES
    indW=which(wdSeries>0.01)
    #CALCULATE VALUE BASED ON PERIOD PARS & STUFF BACK INTO TS AT CORRECT POINT
    for(p in 1:nperiod){
      indWP=intersect(i.pp[[p]],indW)
      genTS[indWP]=calcDayFunc(mean=Hpar_m$W[p],sd=Hpar_sd$W[p],err=resid_ts[indWP])  #Store values on wet days for period
      indDP=outersect(i.pp[[p]],indWP)
      genTS[indDP]=calcDayFunc(mean=Hpar_m$D[p],sd=Hpar_sd$D[p],err=resid_ts[indDP])  #Store values on dry days for period
    }
  }
  return(genTS)
}

#SIMULATE TEMPERATURE SERIES
calcDayFunc<-function(mean=NULL,
                      sd=NULL,
                      err=NULL
){
  val=mean+sd*err
}

#GENERATE RESIDUALS USING SUPPLIED PARAMETERS
residualGenerator<-function(parCor0=NULL,
                            parCor1=NULL,
                            ndays=NULL,
                            nAssocSeries=0,
                            randomVector=NULL
){
  
  #GENERATE RANDOM NUMBERS FROM A STANDARD NORMAL (MEAN=0, STD=1)
  
  if(nAssocSeries==0){  #if just one series
    RN_res<-qnorm(randomVector, mean=0, sd=1)
    res_gen <- residualGenerator_cpp(RN_res, parCor1)
    
  }else{
    #IMPLEMENT A/B MATRIX MATHS HERE
    #
    stop("functionality not yet implemented")
  }
  
  return(res_gen)
}

