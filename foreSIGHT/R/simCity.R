########################
##      SIM CITY      ##
########################

#CONTAINS
  #modSimulator() -
  #argument_check_simulator() - inputs (modelTag=NULL,datStart=NULL,datFinish=NULL,parS=NULL)
#------------------------------------------------
#FUNCTIONS
modSimulator<-function(datStart=NULL,
                       datFinish=NULL,
                       modelTag=NULL,      
                       parS=NULL,      #list that matches modelTag names
                       seed=NULL,
                       file=NULL,
                       IOmode="suppress"
                        ){
  

  #DO CHECKS - need date checks and model tag checks, par checks
  argument_check_simulator(modelTag=modelTag,datStart=datStart,datFinish=datFinish,parS=parS)
  
  #GET ADDITIONAL MODEL INFO, SIMVARS etc
  modelInfo=get.multi.model.info(modelTag=modelTag)
  modelTag=update.simPriority(modelInfo=modelInfo)
  simVar=sapply(X=modelInfo[modelTag],FUN=return.simVar,USE.NAMES=TRUE)       #?CREATE MODEL MASTER INFO - HIGHER LEVEL?
  
  #Manage dates
  dates=makeDates(datStart=datStart,datFinish=datFinish)  #produces dates data.frame (year, month, day)
  datInd=mod.get.date.ind(obs=dates,modelTag=modelTag,modelInfo=modelInfo) #Get datInd for all modelTags
  
  #Culley 2019 new loop to add ar(1) 
  set.seed(seed)
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
      modelInfo[[modelTag[mod]]]$ar1=NULL
    }
  }
  #Culley 2019 new seed generator, but this takes longer
  # seeds<-list()
  # ndays<-datInd[[modelTag[1]]]$ndays
  # set.seed(seed)
  # seeds$wdstatus<-runif((ndays),0,1)
  # seeds$rainamount<-runif((ndays),0,1)
  # seeds$residuals<-rnorm(n=(ndays),mean=0,sd=1)
  
  #LOOP OVER EACH STOCHASTIC MODEL NOMINATED
  out=list()
  
  for(mod in 1:length(modelTag)){
    randomVector <- runif(n=datInd[[modelTag[mod]]]$ndays) # Random vector to be passed into weather generator to reduce runtime
    #IF CONDITIONED ON DRY-WET STATUS, populate wdStatus
    switch(simVar[mod], #
           "P" = {wdStatus=NULL},
           "Temp" = {if(modelInfo[[modelTag[mod]]]$WDcondition==TRUE){
                         wdStatus=out[["P"]]$sim
                       }else{
                         wdStatus=NULL   
                       }
                    },
           {wdStatus=NULL}  #default
    )

    #GRAB PARS RELATED TO modelTag
    parSel=as.double(parS[[modelTag[mod]]])
    
    # write data to model environment
    #----------------------------------
    write_model_env(envir = foreSIGHT_modelEnv, 
                    modelInfo = modelInfo[[modelTag[mod]]], 
                    modelTag = modelTag[mod], 
                    datInd = datInd[[modelTag[mod]]] 
    )
    #-----------------------------------

    out[[simVar[mod]]]=switch_simulator(type=modelInfo[[modelTag[mod]]]$simVar,
                                        parS=parSel,
                                        modelEnv = foreSIGHT_modelEnv,
                                        randomVector = randomVector,
                                        wdSeries=wdStatus,
                                        resid_ts=NULL,
                                        seed=seed)
  }  #end model loop
  
  simDat=makeOutputDataframe(data=out,dates=dates,simVar=simVar,modelTag=modelTag[1])
  
  #WRITE TO FILE
  if(IOmode!="suppress"){
    write.table(simDat,file=file,row.names=FALSE,quote = FALSE,sep=",")
  }

  return(simDat)
}
# data(tankDat); obs=tank_obs
# modelTag=c("P-har12-wgen","Temp-har26-wgen")
# pars=modCalibrator(obs=obs,modelTag=modelTag)  #rethink par configuration
# pdf(file="testPlots_yesShufffle_nw20_nshuffle120.pdf",paper="a4")
# par(mfrow=c(5,1),oma=rep(0,4),mar=c(2,2,1,1))
# time=rep(0,100)
# for(i in 1:100){
#   # if(i==6){windows(); par(mfrow=c(5,1),oma=rep(0,4),mar=rep(1,4))}
#   A=Sys.time()
#   sim=modSimulator(datStart="1950-01-01",
#                    datFinish="1999-12-31",
#                    modelTag=modelTag,
#                    parS=pars,
#                    seed=i,
#                    file=paste0("tester",i,".csv"),
#                    IOmode="verbose")
#   
#   plot(sim$P[1:(365*4)])
#   B=Sys.time()
#   print(B-A)
#   time[i]=as.numeric(B-A)
# }
# mean(time)
# dev.off()

#----------------------------------------------------------------------------------
#Argument checker for modSimulator
argument_check_simulator<-function(modelTag=NULL,
                                   datStart=NULL,
                                   datFinish=NULL,
                                   parS=NULL){
  
  #CHECKS FOR MODELTAGS
  # if (modelTag[1]=="Simple-ann") { stop("Simple scaling does not require calibration - invalid request")}
  if (anyDuplicated(modelTag)!=0) {stop("There are multiple entries of the same model tag")}
  for(i in 1:length(modelTag)){
    if(sum(modelTag[i] %in% modelTaglist)==0){
      stop(paste0("modelTag ",i," unrecognised"))
    }
  }
  
  #Check parameters are specified
  if(is.null(parS)){
    stop("No model parameters specified via parS argument. Parameters are required for each modelTag.")
  }
  
  #Check that each model has parameters specified in a list
  modelPar=ls(parS)             #names of models in par list
  for(i in 1:length(modelTag)){
    if(sum(modelTag[i] %in% modelPar)==0){
      stop(paste0("modelTag ",i,", parameters missing. No parameters supplied in 'pars' list"))
    }
  }
  
  #CHECK length of parameter vectors
  modelInfo=get.multi.model.info(modelTag=modelTag)
  for(i in 1:length(modelTag)){
    #check parS length is equaivalent
    if(length(parS[[modelTag[i]]])!=modelInfo[[modelTag[i]]]$npars){
      stop(paste0("Parameters missing for modelTag: ",modelTag[i],".\nMismatch in length of supplied parameter vector pars$",modelTag[i],".\nModel requires:", paste(modelInfo[[modelTag[i]]]$parNam,collapse=" ")))
    }
  }
  
  #Check dates
  if(datStart>datFinish){
    stop("Check supplied dates. datFinish must occur after datStart. Must use recognised date format: '1990-10-01', '01/10/1990', etc ")
  }
  
}



