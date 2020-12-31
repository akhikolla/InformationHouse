#######################################
##   OPTIMIZER MANAGER FUNCTIONS     ##
#######################################

# CONTAINS
  # gaWrapper() - allows use of suggestions for initial populations
  # screenSuggest() - remove suggestions outside bounds
   # enforceBounds()
      # testBound()
        # outBound()

#-----------------------------------------------
gaWrapper<-function(gaArgs=NULL,        # can specify your own outside
                    modelEnv=NULL,
                    modelInfo=NULL,     # information related to modelTags
                    attSel=NULL,        # attributes selected (vector of strings)
                    attPrim=NULL,       # primary attribute label
                    attInfo=NULL,       # added info regarding attributes
                    datInd=NULL,
                    randomVector = NULL,
                    parSuggest=NULL,    # paramater suggestions
                    target=NULL,        # target locations: desired changes in climate to be simulated, in % relative or abs diff to baseline levels (vector)
                    attObs=NULL,        # observed series attribute values
                    lambda.mult=NULL,   # lambda multiplier for penalty function
                    simSeed=NULL,       # seeds
                    wdSeries=NULL,      
                    resid_ts=NULL
                    ){
  
  #MEMOISE YAY OR NAY - TBC
  #USER SPECIFIED PARALLEL CONTROLS
  timeStart=Sys.time()
  optpar<- ga(type = "real-valued",
              fitness=targetFinder,
              lower = modelInfo$minBound,
              upper = modelInfo$maxBound,
              pcrossover= gaArgs$pcrossover,
              pmutation=gaArgs$pmutation,
              maxiter=gaArgs$maxiter, 
              popSize = gaArgs$popSize, 
              maxFitness = gaArgs$maxFitness, 
              run=gaArgs$run, 
              seed = gaArgs$seed, 
              parallel = gaArgs$parallel,
              keepBest=gaArgs$keepBest, 
              suggestions = parSuggest, 
              monitor = FALSE,             #switchback
              #BELOW RELATED TO TARGETFINDER()
              modelInfo=modelInfo,
              modelEnv=modelEnv,
              attSel=attSel,        
              attPrim=attPrim,       
              attInfo=attInfo, 
              datInd=datInd,
              randomVector = randomVector,
              target=target,        
              attObs=attObs,        
              lambda.mult=lambda.mult,   
              simSeed=simSeed,
              wdSeries=wdSeries,      
              resid_ts=resid_ts
              )
  timeFin=Sys.time()
  timeRun=timeFin-timeStart    #optimisation runtime
  #print(summary(optpar)$fitness)
  
  out=list(par=as.vector(optpar@solution[1,]),
           fitness=as.numeric(optpar@fitnessValue), 
           seed=simSeed,
           opt=optpar,
           runtime=timeRun)
           
  return(out)
  
}

#FUNCTION TO SCREEN SUGGESTED POPULATIONS (screen outside, screen once)
screenSuggest<-function(modelInfo=NULL,
                        modelTag=NULL,
                        parLoc=NULL,        # which pars belong to which model parLoc[[mod]]=c(start,end)
                        suggest=NULL        # suggested pars
){
  nMod=length(modelTag)
  ind=NULL
  for(mod in 1:nMod){
    parSel=suggest[,(parLoc[[mod]][1]:parLoc[[mod]][2])]              #grab par suggestions related to modelTag running
    tmpInd=enforceBounds(suggest=parSel,                              #matrix of suggestions
                         minBound=modelInfo[[modelTag[mod]]]$minBound,
                         maxBound=modelInfo[[modelTag[mod]]]$maxBound)
    #JOIN THE INDICES TOGETHER
    ind=c(ind,tmpInd)
  } 
  
  #REMOVE INDICES DUPLICATES
  ind=unique(ind)
  
  #REMOVE INAPPROPRIATE SUGGESTIONS
  suggest=suggest[ind,]
  return(suggest)
}

enforceBounds<-function(suggest=NULL, #matrix of suggestions
                        minBound=NULL,
                        maxBound=NULL
){
  nSuggest=nrow(suggest)
  ind=vapply(X=seq(1,nSuggest),FUN=testBound,matPar=suggest,minBound=minBound,maxBound=maxBound,FUN.VALUE=numeric(1))
  ind=ind[which(!is.na(ind))]
  ind
}

testBound<-function(ind=NULL,
                    matPar=NULL,
                    minBound=NULL,
                    maxBound=NULL
){
  
  out=outBound(ind=ind,inVector=matPar[ind,],minBound=minBound,maxBound=maxBound)
  out
}
# matTest=c(1,2,3,4,4,4,5,5,5,10,11,12);matTest=matrix(matTest,ncol=3,nrow=4,byrow=TRUE)
# minBound=c(1,1,3); maxBound=c(6,8,13)

outBound<-function(ind=NULL,      #index of pars being evaluated
                   inVector=NULL,
                   minBound=NULL,
                   maxBound=NULL
){
  npar=length(minBound)  #No. of pars that should be inside bound
  ntest=length(which((inVector>minBound)&(inVector<maxBound)))
  if(ntest==npar){ok=ind}else{ok=NA}
  return(ok)
}

# #
# optpar<- ga(type = "real-valued",
#             fitness=targetFinder,
#             min = modelInfo[[2]]$minBound,
#             max = modelInfo[[2]]$maxBound,
#             pcrossover= gaArgs$pcrossover,
#             pmutation=gaArgs$pmutation,
#             maxiter=gaArgs$maxiter, 
#             popSize = gaArgs$popSize, 
#             maxFitness = gaArgs$maxFitness, 
#             run=gaArgs$run, 
#             seed = gaArgs$seed, 
#             parallel = gaArgs$parallel,
#             keepBest=gaArgs$keepBest, 
#             monitor = FALSE, 
#             #BELOW RELATED TO TARGETFINDER()
#             modelTag=modelTag[2],
#             modelInfo=modelInfo[[2]],
#             attSel=attSel[attInd[[2]]],        
#             attPrim=attPrim,       
#             attInfo=attInfo[[2]],     
#             datInd=datInd[[2]],        
#             initCalibPars=initCalibPars, 
#             target=target[attInd[[2]]],        
#             attObs=attObs[attInd[[2]]],        
#             lambda.mult=gaArgs$lambda.mult,   
#             simSeed=simSeed,
#             wdSeries=wdStatus,      
#             resid_ts=resid_ts
# )



#optim=FALSE


# seed1=122
# #.ga.optim.default<-list()
# gaArgs=list(
#   pcrossover= 0.1,
#   pmutation=0.8,
#   maxiter=100,
#   popSize = 100,
#   maxFitness = -0.01,
#   run=300,
#   seed = 1234,
#   parallel = TRUE,
#   keepBest=TRUE
# )

# #TESTER CALL
# optTest=list()
# optTest[[1]]=gaWrapper( gaArgs=gaArgs,  #can specify your own outside
#                       modelTag=modelTag,      # tag to link/identify model
#                       modelInfo=modelInfo,
#                       attSel=attSel,        # attributes selected (vector of strings)
#                       attPrim=attPrim,       # primary attribute label
#                       attInfo=attInfo,       # added info regarding attributes (maybe add in attPrim here)!!!!!!!!!!!!!
#                       datInd=datInd,        # dat ind
#                       initCalibPars=initCalibPars, #vector of pars from initial baseline calibration
#                       target=c(1.05,1.06,1.10),    # target locations: desired changes in climate to be simulated, in % relative or abs diff to baseline levels (vector)
#                       attObs=obs.att,        # observed series attribute values
#                       lambda.mult=2,   # lambda multiplier for penalty function
#                       Nw=100,            # warmup period in days
#                       N=seed1,             # seeds
#                       seed1=seed1,
#                       seed2=seed1,
#                       seed3=seed1
#         )


#Test
# mod=2
# optTest=gaWrapper(gaArgs=.ga.default,          #can specify your own outside
#                     modelTag=modelTag[mod],      # tag to link/identify model
#                     modelInfo=modelInfo[[modelTag[mod]]],
#                     attSel=attSel[attInd[[mod]]],
#                     attPrim=attPrim,
#                     attInfo=attInfo[[modelTag[mod]]],
#                     datInd=datInd[[modelTag[mod]]],
#                     initCalibPars=NULL,
#                     target=targetMat[1,attInd[[2]]],        # target locations: desired changes in climate to be simulated, in % relative or abs diff to baseline levels (vector)
#                     attObs=attObs[attInd[[mod]]],        # observed series attribute values
#                     lambda.mult=1.0,   # lambda multiplier for penalty function
#                     simSeed=1234           # seeds
# )
# 
# plot(optTest$opt)
# 
# out=switch_simulator(type=modelInfo[[modelTag[mod]]]$simVar,parS=optTest$par,
#                     modelTag=modelTag[mod],modelInfo=modelInfo[[modelTag[mod]]],datInd=datInd[[modelTag[mod]]],
#                     initCalibPars=NULL,wdSeries=NULL,resid_ts=NULL,seed=1234)
# 
#   sim.att=attribute.calculator(attSel=attSel[attInd[[mod]]],
#                                    data=out$sim,
#                                    datInd=datInd[[mod]],
#                                    attribute.funcs=attribute.funcs)
