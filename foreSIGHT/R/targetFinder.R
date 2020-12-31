########################################################
## FUNCTION TO SIMULATE SERIES AND DETERMINE FITNESS ###
########################################################

#CONTAINS
  #targetFinder() -
#------------------------------------------------------------------------------------------------ 
targetFinder<- function(parS,               # vector of pars (will change in optim)
                        modelInfo=NULL,
                        modelEnv=NULL,      # tag to link/identify model
                        attSel=NULL,        # attributes selected (vector of strings)
                        attPrim=NULL,       # primary attribute label
                        attInfo=NULL,     # added info regarding attributes (maybe add in attPrim here)!!!!!!!!!!!!!
                        datInd=NULL,
                        randomVector = NULL,
                        target=NULL,        # target locations: desired changes in climate to be simulated, in % relative or abs diff to baseline levels (vector)
                        attObs=NULL,        # observed series attribute values
                        lambda.mult=NULL,   # lambda multiplier for penalty function
                        simSeed=NULL,
                        wdSeries=NULL,      
                        resid_ts=NULL
                        #Nw=NULL,            # warmup period in days
                        # N=NULL,             # seeds
                        # seed1=NULL,
                        # seed2=NULL,
                        # seed3=NULL
){
  
  #SIMULATE SELECTED VARIABLE USING CHOSEN STOCHASTIC MODEL
 
  sim=switch_simulator(type=modelInfo$simVar,          # what vartype is being simulated
                       parS=parS,
                       modelEnv=modelEnv,      
                       randomVector = randomVector,
                       wdSeries=wdSeries,      
                       resid_ts=resid_ts,
                       seed=simSeed)
  
  if(length(which(is.na(sim$sim))) > 0){
    score=-150  #default here
  }else{
    #CALCULATE SELECTED ATTRIBUTE VALUES
    sim.att=attribute.calculator(attSel=attSel,data=sim$sim,datInd=datInd,attribute.funcs=attribute.funcs)
   
    #RELATING TO BASELINE SERIES 
    simPt=unlist(Map(function(type, val,baseVal) simPt.converter.func(type,val,baseVal), attInfo$targetType, as.vector(sim.att),as.vector(attObs)),use.names = FALSE)
    
    #GET OBJECTIVE FUNCTION VALUE ()
    score=objFuncMC(attSel= attSel,     # vector of selected attributes 
                    attPrim=attPrim,      # any primary attributes
                    attInfo=attInfo,
                    simPt=simPt,
                    target=target,
                    penalty.func=penaltyFunc_basic,   #make this changeable (auto calc lambda)
                    lambda=lambda.mult
    )
  }
 
  
  #cat(score)
  return(score)
}









