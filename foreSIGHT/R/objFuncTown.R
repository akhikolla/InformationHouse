#######################################
##    OBJECTIVE FUNCTION LIBRARY     ##
#######################################

# CONTAINS FUNCTIONS
  # objFuncMC() - objective function controller
  # objFuncSelector() - selects version of objective function to use (?)
  # eucDist() - euclidian distance function
  # penaltyFunc_basic() - applies lambda multiplier to distance (target- attrib)
  # penaltyFunc_basicSq() - applies lambda multiplier to distance (target- attrib)^2
    #to come - # penaltyFuncSelector() - selects version of penalty function to use
#-------------------------------------------------------------------
objFuncMC<-function(attSel= NULL,     # vector of selected attributes 
                    attPrim=NULL,     # any primary attributes
                    attInfo=NULL,     # added info regarding attributes (maybe add in attPrim here)
                    simPt=NULL,
                    target=NULL,
                    penalty.func=penaltyFunc_basic,   #also penaltyFunc_basicSq
                    lambda=0          # Anjana: Note that this lambda is unused after the fix for lambda.mult exclusive to attPrim - to be removed
){
  
  # get.ind<-function(x,y){which(x == y)}     # quick function to find which are primary attributes
  
  nAtt=length(attSel)  #how many attributes
  
  #switch from list to vector for eucdist  - should just be able to unlist if only those attribs are calculated
  dist=eucDist(target=as.double(target),simPt=simPt)
  
  primInd=which(attInfo$primType==TRUE)
  if(length(primInd)>0){
    penalty.score=penalty.func(target=as.double(target)[primInd],simPt=simPt[primInd],lambda=attInfo$primMult[primInd]) 
    score=-dist-penalty.score
  }else{
    score=-dist
  }

  
  
  # #FIND RELEVANT ATTPRIMS
  # attPrimSub=intersect(attPrim,attSel) 
  # 
  # # PENALTY FUNCTION
  # if(length(attPrimSub)>0){
  #   #IDENTIFY WHERE PRIMARY ATTRIBUTES ARE LOCATED IN ATTSEL
  #   primInd=vapply(attPrimSub,FUN=get.ind,FUN.VALUE=numeric(1),x=attSel,USE.NAMES = FALSE)  #Indices of primary attributes
  #   
  #   #FIND OUT WHICH ATTPRIM LIVES HERE (AS LAMBDA.MULT STORED IN ORDER)
  #   lamInd=vapply(attPrimSub,FUN=get.ind,FUN.VALUE=numeric(1),x=attPrim,USE.NAMES = FALSE)  #Indices of primary attributes
  # 
  #   #SUPPLY A PENALTY FUNCTION
  #   penalty.score=penalty.func(target=target[primInd],simPt=simPt[primInd],lambda=lambda[lamInd])  #fix this up for different multi-lambda's
  #   score=-dist-penalty.score
  # }else{
  #   score=-dist
  # }
  
  #SCORE IS MADE -'VE AS GA WORKS ON MAXIMISATION
  return(score) 
  
}  
#----------------------------------------------------------
eucDist<-function(target=NULL,  #vector of target locations
                  simPt=NULL    #vector of simulated locations
){
  score <- sqrt(sum((target - simPt)^2L)) 
  
}

#----------------------------------------------------------
penaltyFunc_basic<-function(target=NULL,  #scalar target
                            simPt=NULL,   #scalar sim point
                            lambda=NULL   #multiplier/tuning parameter
                            ){
  penalty=sum(lambda*abs(target-simPt),na.rm=TRUE)
  #return(penalty)
}
#----------------------------------------------------------
penaltyFunc_basicSq<-function(target=NULL,  #scalar target
                            simPt=NULL,   #scalar sim point
                            lambda=NULL   #multiplier/tuning parameter
){
  penalty=sum(lambda*abs(target-simPt)^2L,na.rm=TRUE)
}
#----------------------------------------------------------
simPt.converter.func<-function(type=NULL,      # type of simPt
                          val=NULL,       # value of simulated attribute in normal space
                          baseVal=NULL    # value of observed series attribute in normal space
){
#simPt=NULL
switch(type,
       "frac" = {simPt=(val/baseVal)},
       "pc" = {simPt=(val-baseVal)/baseVal*100},
       "diff" = {simPt=(val-baseVal)},
       -999
       )
  
  return(simPt)
 
}
# do a vector form...(probably not needed)
# simPt.converter.func(type=typ, val=vals, baseVal=baseVals)


