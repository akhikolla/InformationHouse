#######################################################
##    FUNCTIONS TO SAMPLE EXPOSURE SPACE TARGETS     ##
#######################################################


#THIS ALL NEEDS REVISING####

# CONTAINS
  # genSampBounds
  # expSpaceSampManager() - mange exposure space samples
  # expSpaceSampler() - controller function
  # regSample() - sample on regular grid
  # randUniSample() - sample using uniform random func
  # lhsSample() - sample using latin hypercube sampling (requires external package, lhs)
  # checkProvidedMat() - if matrix of targets is provide check it is suitable
  # checkMatGeneric() - checks matrix

#LIBRARIES
#library(lhs)

#generic bounds for each varType
genSampBounds=list(
  P=c(0.7,1.3),
  Temp=c(-6,6),
  uz=c(-10,10),
  RH=c(-10,10),
  PET=c(0.7,1.3)
)

#FUNCTIONS
expSpaceSampManager<-function(exSpArgs=NULL,
                              attInfo=NULL,
                              attSel=NULL,
                              file=NULL
                              # Anjana Commented: decouple seed from exposure space
                              # IOmode=NULL,
                              # nSeed=NULL,
                              # seedID=NULL,
                              # arrayID=NULL
                              ){
  
  if(is.character(exSpArgs)==TRUE) {
    print("READING ATTRIBUTE TARGETS FROM FILE")
    targetMat=read.table(file=exSpArgs,sep=",",header=TRUE)  #IF TARGET MAT IN CSV
    attRot=NULL
  }else{
    if(exSpArgs$type!="OAT"){
      targetMat=expSpaceSampler(type=exSpArgs$type,
                                samp=exSpArgs$samp,
                                bounds=exSpArgs$bounds,
                                varType=attInfo$varType,
                                targetType=attInfo$targetType,
                                attSel=attSel,
                                file)
      attRot=NULL
    }else{
      spaceInfo=expSpaceSampler(type=exSpArgs$type,
                                samp=exSpArgs$samp,
                                bounds=exSpArgs$bounds,
                                varType=attInfo$varType,
                                targetType=attInfo$targetType,
                                attSel=attSel,
                                file)
      targetMat=data.frame(spaceInfo$targets)
      attRot=spaceInfo$attRot
    }
  }
  
  # Anjana Commented: decouple seed from exposure space
  # nTarg=dim(targetMat)[1]
  # 
  # #IF IN DEV MODE AND NSEEDS USED
  # tmpMat=NULL; tmpRot=NULL; seedA=NULL
  # if((seedID == "fixed") & (!is.null(nSeed))){
  #   #CREATE SEED LIST
  #   for(i in 1:nSeed){
  #     tmpSeed=rep((1234+i),nTarg)
  #     tmpMat=rbind(tmpMat,targetMat)
  #     tmpRot=rbind(tmpRot,attRot)
  #     seedA=c(seedA,tmpSeed)
  #   }
  #   targetMat=tmpMat
  #   attRot=tmpRot
  #   seedCatalogue=seedA
  # }else if (seedID=="arrayID"){
  #   seedCatalogue=rep(arrayID,nTarg)
  # } else {
  #   seedCatalogue=rep(seedID,nTarg)
  # }
  # 
  # return(list(targetMat=targetMat,attRot=attRot,seedCatalogue=seedCatalogue))
  
  return(list(targetMat = targetMat,
              attRot = attRot))
  
}

expSpaceSampler<-function(type="regGrid",  #LHS to come
                          samp=10,   #scalar or vector (in order of attributes)
                          bounds=genSampBounds,
                          varType=NULL,
                          targetType=NULL,    #"frac","pc","diff"
                          attSel=NULL,file
){
  
  switch(type,
         "regGrid" = {targets=regSample(samp=samp, 
                                        bounds=bounds,
                                        varType=varType,
                                        targetType=targetType,
                                        attSel=attSel)},
         "OAT" = {targets=oatSample(samp=samp, 
                                    bounds=bounds,
                                    varType=varType,
                                    targetType=targetType,
                                    attSel=attSel)},
         "LHS" = {err="LHS SAMPLING NOT YET IMPLEMENTED"
                  logfile(p("ERROR: ",err),file)
                  stop(err)},
         {err="INVALID SAMPLING OPTION REQUESTED"
         logfile(p("ERROR: ",err),file)
         stop(err)}
         )
  return(targets)
}


oatSample<-function(samp=30,   #scalar or vector (in order of attributes)
                    bounds=genSampBounds,
                    varType=NULL,
                    targetType=NULL,    #"frac","pc","diff"
                    attSel=NULL){

  temp=samp
  temp[which(temp==1)]=0
  n<-sum(temp)
  targets=matrix(NA,nrow=n,ncol=length(attSel))
  for(i in 1:length(attSel)){
    if(varType[i]=="P"){
      targets[,i]=1
    } else if (varType[i]=="Temp") {
      targets[,i]=0
    } else {
      targets[,i]=1
    }
  }
  
  attRot=rep(NA,n) #create vector to store attRot (for OAT sims)
  
  k=1
  for(i in 1:length(attSel)){
    if(samp[i]!=1){
      step=(bounds[[i]][2]-bounds[[i]][1])/(samp[i]-1)
      for(j in 1:samp[i]){
        targets[k,i]=bounds[[i]][1]+(j-1)*step 
        attRot[k]=attSel[i]
        k=k+1
      }
    } #else {
      # if (length(bounds[[i]]) == 1) {
      #   targets[k,i]=bounds[[i]][1]
      # } else {
      #   targets[k,i]=bounds[[i]][1]+((bounds[[i]][2]-bounds[[i]][1])/2)
      # }
      # attRot[k]=attSel[i]
    # }
  }
  colnames(targets)=attSel
  return(list(targets=targets,attRot=attRot))
}


#TESTER
# attSel_test=c("P_ann_tot_m","Temp_ann_dyAll_m","P_ann_tot_med")
# varType_test=c("P","Temp","P")
# targetType_test=c("frac","diff","frac")
# test=expSpaceSampler(type="regGrid",samp=10,bounds=genSampBounds,varType=varType_test,targetType=targetType_test,attSel=attSel_test)

regSample<-function(samp=30,   #scalar or vector (in order of attributes)
                    bounds=genSampBounds,
                    varType=NULL,
                    targetType=NULL,    #"frac","pc","diff"
                    attSel=NULL
                    
){
  #SAMPLE HOW MANY TIMES IN EACH DIM
  if((length(samp)!=length(varType)) & (length(samp)!=1)){
    samp=samp[1]
    nSamp=1
  }else{
    nSamp=length(samp)
  }
  
  #DETERMINE IF USEFUL BOUNDS SUPPLIED PER VARTYPE OR ATTSEL (IF NOT TERMINATE)
  checkedBounds=checkBounds(bounds=bounds,varType=varType,attSel=attSel)

  #Atts bounds specified
  dimSet=list()
  dimNam=ls(checkedBounds$bound)
  if(length(dimNam)==length(attSel)){
    for(i in 1:length(checkedBounds$bound)){
      if(length(checkedBounds$bound[[i]])==1){
        dimSet[[attSel[i]]]=c(checkedBounds$bound[[i]])
      }else{
        if(nSamp>1){int=(samp[i]-1)}else{int=(samp-1)}
        rng=checkedBounds$bound[[i]][2]-checkedBounds$bound[[i]][1]
        BY=rng/int
        dimSet[[attSel[i]]]=seq(checkedBounds$bound[[i]][1],checkedBounds$bound[[i]][2],by=BY) #store in list
      }
    }
  }else{
    err=p("Incorrect sampling bounds supplied. Program terminated")
    logfile(p("ERROR: ",err),file)
    stop(err)
  }
  
  #MAKE SETS INTO TARGETS
  targets=expand.grid(dimSet,KEEP.OUT.ATTRS = FALSE)  #keeps names $P,$Temp $uz etc
  
  #CHECK FOR ZERO ROW (MAKE SURE BASELINE CLIMATE IS SIMULATED)
  #targets=checkProvidedMat(mat=targets,targetType=targetType)
  
  return(targets)
  
}

#TESTER
# attBounds=list(
#   P_ann_tot_m=c(0.7,1.3),
#   Temp_ann_dyAll_m=c(0),
#    P_ann_tot_med=c(0.5,1.3)
# )
# attSel_test=c("P_ann_tot_m","Temp_ann_dyAll_m")
# varType_test=c("P","Temp","P")
# targetType_test=c("frac","diff","frac")
# regSample(samp=5,bounds=attBounds,varType=varType_test,targetType=targetType_test,attSel=attSel_test)

# attSel_test=c("P_ann_tot_m","Temp_ann_dyAll_m","P_ann_tot_med")
# varType_test=c("P","Temp","P")
# targetType_test=c("frac","diff","frac")
# regSample(samp=5,bounds=genSampBounds,varType=varType_test,targetType=targetType_test,attSel=attSel_test)

#------------------------------
checkBounds<-function(bounds=NULL,
                      varType=NULL,
                      attSel=NULL
){
  
  vars=varType[which(!duplicated(varType))] #remove duplicates from vector

  tmpAtt=intersect(ls(bounds),attSel)      # determine bounds
  tmpVar=intersect(ls(bounds),varType) 
  
  nAttBound=length(tmpAtt)
  nVarBound=length(tmpVar)
  
  #IF NEITHER FIT
  if((nAttBound==0) & (nVarBound==0)){
    err=p("Incorrect sampling bounds supplied. Program terminated")
    logfile(p("ERROR: ",err),file)
    stop(err)
  } 
  
  #LIKELY ATT BOUNDED
  if(nAttBound>nVarBound){
      if(nAttBound==length(attSel)){  #check if enough bounds specified
        boundType="Att"
        out.bound=list()
        for(i in 1:length(attSel)){out.bound[[attSel[i]]]=bounds[[attSel[i]]]} #copy needed bounds
      }else{
        err=p("Incorrect sampling bounds supplied. Program terminated")
        logfile(p("ERROR: ",err),file)
        stop(err)
      }
  }else{  #LIKELY VAR BOUNDED
    if(nVarBound==length(vars)){ 
        boundType="Var"
        out.bound=list()
        for(i in 1:length(varType)){
          out.bound[[attSel[i]]]=bounds[[varType[i]]] #copy needed bounds
          } 
    }else{
      err=p("Incorrect sampling bounds supplied. Program terminated")
      logfile(p("ERROR: ",err),file)
      stop(err)
    }
  }
  
  out=list(boundType=boundType,
           bound=out.bound,
           tmpAtt=tmpAtt
           )
  
  return(out)
}


#checkBounds(bounds=genSampBounds,varType=varType,attSel=attSel)

#ALL REQUIRE NO. OF ATTRIBUTES, BOUNDS
checkProvidedMat<-function(mat=NULL,                   # supplied matrix - each row is a target location
                           targetType=NULL    # type (will be vector of "pc","frac" or "diff')
                           ){

  #row to join to 
  nr=dim(mat)[2]  #nrows
  zero.row=rep(0,nr)
  
  #Make baseline climate row
  ind.pc=which(targetType=="pc")
  ind.frac=which(targetType=="frac")
  ind.diff=which(targetType=="diff")
  zero.row[ind.pc]=100
  zero.row[ind.frac]=1
  zero.row[ind.diff]=0
  
  mat=checkMatGeneric(mat=mat,zero.row=zero.row)

  return(mat)
}

checkMatGeneric<-function(mat=NULL,       # matrix - each row is a target location
                          zero.row=NULL   # row that signifies no change in any attribute e.g. 0, 1, 100
                          ){
  base.row=NULL
  
  for(i in 1:nrow(mat)){
    temp=(mat[i,] == zero.row)
    temp=which(temp == FALSE)
    if(length(temp)==0){
      base.row=c(base.row,i)
    }
  }
  
  if(length(base.row)==0){
    mat=rbind(zero.row,mat,deparse.level=0)                 #if no row of vals add it in
  }else{
    mat=rbind(zero.row,mat[-base.row,],deparse.level = 0) # make base row the 1st entry and remove any duplicate
  }
  return(mat)
}

#TESTER
# mat=matrix(c(1,2,4,1,2,3,4,4,4),nrow=3,ncol=3)
# checkMatGeneric(mat=mat,zero.row=c(1,0,0))
# checkProvidedMat(mat=mat,targetType=targetType)
