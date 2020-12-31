###############################
##      OUTPUT MANAGER      ##
###############################

#CONTAINS
  # saveTarget() 
    # nameMaker()  
    # makeOutputDataframe()
  # csvMaker()
    # writeToCSV()


#------------------------------------------
#FUNCTIONS

nameMaker<-function(attSel=NULL,  # vector of selected attributes
                    target=NULL   # vector coordinates of targets
){
  #SQUASH ATTSEL CHARACTERS
  squashfunc<-function(x) {paste(x,sep="",collapse="")}
  nAtt=length(attSel)
  tmp.str=strsplit(x=attSel,split="_");tmp.str=vapply(tmp.str,FUN=squashfunc,FUN.VALUE=character(1))
  #ATTACH TARGET COORDINATES
  tarVal=vapply(X=as.numeric(target),FUN=signif,digits=5,FUN.VALUE=numeric(1))
  tarVal=substr(tarVal,start=1,stop=5)
  fnam=paste(paste(tarVal[-nAtt],tmp.str[-nAtt],"_",sep="",collapse=""),paste(tarVal[nAtt],tmp.str[nAtt],sep="",collapse=""),sep="")
  #RETURN NAME MINUS FILE EXTENSION
  return(fnam)
}

makeOutputDataframe<-function(data=NULL,
                              dates=NULL,
                              simVar=NULL,
                              modelTag=NULL
){
  #FUNCTION TO RETURN A SUBSECTION OF THE LIST
  switch(modelTag,
         "Simple-ann" = {returnTS=function(data,simVar){data[[simVar]]}},
                        {returnTS=function(data,simVar){data[[simVar]]$sim}}
         )

  #MAKE OUTPUTTED DATAFRAME
  outDat=vapply(simVar,FUN=returnTS,data=data,FUN.VALUE=numeric(length(dates[,1])))
  outDat=cbind(dates,outDat)
  names(outDat)=c(names(dates),simVar)
  
  return(outDat)
  
}

saveTarget<-function(data=NULL,       # data[[i]]  ->    $P, $Temp $attSim $targetSim
                     dates=NULL,      # data frame dates info mm,dd,yy
                     modelTag=NULL,   # vector of modelTags
                     modelInfo=NULL,  # list of modelInfo related to modelTag
                     simVar=NULL,     # vector of variables simulated
                     target=NULL,     # vector coords of target location
                     attSel=NULL,     # vector of attributes Selected
                     attPrim=NULL,     # vector of primary attributes nominated
                     paths=NULL
){
  #MAKE FILENAME
  fnam=nameMaker(attSel=attSel,target=target)
  
  #MAKE OUTPUTTED DATAFRAME
  simDat=makeOutputDataframe(data=data,dates=dates,simVar=simVar,modelTag=modelTag[1])
  
  #RENAME WHAT YOU WANT TO SAVE AT TOP LEVEL
  attSim=data$attSim
  par=data$parS
  targetSimulated=data$targetSim
  targetRequested=target
  seed=data[[1]]$seed   #this assume same seed for all
  
  #WRITE TO .RDATA FILE
  save(simDat,modelTag,modelInfo,attSel,attPrim,attSim,targetSimulated,targetRequested,seed,par,file=paste(paths$RData,"/",fnam,".RData",sep=""))
  
  #WRITE OUTPUT CSV
  write.table(simDat,file=paste(paths$CSV,"/",fnam,".csv",sep=""),row.names=FALSE,quote = FALSE,sep=",")  
    
  return(fnam)
}

# #simpleSaveTarget
# simpleSaveTarget<-function(data=NULL,       # data[[i]]$P, $Temp $attSim $targetSim
#                            dates=NULL,      # data frame dates info mm,dd,yy
#                            simVar=NULL,     # vector of variables simulated
#                            attSel=NULL,     # vector of selected attributes
#                            target=NULL,      # vector coordinates of targets
#                            modelTag=NULL,
#                            modelInfo=NULL,
#                            paths=NULL,
#                            suffixFileName = NULL
# ){
# 
#   # MAKE FILENAME (SANS .CSV EXTENSION)
#   # fnam=nameMaker(attSel=attSel,target=target)
#   fnam = paste0("sim_", suffixFilename)
# 
#   #MAKE OUTPUTTED DATAFRAME
#   simDat=makeOutputDataframe(data=data,dates=dates,simVar=simVar,modelTag=modelTag)
#   
#   #WRITE OUTPUT CSV
#   write.table(simDat,file=paste(paths$CSV,"/",fnam,".csv",sep=""),row.names=FALSE,quote = FALSE,sep=",")
#   
#   #WRITE TO RDATA
#   save(simDat,modelTag,modelInfo,attSel,target,file=paste(paths$RData,"/",fnam,".RData",sep=""))
#   return(fnam)
# }

#WRITING TO CSV
writeToCSV<-function(data=NULL,
                     dates=NULL,
                     fnam=NULL,
                     simVar=NULL,
                     modelTag=NULL
){
  
  #ADD CSV TO FILENAME
  filename=paste(fnam,".csv",sep="")
  
  #MAKE OUTPUTTED DATAFRAME
  outDat=makeOutputDataframe(data=data,dates=dates,simVar=simVar,modelTag=modelTag)
  
  #WRITE OUTPUT CSV
  write.table(outDat,file=filename,row.names=FALSE,quote = FALSE,sep=",")
  
}

# saveTargets<-function(data=NULL,       #data[[i]]$P, $Temp $attSim $targetSim
#                       dates=NULL,      #data frame dates info mm,dd,yy
#                       modelTag=NULL,   #vector of modelTags
#                       target=NULL,     #matrix of targets
#                       attSel=NULL,     #vector of attributes Selected
#                       attPrim=NULL     #vector of primary attributes nominated
#                           ){
#   
#   nTarget=nrow(target)
#   flist=rep(NA,nTarget)
#   
#   for(i in 1:nTarget){
#     
#     sim=data[[i]];sim$targetRequested=target[i,]
#     #make filename
#     # tarVal=vapply(X=as.numeric(sim$targetRequested),FUN=signif,digits=2,FUN.VALUE=numeric(1)); tarVal=substr(tarVal,start=1,stop=4)
#     # fnam=paste(paste(tarVal[-nAtt],tmp.str[-nAtt],"_",sep="",collapse=""),paste(tarVal[nAtt],tmp.str[nAtt],sep="",collapse=""),sep="")
#     
#     fnam=nameMaker(attSel=attSel,target=target[i,])
#     flist[i]=paste(fnam,".RData",sep="")
#     
# 
#     #Save what you want
#     save(sim,dates,modelTag,attSel,attPrim,file=flist[i])
#   }
# 
#   #OUTPUT FILENAME LIST
#   write.table(flist,file = "filenames.txt",row.names=FALSE,quote = FALSE,col.names=c("Filename"))
# }


#TEST
#saveTargets(data=data,dates=obs[,c("year","month","day")],modelTag=modelTag,target=targetMat[1:1,], attSel=attSel,attPrim=attPrim)
# obs[,c("year","month","day")]
# targetMat[]




  

# csvTargets<-function(data=NULL,       #data[[i]]$P, $Temp $attSim $targetSim
#                      dates=NULL,      #data frame dates info mm,dd,yy
#                      modelTag=NULL,   #vector of modelTags
#                      target=NULL,     #matrix of targets
#                      attSel=NULL,     #vector of attributes Selected
#                      simVar=NULL
# ){
#   
#   
#   squashfunc<-function(x) {paste(x,sep="",collapse="")}
#   nAtt=length(attSel)
#   tmp.str=strsplit(x=attSel,split="_");tmp.str=vapply(tmp.str,FUN=squashfunc,FUN.VALUE=character(1))
#   nTarget=nrow(target)
#   
#   for(i in 1:nTarget){
#     #MAKE FILENAME FOR TARGET
#     tarVal=vapply(X=as.numeric(data[[i]]$targetSim),FUN=signif,digits=2,FUN.VALUE=numeric(1)); tarVal=substr(tarVal,start=1,stop=4)
#     fnam=paste(paste(tarVal[-nAtt],tmp.str[-nAtt],"_",sep="",collapse=""),paste(tarVal[nAtt],tmp.str[nAtt],sep="",collapse=""),sep="")
#     
#     #WRITE TO CSV
#     writeToCSV(data=data[[i]],dates=dates,fnam=fnam,simVar=simVar)
#   }
# 
# }



# dates=tmp[,c("year","month","day")]
# writeToCSV(data=sim[[i]],dates=dates,fnam=fnam,simVar=simVar)
# 
# saveTarget<-function(data=NULL,       #data[[i]]$P, $Temp $attSim $targetSim
#                       dates=NULL,      #data frame dates info mm,dd,yy
#                       modelTag=NULL,   #vector of modelTags
#                       target=NULL,     #matrix of targets
#                       attSel=NULL,     #vector of attributes Selected
#                       attPrim=NULL     #vector of primary attributes nominated
# ){
#   
#   
#   squashfunc<-function(x) {paste(x,sep="",collapse="")}
#   nAtt=length(attSel)
#   tmp.str=strsplit(x=attSel,split="_");tmp.str=vapply(tmp.str,FUN=squashfunc,FUN.VALUE=character(1))
#   
#   nTarget=nrow(target)
#   flist=rep(NA,nTarget)
#   
#   for(i in 1:nTarget){
#     
#     sim=data[[i]];sim$targetRequested=target[i,]
#     #make filename
#     tarVal=vapply(X=as.numeric(sim$targetRequested),FUN=signif,digits=2,FUN.VALUE=numeric(1)); tarVal=substr(tarVal,start=1,stop=4)
#     fnam=paste(paste(tarVal[-nAtt],tmp.str[-nAtt],"_",sep="",collapse=""),paste(tarVal[nAtt],tmp.str[nAtt],sep="",collapse=""),sep="")
#     flist[i]=paste(fnam,".RData",sep="")
#     
#     
#     #Save what you want
#     save(sim,dates,modelTag,attSel,attPrim,file=flist[i])
#   }
#   
#   #OUTPUT FILENAME LIST
#   write.table(flist,file = "filenames.txt",row.names=FALSE,quote = FALSE,col.names=c("Filename"))
# }


