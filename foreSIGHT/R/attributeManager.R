#######################################
##      ATTRIBUTE MANAGER            ##
#######################################

#CONTAINS
  #attribute.funcs - LIST (STORED in...)
  #attribute.calculator() - calculate values of attributes
  #attribute.info.check() - get targetType, varType and identify any invalid model selections
    #check.attribute.model.combo() - check if any attribute-model combos are invalid
    #get.attribute.info() - Identify invalid models
    #get.target.type() - treated as fractions, percent or abs value
    #get.attribute.varType() - "P", "Temp" OR ...


#---------------------------------------------------------------------------------------


threshWD=0.999

#ATTRIBUTE FUNCTION LIST - CAN CALL SPECIFIC ATTRIBUTE CALCS FROM THE FUNCTION LIST
#list of attribute calculation functions - NB: all must have similar/same arguments
#used via lappply (attribute.funcs[attSel], function(f) f(data,datInd)) - return labelled list of outputs
#NEED TO THINK ABOUT WHAT ENVIRONMENT THIS IS STORED IN. #
attribute.funcs=list(
  P_ann_tot_m=function(data,datInd) extractor.summaryMean(func=sum,data=data,indx=datInd$i.yy,nperiod=datInd$nyr),                       # function labelled "Ptot_m" in list #get.tot(data)/datInd$nyr
  
  P_ann_dyWet_m=function(data,datInd) get.wet.average(data,threshold=threshWD),            # function labelled "dyWet_m" in list 

  P_ann_nWet_m=function(data,datInd) extractor.summaryMean(func=get.nwet,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,threshold=threshWD),  
  
  #P_ann_P99_m=function(data,datInd) extractor.summaryMean(func=get.quantile,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,quant=0.99),
  P_ann_P99_m=function(data,datInd) get.quantile(data,quant=0.99),
  
  P_ann_dyWet99p_m=function(data,datInd) extractor.summaryMean(func=get.quantile.wet,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,quant=0.99),
  
  P_ann_avgWSD_m=function(data,datInd) mean(get.spell.lengths(data=data,thresh=0,type="wet"),na.rm=TRUE),
  
  P_ann_avgDSD_m=function(data,datInd) mean(get.spell.lengths(data=data,thresh=0,type="dry"),na.rm=TRUE),
  
  P_JJA_avgWSD_m=function(data,datInd) mean(get.spell.lengths(data=data[datInd$i.ss[[4]]],thresh=0,type="wet"),na.rm=TRUE),
  
  P_MAM_avgWSD_m=function(data,datInd) mean(get.spell.lengths(data=data[datInd$i.ss[[3]]],thresh=0,type="wet"),na.rm=TRUE),
  
  P_DJF_avgWSD_m=function(data,datInd) mean(get.spell.lengths(data=data[datInd$i.ss[[2]]],thresh=0,type="wet"),na.rm=TRUE),
  
  P_SON_avgWSD_m=function(data,datInd) mean(get.spell.lengths(data=data[datInd$i.ss[[1]]],thresh=0,type="wet"),na.rm=TRUE),
  
  P_JJA_avgDSD_m=function(data,datInd) mean(get.spell.lengths(data=data[datInd$i.ss[[4]]],thresh=0,type="dry"),na.rm=TRUE),
  
  P_MAM_avgDSD_m=function(data,datInd) mean(get.spell.lengths(data=data[datInd$i.ss[[3]]],thresh=0,type="dry"),na.rm=TRUE),
  
  P_DJF_avgDSD_m=function(data,datInd) mean(get.spell.lengths(data=data[datInd$i.ss[[2]]],thresh=0,type="dry"),na.rm=TRUE),
  
  P_SON_avgDSD_m=function(data,datInd) mean(get.spell.lengths(data=data[datInd$i.ss[[1]]],thresh=0,type="dry"),na.rm=TRUE),
  
  P_JJA_dyWet_m=function(data,datInd) extractor(func=get.wet.average,data=data,indx=datInd$i.ss[[4]],threshold=0),
  
  P_MAM_dyWet_m=function(data,datInd) extractor(func=get.wet.average,data=data,indx=datInd$i.ss[[3]],threshold=0),
  
  P_DJF_dyWet_m=function(data,datInd) extractor(func=get.wet.average,data=data,indx=datInd$i.ss[[2]],threshold=0),
  
  P_SON_dyWet_m=function(data,datInd) extractor(func=get.wet.average,data=data,indx=datInd$i.ss[[1]],threshold=0),
  
  P_JJA_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.ss[[4]],nblock=datInd$nyr),
  
  P_MAM_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.ss[[3]],nblock=datInd$nyr),
  
  P_DJF_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.ss[[2]],nblock=datInd$nyr),
  
  P_SON_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.ss[[1]],nblock=datInd$nyr),
  
  P_ann_seasRatio_m=function(data,datInd) (extractor(func=get.avg.tot,data=data,indx=c(datInd$i.ss[[3]],datInd$i.ss[[4]]),nblock=datInd$nyr)/extractor(func=get.avg.tot,data=data,indx=c(datInd$i.ss[[1]],datInd$i.ss[[2]]),nblock=datInd$nyr)),

  P_ann_maxWSD_m=function(data,datInd) extractor.summaryMean(func=get.spell.lengths.max,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,thresh=0,type="wet"), #no fortran
  
  P_ann_maxDSD_m=function(data,datInd) extractor.summaryMean(func=get.spell.lengths.max,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,thresh=0,type="dry"), #no fortran
  
  P_ann_R10_m=function(data,datInd) extractor.summaryMean(func=R10calc,data=data,indx=datInd$i.yy,nperiod=datInd$nyr),    #removed fortran .so
  
  Temp_ann_GSL_m=function(data,datInd) extractor.summaryMean(func=GSLcalc,data=data,indx=datInd$i.yy,nperiod=datInd$nyr),    #removed fortran .so
  
  Temp_ann_CSL_m=function(data,datInd) extractor.summaryMean(func=CSLcalc,data=data,indx=datInd$i.yy,nperiod=datInd$nyr),    #remove fortran .so
  
  Temp_ann_avg_m=function(data,datInd) extractor.summaryMean(func=mean,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,na.rm=TRUE),
  
  Temp_ann_P5_m=function(data,datInd) extractor.summaryMean(func=get.quantile,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,quant=0.05),
  
  Temp_ann_P95_m=function(data,datInd) extractor.summaryMean(func=get.quantile,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,quant=0.95),
  
  Temp_ann_F0_m=function(data,datInd) extractor.summaryMean(func=F0calc,data=data,indx=datInd$i.yy,nperiod=datInd$nyr), #add fortran .so
  
  Temp_ann_rng_m=function(data,datInd) extractor.summaryMean(func=get.quantile.rng,data=data,indx=datInd$i.yy,nperiod=datInd$nyr),
  
  PET_ann_avg_m=function(data,datInd) extractor.summaryMean(func=mean,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,na.rm=TRUE),
 
  #PET_ann_rng_m=function(data,datInd) extractor.summaryMean(func=get.quantile.rng,data=data,indx=datInd$i.yy,nperiod=datInd$nyr),
  PET_ann_rng_m=function(data,datInd) get.quantile.rng(data),
  
  PET_ann_P5_m=function(data,datInd) extractor.summaryMean(func=get.quantile,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,quant=0.05),
  
  PET_ann_P95_m=function(data,datInd) extractor.summaryMean(func=get.quantile,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,quant=0.95),
  
  P_Jan_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.mm[[1]],nblock=datInd$nyr),
  
  P_Feb_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.mm[[2]],nblock=datInd$nyr),
  
  P_Mar_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.mm[[3]],nblock=datInd$nyr),
  
  P_Apr_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.mm[[4]],nblock=datInd$nyr),
  
  P_May_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.mm[[5]],nblock=datInd$nyr),
  
  P_Jun_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.mm[[6]],nblock=datInd$nyr),
  
  P_Jul_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.mm[[7]],nblock=datInd$nyr),
  
  P_Aug_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.mm[[8]],nblock=datInd$nyr),
  
  P_Sep_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.mm[[9]],nblock=datInd$nyr),
  
  P_Oct_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.mm[[10]],nblock=datInd$nyr),
  
  P_Nov_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.mm[[11]],nblock=datInd$nyr),
  
  P_Dec_tot_m=function(data,datInd) extractor(func=get.avg.tot,data=data,indx=datInd$i.mm[[12]],nblock=datInd$nyr),
  
  Temp_JJA_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.ss[[4]],na.rm=TRUE),
  
  Temp_MAM_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.ss[[3]],na.rm=TRUE),
  
  Temp_DJF_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.ss[[2]],na.rm=TRUE),
  
  Temp_SON_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.ss[[1]],na.rm=TRUE),
  
  Temp_Jan_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[1]],nblock=datInd$nyr),
  
  Temp_Feb_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[2]],nblock=datInd$nyr),
  
  Temp_Mar_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[3]],nblock=datInd$nyr),
  
  Temp_Apr_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[4]],nblock=datInd$nyr),
  
  Temp_May_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[5]],nblock=datInd$nyr),
  
  Temp_Jun_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[6]],nblock=datInd$nyr),
  
  Temp_Jul_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[7]],nblock=datInd$nyr),
  
  Temp_Aug_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[8]],nblock=datInd$nyr),
  
  Temp_Sep_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[9]],nblock=datInd$nyr),
  
  Temp_Oct_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[10]],nblock=datInd$nyr),
  
  Temp_Nov_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[11]],nblock=datInd$nyr),
  
  Temp_Dec_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[12]],nblock=datInd$nyr),
  
  PET_JJA_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.ss[[4]],na.rm=TRUE),
  
  PET_MAM_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.ss[[3]],na.rm=TRUE),
  
  PET_DJF_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.ss[[2]],na.rm=TRUE),
  
  PET_SON_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.ss[[1]],na.rm=TRUE),
  
  PET_Jan_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[1]],nblock=datInd$nyr),
  
  PET_Feb_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[2]],nblock=datInd$nyr),
  
  PET_Mar_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[3]],nblock=datInd$nyr),
  
  PET_Apr_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[4]],nblock=datInd$nyr),
  
  PET_May_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[5]],nblock=datInd$nyr),
  
  PET_Jun_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[6]],nblock=datInd$nyr),
  
  PET_Jul_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[7]],nblock=datInd$nyr),
  
  PET_Aug_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[8]],nblock=datInd$nyr),
  
  PET_Sep_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[9]],nblock=datInd$nyr),
  
  PET_Oct_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[10]],nblock=datInd$nyr),
  
  PET_Nov_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[11]],nblock=datInd$nyr),
  
  PET_Dec_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[12]],nblock=datInd$nyr),

  Radn_ann_avg_m=function(data,datInd) extractor.summaryMean(func=mean,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,na.rm=TRUE),

  Radn_ann_rng_m=function(data,datInd) extractor.summaryMean(func=get.quantile.rng,data=data,indx=datInd$i.yy,nperiod=datInd$nyr),
  
  Radn_ann_P5_m=function(data,datInd) extractor.summaryMean(func=get.quantile,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,quant=0.05),
  
  Radn_ann_P95_m=function(data,datInd) extractor.summaryMean(func=get.quantile,data=data,indx=datInd$i.yy,nperiod=datInd$nyr,quant=0.95),
  
  Radn_JJA_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.ss[[4]],na.rm=TRUE),
  
  Radn_MAM_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.ss[[3]],na.rm=TRUE),
  
  Radn_DJF_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.ss[[2]],na.rm=TRUE),
  
  Radn_SON_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.ss[[1]],na.rm=TRUE),
  
  Radn_Jan_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[1]],nblock=datInd$nyr),
  
  Radn_Feb_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[2]],nblock=datInd$nyr),
  
  Radn_Mar_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[3]],nblock=datInd$nyr),
  
  Radn_Apr_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[4]],nblock=datInd$nyr),
  
  Radn_May_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[5]],nblock=datInd$nyr),
  
  Radn_Jun_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[6]],nblock=datInd$nyr),
  
  Radn_Jul_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[7]],nblock=datInd$nyr),
  
  Radn_Aug_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[8]],nblock=datInd$nyr),
  
  Radn_Sep_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[9]],nblock=datInd$nyr),
  
  Radn_Oct_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[10]],nblock=datInd$nyr),
  
  Radn_Nov_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[11]],nblock=datInd$nyr),
  
  Radn_Dec_avg_m=function(data,datInd) extractor(func=mean,data=data,indx=datInd$i.mm[[12]],nblock=datInd$nyr)
  
  #------------- add below
)

#ATTRIBUTE CALCULATOR FUNCTION
attribute.calculator<-function(attSel=NULL,         #list of evaluated attribute names
                               data=NULL,           #timeseries data
                               datInd=NULL,         #dat indices and properties (e.g. datInd$nyr, datInd$i.yy)
                               attribute.funcs=NULL #list of attribute calculating functions
                               ){
  
  out=lapply(attribute.funcs[attSel], function(f) f(data,datInd)) #returns labelled list of outputs
  return(out)
 
}

#ATTRIBUTE AUX INFO (determine attribute type and if approved combo with model used)
attribute.info.check<-function(attSel=NULL,  # vector of selected attributes (strings)
                               attPrim=NULL,
                               lambda.mult=NULL
                              #simVar=NULL    # vector of variables simulated using models e.g. c("P","Temp")
                              # modelTag=NULL # model selected
){
  nAtt=length(attSel) # no. of attributes nominated
  attInfo=list()      #create blank list for storage
  
  #attribute name chopper function 
  attInfo$varType=vapply(attSel,FUN = get.attribute.varType,FUN.VALUE=character(1),USE.NAMES = FALSE) #drop use of names as comes ordered anyway
  
  #ASSIGN TARGET TYPE (IF P USE "FRAC", IF T USE "DIFF")
  attInfo$targetType=vapply(attInfo$varType,FUN=get.target.type,FUN.VALUE=character(1),USE.NAMES=FALSE)
  
  #FIND WHICH ARE PRIMARY
  if(is.null(attPrim)){
    attInfo$primType=rep(FALSE,nAtt)
    attInfo$primMult=rep(0,nAtt)
  }else{
    get.ind<-function(x,y){which(x == y)}     # quick function to find which are primary attributes
    primInd=vapply(attPrim,FUN=get.ind,FUN.VALUE=numeric(1),x=attSel,USE.NAMES = FALSE)  #Indices of primary attributes
    attInfo$primType=rep(FALSE,nAtt)
    attInfo$primType[primInd]=TRUE   #mark primary ones 'TRUE'

    # Get lambda.mult (corresponding to attPrim) in the correct locations wrt attSel
    # Infomation in lamda.mult will be transferred to primMult & used in the penalty score calculation for the objective function
    # Note: attPrim is a member of attSel - but the order need not be exact, so assignment has to be done inside the loop
    attInfo$primMult=rep(0,nAtt)
    for(i in 1:length(attPrim)){
      # Locate attPrim
      indPrim <- which(attSel == attPrim[i])
      # Assign lambda.mult corresponding to attPrim
      attInfo$primMult[indPrim] <- lambda.mult[i]
    }
  }

  #CHECK FOR INVALID MODEL CHOICE - returns a logical for each attribute
  # attInfo$modelInvalid=vapply(attSel,FUN=check.attribute.model.combo,FUN.VALUE=logical(1),USE.NAMES=FALSE,modelTag=modelTag)
  
  return(attInfo)
}

get.att.ind<-function(attInfo=NULL,
                      simVar=NULL
){
  #DETERMINE WHICH ATTRIBUTE RELATES TO WHICH SIMULATOR
  attInd=list()
  if(simVar[1] != "All"){                    # ONLY DO IF STOCHASTIC GENERATION IS SELECTED (not simple scaling)
    for(i in 1:length(simVar)){
      attInd[[simVar[i]]]= which(attInfo$varType==simVar[i])
    }
  }
  return(attInd)
}



update.att.Info<-function(attInfo=NULL,
                          attInd=NULL,
                          modelTag=NULL,
                          simVar=NULL
){
  #divide up attInfo to different models
    for(i in 1:length(modelTag)){
      attInfo[[modelTag[i]]]$varType=attInfo$varType[attInd[[simVar[i]]]]
      attInfo[[modelTag[i]]]$targetType=attInfo$targetType[attInd[[simVar[i]]]]
      attInfo[[modelTag[i]]]$primType=attInfo$primType[attInd[[simVar[i]]]]
      attInfo[[modelTag[i]]]$primMult=attInfo$primMult[attInd[[simVar[i]]]]
    }
  return(attInfo)
}

#GETS VARTYPE BY READING FIRST ELEMENT OF ATTRIBUTE STRING
get.attribute.varType<-function(attrib=NULL, # attribute name
                                 sep="_"){
  varType=strsplit(x = attrib,split=sep)[[1]][1]
  return(varType)
}

#get.attribute.varType(attrib=attSel[1], sep="_")

#TARGET TYPE CLASSIFIER
get.target.type<-function(varType=NULL){
  if(varType == "P"){
    targetType="frac"
  }else{
      if(varType == "Temp"){
        targetType="diff"
      }else{
        targetType="frac"
      }
  }
  return(targetType)
}

tagBlender<-function(attLab=NULL
){
  
  chopped=strsplit(x = attLab,split="_")[[1]]
  
  #variable type
  if(chopped[1]== "P"){
    vtype="rainfall (fraction)"
  }else if(chopped[1]== "Temp"){
    vtype="temperature (additive change)"
  }else if(chopped[1]== "PET"){
    vtype="PET (fraction)"
  }else if(chopped[1]=="Radn"){
    vtype="Radn (fraction)"
  }
  
  #aggregation type
  if(chopped[2]== "ann"){
    atype="annual"
  }else if(chopped[2]== "JJA"){
    atype="JJA"
  }else if(chopped[2]== "MAM"){
    atype="MAM"
  }else if(chopped[2]== "DJF"){
    atype="DJF"
  }else if(chopped[2]== "SON"){
    atype="SON"
  }else if(chopped[2]== "Jan"){
    atype="Jan"
  }else if(chopped[2]== "Feb"){
    atype="Feb"
  }else if(chopped[2]== "Mar"){
    atype="Mar" 
  }else if(chopped[2]== "Apr"){
    atype="Apr"     
  }else if(chopped[2]== "May"){
    atype="May" 
  }else if(chopped[2]== "Jun"){
    atype="Jun"  
  }else if(chopped[2]== "Jul"){
    atype="Jul"  
  }else if(chopped[2]== "Aug"){
    atype="Aug"      
  }else if(chopped[2]== "Sep"){
    atype="Sep"  
  }else if(chopped[2]== "Oct"){
    atype="Oct"
  }else if(chopped[2]== "Nov"){
    atype="Nov"
  }else if(chopped[2]== "Dec"){
    atype="Dec"
  }
  
  #metricType
  if(chopped[3]== "nWet"){
    mtype="no. wet days"
  }else if(chopped[3]== "dyWet"){
    mtype="wet day amount"
  }else if(chopped[3]== "DSD"){
    mtype="dryspell duration"
  }else if(chopped[3]== "P99"){
    mtype="99th percentile day amount"
  }else if(chopped[3]== "dyWet99p"){
    mtype="99th percentile wet day amount"
  }else if(chopped[3]== "avgWSD"){
    mtype="average wetspell duration"
  }else if(chopped[3]== "avgDSD"){
    mtype="average dryspell duration"
  }else if(chopped[3]== "maxDSD"){
    mtype="max dryspell duration"
  }else if(chopped[3]== "maxWSD"){
    mtype="max wetspell duration"
  }else if(chopped[3]== "tot"){
    mtype="total"
  }else if(chopped[3]== "R10"){
    mtype="no. days above 10mm"
  }else if(chopped[3]== "GSL"){
    mtype="growing season length"
  }else if(chopped[3]== "CSL"){
    mtype="cold season length"
  }else if(chopped[3]== "avg"){
    mtype="average"
  }else if(chopped[3]== "P5"){
    mtype="5th percentile"
  }else if(chopped[3]== "P95"){
    mtype="95th percentile"
  }else if(chopped[3]== "F0"){
    mtype="frost days"
  }else if(chopped[3]== "rng"){
    mtype="range"
  }else if(chopped[3]== "90pX"){
    mtype="percent above historical 90th percentile"
  }else if(chopped[3]== "90X"){
    mtype="volume above historical 90th percentile"
  }else if(chopped[3]== "seasRatio"){
    mtype="ratio of wet to dry seasonal volume"
  }
  
  #statType
  if(chopped[4]== "m"){
    stype="Mean"  #as yet un-used
  }
  
  #stitch togther
  phrase=paste(stype,atype,mtype,vtype)
  phrase  
}


tagBlender_noUnits<-function(attLab=NULL
){
  
  chopped=strsplit(x = attLab,split="_")[[1]]
  
  #variable type
  if(chopped[1]== "P"){
    vtype="rainfall"
  }else if(chopped[1]== "Temp"){
    vtype="temperature"
  }else if(chopped[1]== "PET"){
    vtype="PET"
  }else if(chopped[1]=="Radn"){
    vtype="Radn"
  }
  
  #aggregation type
  if(chopped[2]== "ann"){
    atype="annual"
  }else if(chopped[2]== "JJA"){
    atype="JJA"
  }else if(chopped[2]== "MAM"){
    atype="MAM"
  }else if(chopped[2]== "DJF"){
    atype="DJF"
  }else if(chopped[2]== "SON"){
    atype="SON"
  }else if(chopped[2]== "Jan"){
    atype="Jan"
  }else if(chopped[2]== "Feb"){
    atype="Feb"
  }else if(chopped[2]== "Mar"){
    atype="Mar" 
  }else if(chopped[2]== "Apr"){
    atype="Apr"     
  }else if(chopped[2]== "May"){
    atype="May" 
  }else if(chopped[2]== "Jun"){
    atype="Jun"  
  }else if(chopped[2]== "Jul"){
    atype="Jul"  
  }else if(chopped[2]== "Aug"){
    atype="Aug"      
  }else if(chopped[2]== "Sep"){
    atype="Sep"  
  }else if(chopped[2]== "Oct"){
    atype="Oct"
  }else if(chopped[2]== "Nov"){
    atype="Nov"
  }else if(chopped[2]== "Dec"){
    atype="Dec"
  }
  
  #metricType
  # Do not need to add rainfall to the end
  # Added "rainfall" already here for the P99 indices & dyWet
  if(chopped[3]== "nWet"){
    mtype="no. wet days"
  }else if(chopped[3]== "dyWet"){
    mtype="wet day rainfall"
  }else if(chopped[3]== "DSD"){
    mtype="dryspell duration"
  }else if(chopped[3]== "P99"){
    mtype="99th percentile rainfall"
  }else if(chopped[3]== "dyWet99p"){
    mtype="99th percentile wet day rainfall"
  }else if(chopped[3]== "avgWSD"){
    mtype="average wetspell duration"
  }else if(chopped[3]== "avgDSD"){
    mtype="average dryspell duration"
  }else if(chopped[3]== "maxDSD"){
    mtype="max dryspell duration"
  }else if(chopped[3]== "maxWSD"){
    mtype="max wetspell duration"
  
  #*********NEED to add vtype 
  }else if(chopped[3]== "tot"){
    mtype="total"
  }else if(chopped[3]== "R10"){
    mtype="no. days above 10mm"
  
  # do not need
  }else if(chopped[3]== "GSL"){
    mtype="growing season length"
  }else if(chopped[3]== "CSL"){
    mtype="cold season length"
  
  #********NEED to add vtype
  }else if(chopped[3]== "avg"){
    mtype="average"
  }else if(chopped[3]== "P5"){
    mtype="5th percentile"
  }else if(chopped[3]== "P95"){
    mtype="95th percentile"
  
  # do not need
  }else if(chopped[3]== "F0"){
    mtype="no. of frost days"
  }else if(chopped[3]== "rng"){
    mtype="range"
  
  # Not used currently
  }else if(chopped[3]== "90pX"){
    mtype="percent above historical 90th percentile"
  }else if(chopped[3]== "90X"){
    mtype="volume above historical 90th percentile"
  
  # Already added "rainfall" here
  }else if(chopped[3]== "seasRatio"){
    mtype="ratio of wet to dry season rainfall"
  }
  
  #statType
  if(chopped[4]== "m"){
    stype="Mean"  #as yet un-used
  }
  
  #stitch togther
  if (chopped[3] %in% c("tot", "R10", "avg", "P5", "P95", "rng")) {
    phrase=paste(stype,atype,mtype,vtype)
  } else {
    phrase=paste(stype,atype,mtype)
  }
  return(phrase)  
}


