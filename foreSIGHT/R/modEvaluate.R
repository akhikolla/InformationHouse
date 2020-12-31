# #Created by Sam Culley 14/11/2019
# 
# #Functions for evaluating WGEN models fit to observed data
# 
# modEvaluate<-function(data,path){
#  
#   plots<-compare.WGEN.SILO(data=data,path=path,col=4,ylim=c(300,1100),var="Total P (mm)",ylim2=c(0,150),ylim3=c(0,300),ylim4=c(0,30))
#   ggsave("HistWGENP.pdf",plots[[1]],width=12,height=10)
#   ggsave("HistWGENPMonthly.pdf",plots[[2]],width=12,height=10)
#   ggsave("HistWGENNWET.pdf",plots[[3]],width=12,height=10)
#   ggsave("HistWGENNWETMonthly.pdf",plots[[4]],width=12,height=10)
# 
#   plots2<-compare.WGEN.SILO(data=data,path=path,col=5,ylim=c(1100,1600),var="Total PET (mm)",ylim2=c(0,300))
#   ggsave("HistWGENPET.pdf",plots2[[1]],width=12,height=10)
#   ggsave("HistWGENPETMonthly.pdf",plots2[[2]],width=12,height=10)
#   
#   plots3<-compare.WGEN.SILO.sd(data=data,path=path,col=4,ylim=c(0,10),var="P StDev (mm)",ylim2=c(0,10),ylim3=c(0,1),ylim4=c(0,1))
#   ggsave("HistWGENPsd.pdf",plots3[[1]],width=12,height=10)
#   ggsave("HistWGENPMonthlysd.pdf",plots3[[2]],width=12,height=10)
#   ggsave("HistWGENNWETsd.pdf",plots3[[3]],width=12,height=10)
#   ggsave("HistWGENNWETMonthlysd.pdf",plots3[[4]],width=12,height=10)
# 
#   plots4<-compare.WGEN.SILO.sd(data=data,path=path,col=5,ylim=c(0,5),var="PET StDev (mm)",ylim2=c(0,5))
#   ggsave("HistWGENPETsd.pdf",plots4[[1]],width=12,height=10)
#   ggsave("HistWGENPETMonthlysd.pdf",plots4[[2]],width=12,height=10)
# 
# }
# 
# 
# compare.WGEN.SILO<-function(data,path,col,ylim,var,ylim2,ylim3=NULL,ylim4=NULL){
#   
#   
#   AnnualPRaw<-data.frame(as.numeric(data$year),as.numeric(data$month),data[col])
#   names(AnnualPRaw)<-c("Year","Month","Data")
#   if(col==4 || col==5){
#     AnnualSummary<-aggregate(x=AnnualPRaw$Data,by=list(AnnualPRaw$Year),FUN=sum)[2]
#   }
#   AnnualPlot<-data.frame(seq(1961,2005),AnnualSummary)
#   names(AnnualPlot)<-c("Year","DataV")
# 
#   if(col==4 || col==5){
#     AnnualSummary<-aggregate(x=AnnualPRaw$Data,by=list(AnnualPRaw$Month),FUN=sum)[2]
#   }
#   MonthlyPlot<-data.frame(seq(1,12),AnnualSummary/45)
#   names(MonthlyPlot)<-c("Month","DataV")
#   
#   p1<-forward.WGEN.plot(path=path,obs=AnnualPlot,col=col,ylim=ylim,var=var,func="sum",nwet=FALSE)
#   
#   p2<-forward.WGEN.plot.monthly(path=path,obs=MonthlyPlot,col=col,ylim=ylim2,var=var,func="sum",nwet=FALSE)
#   
#   if(var=="Total P (mm)" || var =="P StDev (mm)"){
#     temp<-data[4]
#     temp[temp>0]<-1
#     AnnualPRaw<-data.frame(as.numeric(data$year),as.numeric(data$month),temp)
#     names(AnnualPRaw)<-c("Year","Month","Data")
#     AnnualSummary<-aggregate(x=AnnualPRaw$Data,by=list(AnnualPRaw$Year),FUN=sum)[2]
#     AnnualPlot<-data.frame(seq(1961,2005),AnnualSummary)
#     names(AnnualPlot)<-c("Year","DataV")
#     AnnualSummary<-aggregate(x=AnnualPRaw$Data,by=list(AnnualPRaw$Month),FUN=sum)[2]
#     MonthlyPlot<-data.frame(seq(1,12),AnnualSummary/45)
#     names(MonthlyPlot)<-c("Month","DataV")
#     
#     p3<-forward.WGEN.plot(path=path,obs=AnnualPlot,col=col,ylim=ylim3,var="Number of Wet Days",func="sum",nwet=TRUE)
#     
#     p4<-forward.WGEN.plot.monthly(path=path,obs=MonthlyPlot,col=col,ylim=ylim4,var="Number of Wet Days",func="sum",nwet=TRUE)
#     
#     out<-list()
#     out[[1]]<-p1
#     out[[2]]<-p2
#     out[[3]]<-p3
#     out[[4]]<-p4
#   } else {
#     out<-list()
#     out[[1]]<-p1
#     out[[2]]<-p2
#   }
#   return(out)
# }
# 
# compare.WGEN.SILO.sd<-function(data,path,col,ylim,var,ylim2,ylim3=NULL,ylim4=NULL){
#   
#   
#   AnnualPRaw<-data.frame(as.numeric(data$year),as.numeric(data$month),data[col])
#   names(AnnualPRaw)<-c("Year","Month","Data")
#   if(col==4 || col==5){
#     AnnualSummary<-aggregate(x=AnnualPRaw$Data,by=list(AnnualPRaw$Year),FUN=sd)[2]
#   }
#   AnnualPlot<-data.frame(seq(1961,2005),AnnualSummary)
#   names(AnnualPlot)<-c("Year","DataV")
#   
#   if(col==4 || col==5){
#     AnnualSummary<-aggregate(x=AnnualPRaw$Data,by=list(AnnualPRaw$Month),FUN=sd)[2]
#   }
#   MonthlyPlot<-data.frame(seq(1,12),AnnualSummary)
#   names(MonthlyPlot)<-c("Month","DataV")
#   
#   p1<-forward.WGEN.plot(path=path,obs=AnnualPlot,col=col,ylim=ylim,var=var,func="sd",nwet=FALSE)
#   
#   p2<-forward.WGEN.plot.monthly(path=path,obs=MonthlyPlot,col=col,ylim=ylim2,var=var,func="sd",nwet=FALSE)
#   
#   if(var=="P StDev (mm)"){
#     temp<-data[4]
#     temp[temp>0]<-1
#     AnnualPRaw<-data.frame(as.numeric(data$year),as.numeric(data$month),temp)
#     names(AnnualPRaw)<-c("Year","Month","Data")
#     AnnualSummary<-aggregate(x=AnnualPRaw$Data,by=list(AnnualPRaw$Year),FUN=sd)[2]
#     AnnualPlot<-data.frame(seq(1961,2005),AnnualSummary)
#     names(AnnualPlot)<-c("Year","DataV")
#     AnnualSummary<-aggregate(x=AnnualPRaw$Data,by=list(AnnualPRaw$Month),FUN=sd)[2]
#     MonthlyPlot<-data.frame(seq(1,12),AnnualSummary)
#     names(MonthlyPlot)<-c("Month","DataV")
#     
#     p3<-forward.WGEN.plot(path=path,obs=AnnualPlot,col=col,ylim=ylim3,var="Number of Wet Days",func="sd",nwet=TRUE)
#     
#     p4<-forward.WGEN.plot.monthly(path=path,obs=MonthlyPlot,col=col,ylim=ylim4,var="Number of Wet Days",func="sd",nwet=TRUE)
#     
#     out<-list()
#     out[[1]]<-p1
#     out[[2]]<-p2
#     out[[3]]<-p3
#     out[[4]]<-p4
#   } else {
#     out<-list()
#     out[[1]]<-p1
#     out[[2]]<-p2
#   }
#   return(out)
# }
# 
# forward.WGEN.plot<-function(path,obs,col,ylim,var,func,nwet){
#   Masterdataplot<-data.frame(matrix(ncol=3,nrow=0))
#   names(Masterdataplot)<-c("Year","Data","Replicate")
#   temp<-seq(1:100)
#   rep<-sprintf("%03d",temp)
#   for(i in 1:100){
#     print(i)
#     file<-paste0(path,"out",i,".csv")
#     
#     data<-read.csv(file,header=TRUE)
#     names(data)<-c("Year","Month","Day","P","PET")
#     if(nwet==TRUE){
#       temp<-data[4]
#       temp[temp>0]<-1
#       tempdata<-data.frame(data$Year,temp)
#       names(tempdata)<-c("Year","Data")
#     } else {
#       tempdata<-data.frame(data$Year,data[col])
#       names(tempdata)<-c("Year","Data")
#     }
#     
#     if(func=="sum"){
#     AnnualSummary<-aggregate(x=tempdata$Data,by=list(tempdata$Year),FUN=sum)[2]
#     } else if (func=="sd"){
#     AnnualSummary<-aggregate(x=tempdata$Data,by=list(tempdata$Year),FUN=sd)[2]
#     }
#     dataplot<-data.frame(seq(1961,2005),AnnualSummary,rep(rep[i],length(AnnualSummary)))
#     names(dataplot)<-c("Year","Data","Replicate")
#     Masterdataplot<-rbind(Masterdataplot,dataplot)
#     #create data frame
#     
#   }
#   p1<-ggplot()+
#     geom_line(data=Masterdataplot,aes(x=Year,y=Data,group=Replicate),size=0.5,alpha=0.1)+
#     ylim(ylim[1],ylim[2])+
#     geom_line(data=obs,aes(x=Year,y=DataV),colour="#F8766D",size=1.5)+
#     labs(y=var)
#   return(p1)
# }
# 
# forward.WGEN.plot.monthly<-function(path,obs,col,ylim,var,func,nwet){
#   Masterdataplot<-data.frame(matrix(ncol=3,nrow=0))
#   names(Masterdataplot)<-c("Month","Data","Replicate")
#   temp<-seq(1:100)
#   rep<-sprintf("%03d",temp)
#   for(i in 1:100){
#     print(i)
#     file<-paste0(path,"out",i,".csv")
#     
#     data<-read.csv(file,header=TRUE)
#     names(data)<-c("Year","Month","Day","P","PET")
#     if(nwet==TRUE){
#       temp<-data[4]
#       temp[temp>0]<-1
#       tempdata<-data.frame(data$Month,temp)
#       names(tempdata)<-c("Month","Data")
#     } else {
#       tempdata<-data.frame(data$Month,data[col])
#       names(tempdata)<-c("Month","Data")
#     }
#     
#     
#     if(func=="sum"){
#       AnnualSummary<-(aggregate(x=tempdata$Data,by=list(tempdata$Month),FUN=sum)[2])/45
#     } else if (func=="sd"){
#       AnnualSummary<-aggregate(x=tempdata$Data,by=list(tempdata$Month),FUN=sd)[2]
#     }
#     dataplot<-data.frame(seq(1,12),AnnualSummary,rep(rep[i],length(AnnualSummary)))
#     
#     
#     names(dataplot)<-c("Month","Data","Replicate")
#     Masterdataplot<-rbind(Masterdataplot,dataplot)
#     #create data frame
#     
#   }
#   p1<-ggplot()+
#     geom_line(data=Masterdataplot,aes(x=Month,y=Data,group=Replicate),size=0.5,alpha=0.1)+
#     ylim(ylim[1],ylim[2])+
#     geom_line(data=obs,aes(x=Month,y=DataV),colour="#F8766D",size=1.5)+
#     scale_x_discrete(limits=seq(1,12),labels=c("1","2","3","4","5","6","7","8","9","10","11","12"))+
#     labs(y=var)
#   return(p1)
# }





# modEvaluate2<-function(){
#   
#   
#   inputcheck<-input_check(obs,file,simLengthNyrs)
#   obs=inputcheck$data                                      # USE NEW APPENDED/CHECKED DATA
#   
#   
#   #GET ADDITIONAL MODEL INFO, ATT INFO & SORT (make into separate script/functions)
#   nMod=length(modelTag)   
#   modelInfo=get.multi.model.info(modelTag=modelTag)
#   
#   #UPDATE MODELINFO IF NEEDED
#   for(mod in 1:nMod){
#     if(!is.null(modelInfoMod[[modelTag[mod]]])){
#       #modifyList
#       if(mod==1) progress("Updating model info...",file)
#       defaultMods=list(minBound=NULL,maxBound=NULL,fixedPars=NULL)
#       modPars=modifyList(defaultMods,modelInfoMod[[modelTag[mod]]])
#       modelInfo[[modelTag[mod]]]=update.model.info(modelTag=modelTag[mod],
#                                                    modelInfo=modelInfo[[modelTag[mod]]],
#                                                    fixedPars=modPars$fixedPars,
#                                                    minUserBound=modPars$minBound,
#                                                    maxUserBound=modPars$maxBound,
#                                                    file=file)  #need to build in checks for this
#       # if(!is.na(modelInfo[[modelTag[mod]]]$fixedPars)
#     }
#   }
#   
#   modelTag=update.simPriority(modelInfo=modelInfo)
#   simVar=sapply(X=modelInfo[modelTag],FUN=return.simVar,USE.NAMES=TRUE)       #?CREATE MODEL MASTER INFO - HIGHER LEVEL?
#   attInfo=attribute.info.check(attSel=attSel,attPrim=attPrim)                                 # vector of selected attributes (strings)
#   if(modelTag[1] == "Simple-ann"){simVar=attInfo$varType}
#   attInd=get.att.ind(attInfo=attInfo,simVar=simVar)
#   attInfo=update.att.Info(attInfo=attInfo,attInd=attInd,modelTag=modelTag,simVar=simVar) #add extra level for easier model mangmt
#   if(modelTag[1] != "Simple-ann"){nParTot=0;for(i in 1:nMod){nParTot=nParTot+modelInfo[[i]]$npars}}                      #total number of pars
#   
#   
#   #GET DATES DATA (and indexes for harmonic periods)
#   datInd=mod.get.date.ind(obs=obs,modelTag=modelTag,modelInfo=modelInfo,southHemi=TRUE)
#   
#   
#   devPlotSummaryReps
#   
# }
# 
# 
# devPlotSummaryReps<-function(obs=NULL,
#                              sim=NULL,
#                              simVar=NULL,
#                              datInd=NULL,
#                              modelTag=NULL,
#                              IOmode=NULL
#                             # paths=paths
# ){
#   
#   
#   
#   
#   
#   # DO EVERYTHING FOR OBSERVED SERIES ONCE
#   obsDat=list()
#   for(i in 1:length(simVar)){
#     plotVar=simVar[i]            #select variable to evaluate
#     obsDat[[plotVar]]=collateDat(TS=obs[[plotVar]],datInd=datInd[["obs"]],plotVar=plotVar)
#   }
#   
#   #LOOP OVER TARGETS
#   simDat=list()
#   
#   
#   #LABEL A PDF
#   fnam="WGENcheck.pdf"
#   #PLOT STUFF TO A PDF
#   pdf(file=fnam,height=8.27,width=11.69)   #landscape a4 page
#   par(mar=c(3,5,3,3),oma=c(2,2,2,2))
#   
#   #FRONT BOILERPLATE INFO
#   
#   #TRAFFIC LIGHT PLOT HERE
#   
#   #SET LAYOUT - 2 ROWS, 1 COLUMN
#   
#   for(mod in 1:length(simVar)){
#     
#     plotVar=simVar[mod]            #select variable to evaluate
#     simTest=sim[[plotVar]]$sim
#     
#     simDat[[plotVar]]=collateDat(TS=simTest,datInd=datInd[[modelTag[mod]]],plotVar=plotVar)
#     
#     #print(obsDat[[plotVar]][["ann_count_nWet"]])
#     # print(simDat[[plotVar]][["ann_count_nWet"]])
#     
#     #SET LAYOUT - 2 ROWS, 1 COLUMN
#     par(mfrow=c(2,1),xaxs="i",xpd=FALSE)
#     par(mar=c(3,5,3,3),oma=c(3,5,3,3))
#     #REGULAR PLOTS
#     #monthwise batch
#     runTag="mon_mean_dyAll"; lab=paste(plotVar,": daily mean", sep="")
#     monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#     
#     runTag="mon_sd_dyAll"; lab=paste(plotVar,": daily sd", sep="")
#     monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#     
#     runTag="mon_sum_dyAll"; lab=paste(plotVar,": total", sep="")
#     monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#     
#     if(plotVar == "P"){
#       runTag="mon_mean_dyWet"; lab=paste(plotVar,": daily wet mean", sep="")
#       monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       
#       runTag="mon_sd_dyWet"; lab=paste(plotVar,": daily wet sd", sep="")
#       monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       
#       # runTag="mon_mean_nWet"; lab=paste(plotVar,": no. wet days mean", sep="")
#       # monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       
#       runTag="mon_count_nWet"; lab=paste(plotVar,": no. wet days", sep="")
#       monthwise.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#     }
#     
#     
#     #seasonal batch
#     runTag="seas_mean_dyAll"; lab=paste(plotVar,": daily mean", sep="")
#     seasonal.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#     
#     runTag="seas_sd_dyAll"; lab=paste(plotVar,": daily sd", sep="")
#     seasonal.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#     
#     runTag="seas_sum_dyAll"; lab=paste(plotVar,": total", sep="")
#     seasonal.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#     
#     #annual batch
#     #change
#     par(mfrow=c(2,3),xaxs="i",xpd=FALSE)  #assuming a4 landscape layout
#     par(mar=c(3,5,3,3),oma=c(3,5,3,3))
#     
#     runTag="ann_mean_dyAll"; lab=paste(plotVar,": daily mean", sep="")
#     annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#     
#     runTag="ann_sd_dyAll"; lab=paste(plotVar,": daily sd", sep="")
#     annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#     
#     runTag="ann_sum_dyAll"; lab=paste(plotVar,": total", sep="")
#     annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#     
#     if(plotVar == "P"){
#       runTag="ann_mean_dyWet"; lab=paste(plotVar,": daily wet mean", sep="")
#       annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       runTag="ann_sd_dyWet"; lab=paste(plotVar,": daily wet sd", sep="")
#       annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       # runTag="ann_mean_nWet"; lab=paste(plotVar,": no. wet days mean", sep="")
#       # annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       runTag="ann_count_nWet"; lab=paste(plotVar,": no. wet days", sep="")
#       annual.boxplots(simDat=simDat[[plotVar]][[runTag]],obsDat=obsDat[[plotVar]][[runTag]],compObs=TRUE, metricTag=lab)
#       
#     }
#     
#   }
#   dev.off()  #STOP PLOTTING TO PDF
#   
#   
# }