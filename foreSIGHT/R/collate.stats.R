#################################################
##  MAKE THESE  GENERIC  EXTRACTOR FUNCTIONS   ##
#################################################

#CONTAINS
  #CollateDat - wrapper for extract.funcs
  #CollateDatReps - wrapper for extract.funcs when comparing across replicates
  #stat.func - list of functions
  #funcSel - vector of selected functions
  #monthwise.extract.func()
  #annual.extract.func()
  #seasonal.extract.func()


collateDat<-function(TS=NULL,
                     datInd=NULL,
                     plotVar=NULL
){
  statTag=c("mean","sd","max","min","sum","count")
  aggTag=c("mon","ann","seas")
  metricTag=c("dyAll","dyWet","nWet")  
  
  outDat=list()
  
  #GENERAL STATISTICS
  k=1
  for(i in 1:3){
    for(j in 1:5){
      runTag=paste(aggTag[i],statTag[j],metricTag[k],sep="_")  #will need amending
      evalFunc=eval(parse(text = statTag[j]))
      switch(aggTag[i],
             "mon" = {outDat[[runTag]]=monthwise.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],na.rm=TRUE)},
             "ann" = {outDat[[runTag]]=annual.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],na.rm=TRUE)},
             "seas" = {outDat[[runTag]]=seasonal.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],na.rm=TRUE)},
             -999.00 )
    }  #stat loop
  } #agg loop
  
  #RAINFALL SPECIFIC STATISTICS
  if(plotVar == "P"){   #only if rainfall
    k=2 
    for(i in 1:3){
      for(j in 1:2){
        runTag=paste(aggTag[i],statTag[j],metricTag[k],sep="_")  #will need amending
        switch(statTag[j],
               "mean" = {evalFunc=get.wet.average},
               "sd"  = {evalFunc=get.wet.sd},
               {evalFunc=eval(parse(text = statTag[j]))}
        )
        switch(aggTag[i],
               "mon" = {outDat[[runTag]]=monthwise.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],threshold=0)},
               "ann" = {outDat[[runTag]]=annual.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],threshold=0)},
               "seas" = {outDat[[runTag]]=seasonal.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],threshold=0)},
               -999.00 )
      }  #stat loop
    } #agg loop
    
    k=3; j=6
    for(i in 1:3){
      runTag=paste(aggTag[i],statTag[j],metricTag[k],sep="_")  #will need amending
      evalFunc=get.nwet
      switch(aggTag[i],
             "mon" = {outDat[[runTag]]=monthwise.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],threshold=0)},
             "ann" = {outDat[[runTag]]=annual.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],threshold=0)},
             "seas" = {outDat[[runTag]]=seasonal.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],threshold=0)},
             -999.00 )
    }
    
  } #end rain loop
  return(outDat)
}

collateDatReps<-function(TS=NULL,
                         datInd=NULL,
                         plotVar=NULL
){
  statTag=c("mean","sd","max","min","sum","count")
  aggTag=c("mon","ann","seas")
  metricTag=c("dyAll","dyWet","nWet")  
  
  outDat=list()
  
  #GENERAL STATISTICS
  k=1
  for(i in 1:3){
    for(j in 1:5){
      runTag=paste(aggTag[i],statTag[j],metricTag[k],sep="_")  #will need amending
      evalFunc=eval(parse(text = statTag[j]))
      switch(aggTag[i],
             "mon" = {outDat[[runTag]]=monthwise.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],na.rm=TRUE)},
             "ann" = {outDat[[runTag]]=annual.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],na.rm=TRUE)},
             "seas" = {outDat[[runTag]]=seasonal.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],na.rm=TRUE)},
             -999.00 )
    }  #stat loop
  } #agg loop
  
  #RAINFALL SPECIFIC STATISTICS
  if(plotVar == "P"){   #only if rainfall
    k=2 
    for(i in 1:3){
      for(j in 1:2){
        runTag=paste(aggTag[i],statTag[j],metricTag[k],sep="_")  #will need amending
        switch(statTag[j],
               "mean" = {evalFunc=get.wet.average},
               "sd"  = {evalFunc=get.wet.sd},
               {evalFunc=eval(parse(text = statTag[j]))}
        )
        switch(aggTag[i],
               "mon" = {outDat[[runTag]]=monthwise.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],threshold=0)},
               "ann" = {outDat[[runTag]]=annual.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],threshold=0)},
               "seas" = {outDat[[runTag]]=seasonal.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],threshold=0)},
               -999.00 )
      }  #stat loop
    } #agg loop
    
    k=3; j=6
    for(i in 1:3){
      runTag=paste(aggTag[i],statTag[j],metricTag[k],sep="_")  #will need amending
      evalFunc=get.nwet
      switch(aggTag[i],
             "mon" = {outDat[[runTag]]=monthwise.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],threshold=0)},
             "ann" = {outDat[[runTag]]=annual.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],threshold=0)},
             "seas" = {outDat[[runTag]]=seasonal.extract.func(f=evalFunc,dat=TS,datInd=datInd,stat.func=stat.func,funcSel=funcSel[1:5],threshold=0)},
             -999.00 )
    }
    
  } #end rain loop
  return(outDat)
}
#--------------------------------------------------------------------------------
#WRITE WRAPPER FUNC THAT USES EACH OF THESE 
  #direct different tags to different functions...

# cat=c("tot","max","nwet","dsd","wsd","mean","sd","nwetsd","skew","skewWet")
# subcat=c("ann","month","seas","all")

#Function lists to store
stat.func<-list(max=function(data) max(x=data,na.rm=TRUE),
                min=function(data) min(x=data,na.rm=TRUE),
                median=function(data) median(x=data,na.rm=TRUE),
                mean=function(data) mean(x=data,na.rm=TRUE),
                sd=function(data) sd(x=data,na.rm=TRUE),
                skew=function(data) skewness(x=data,na.rm=TRUE)
                # p5=function(data) quantile(x=data,quant=0.05,na.rm=TRUE),  #maybe optional
                # p95=function(data) quantile(x=data,quant=0.95,na.rm=TRUE)  #maybe optional
                #wetspell dryspell function versions
)

funcSel=c("max","min","median","mean","sd","skew")  #to match stat.func/and make some optional
#----------------------------------------------------------------------------------------------
#PART 1: GET STATISTICS ON A MONTHLY DIVISION (should be a generic function style with on off toggles)
monthwise.extract.func=function(f,       #function by which data will be extracted sum, mean, max, min etc
                                dat=NULL,
                                datInd=NULL,
                                stat.func=NULL, # whole stat.func list
                                funcSel=NULL,   # selected functions
                                ...
    ){
  
  
  #may need some on/off toggles here
  nMonths=12*(datInd$nyr)
  
  sim.month=rep(0,nMonths); trans.ind=rep(0,nMonths)
  counter=0                                                # zero month counter
  for(Y in 1:datInd$nyr){
    for(m in 1:12){
      counter=counter+1                                    # move to next month
      ind=intersect(x=datInd$i.mm[[m]],y=datInd$i.yy[[Y]]) # extract indices for which the months in the examined year correspond
      trans.ind[counter]=m
      sim.month[counter]=extractor(func=f,data=dat,indx=ind,...)
    }
  }
  
  #SORT INTO MONTHLY ARRAYS TO STORE
  sim.month.sort=matrix(0,nrow=12,ncol=datInd$nyr)
  for(m in 1:12){  
    ind=which(trans.ind==m)
    sim.month.sort[m,]=sim.month[ind]    #stacked in rows
  }
  
  #GRAB INDICATOR STATS (MEAN, MAX, MIN, MED, SKEW,SD)   
   short.func<-function(f,dat){apply(dat,MARGIN=1,FUN=f)}   # create local apply function
   out=lapply(stat.func[funcSel],function(f) short.func(f,dat=sim.month.sort)) #get stats for each month division
  
  # COR1'S REQUIRE DIFF APPROACH for storage(potentially)

  #COMPILE OUT LIST
  out$TS=sim.month.sort
  #already contains max,min etc per stat.funcs
  return(out)
}
#Tester
#tmp=monthwise.extract.func(f=mean,dat=obs,datInd=datInd,stat.func=stat.func,funcSel=funcSel)


#PART 2: GET STATISTICS ON AN ANNUAL DIVISION
annual.extract.func=function(f,       #function by which data will be extracted sum, mean, max, min etc
                             dat=NULL,
                             datInd=NULL,
                             stat.func=NULL, # whole stat.func list
                             funcSel=NULL,   # selected functions
                             ...
){

  sim.ann=rep(0,datInd$nyr)
  for(Y in 1:datInd$nyr){
      sim.ann[Y]=extractor(func=f,data=dat,indx=datInd$i.yy[[Y]],...)
  }
  
  #GRAB INDICATOR STATS (MEAN, MAX, MIN, MED, SKEW,SD, etc as needed)   
  out=lapply(stat.func[funcSel],function(f) f(sim.ann)) #get stats for each month division

  # COR1'S REQUIRE DIFF APPROACH
  
  #ADD TS to out list
  out$TS=sim.ann
  
  return(out)
}
#Tester
#tmp=annual.extract.func(f=mean,dat=obs,datInd=datInd,stat.func=stat.func,funcSel=funcSel)

#PART 3: GET STATISTICS ON A SEASONAL DIVISION (similar to monthly set... or annual...)
seasonal.extract.func=function(f,       #function by which data will be extracted sum, mean, max, min etc
                                dat=NULL,
                                datInd=NULL,
                                stat.func=NULL, # whole stat.func list
                                funcSel=NULL,   # selected functions
                                ...
){

  sim.seas=rep(0,4*datInd$nyr); trans.ind=rep(0,4*datInd$nyr)
  counter=0                                                # zero month counter
  for(Y in 1:datInd$nyr){
    for(s in 1:4){
      counter=counter+1                                    # move to next month
      ind=intersect(x=datInd$i.ss[[s]],y=datInd$i.yy[[Y]]) # extract indices for which the months in the examined year correspond
      trans.ind[counter]=s
      sim.seas[counter]=extractor(func=f,data=dat,indx=ind,...)
    }
  }
  
  #SORT INTO MONTHLY ARRAYS TO STORE
  sim.seas.sort=matrix(0,nrow=4,ncol=datInd$nyr)
  for(s in 1:4){  
    ind=which(trans.ind==s)
    sim.seas.sort[s,]=sim.seas[ind]    #stacked in rows
  }
  
  #GRAB INDICATOR STATS (MEAN, MAX, MIN, MED, SKEW,SD)   
  short.func<-function(f,dat){apply(dat,MARGIN=1,FUN=f)}   # create local apply function
  out=lapply(stat.func[funcSel],function(f) short.func(f,dat=sim.seas.sort)) #get stats for each month division
  
  # COR1'S REQUIRE DIFF APPROACH for storage(potentially)
  
  #COMPILE OUT LIST
  out$TS=sim.seas.sort
  #already contains max,min etc per stat.funcs
  return(out)
}
#Tester
#tmp=seasonal.extract.func(f=mean,dat=obs,datInd=datInd,stat.func=stat.func,funcSel=funcSel)
