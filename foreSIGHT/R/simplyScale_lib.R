#######################################
##   SIMPLE SCALING FUNCTIONS     ##
#######################################

# CONTAINS
  #function for scaling over a period

simple.scaling<-function(target=NULL,          #extracted from matrix output of exposure space maker. annual format is a vector, passed from wrapper.
                         targetType=NULL,      #diff or frac
                         data=NULL,            #data series stripped of dates
                         varType=NULL,         #names of data series
                         period=NULL,          #number of separations for scaling e.g. annual all at once, seasonal in four stages.
                         i.pp=NULL            #index to control what time step entries are changed at a time
) {

  temp=matrix(NA,nrow=nrow(data),ncol=ncol(data)) #create matrix to pass out for each target
  
  for (p in 1:period) {           #for annual this is equal to 1
    #select target subset if doing seasonal
    #******need method of dealing with different targets
    
    for (j in 1:ncol(data)) {              #loop for each variable, to check its scale type.
      switch(targetType[j],
             "frac" = {temp[i.pp[[p]],j]=data[i.pp[[p]],j]*target[j]},
             "diff" = {temp[i.pp[[p]],j]=data[i.pp[[p]],j]+target[j]},
             -99.00
             )
    }
  }
  temp=data.frame(temp)
  names(temp)=varType
  return(temp)
}

#TEST
#simple.scaling(target=target,targetType=targetType, data=mat,varType=varType,period=1,index=i.pp )

#----------------------------------------------
# try to imitate format seen in generic generation functions in WGEN_lib.R