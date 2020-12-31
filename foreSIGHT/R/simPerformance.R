# AD: This file can be removed
#FUNCTION TO JUST MANAGE THE SIMULATION OF PERFORMANCE
simPerformance<-function(data=NULL,       #output from scenario generator
                         systemModel=NULL,
                         systemArgs=NULL,
                         simDirectory=NULL,
                         performance=NULL,
                         IOmode="suppress"
                         #add in option for multiple metrics
  
){
  
  if(IOmode=="verbose"){
    path<-file.path(".",simDirectory)
    if(!isTRUE(dir.exists(path))){
      dir.create(path)
    }
    
    path<-file.path(".",simDirectory,"systemPerformance")
    if(!isTRUE(dir.exists(path))){
      dir.create(path)
    }
  }
  
  #Running system model 
  if(is.null(performance)){
    performance=matrix(NA,nrow=data$nRep,ncol=1) # need to update to have more than one performance metric (more than 1 col)
    for(i in 1:data$nRep){
      performance[i]=systemModel(data=data$data[[i]],systemArgs=systemArgs,repID=i)
    }
  }
  
  if(length(performance)==data$nRep){
    data$performance=performance
  }else{
      cat("Length of supplied system performance vector and number of scenarios does not match.")
  }
    
  return(data)   #return data object with additional 'performance' field
    
}
# #write out performance if ok'd by user
# perfMap=cbind(data$target,performance)


