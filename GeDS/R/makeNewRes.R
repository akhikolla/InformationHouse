makeNewRes <- function(resold, recurr){
  ids <- cumsum(recurr)
  ids<- c(0,ids)
  newres <- numeric((length(ids)-1))
  for(i in 1:(length(ids)-1)){
    newres[i] <- sum(resold[(ids[i]+1):ids[i+1]])/recurr[i]
  }
  return(newres)
}


makeNewRes2 <- function(resold, recurr, weights){
  ids <- cumsum(recurr)
  ids<- c(0,ids)
  newres <- numeric((length(ids)-1))
  newweights <- numeric((length(ids)-1))

  for(i in 1:(length(ids)-1)){
    newres[i] <- sum(resold[(ids[i]+1):ids[i+1]])/recurr[i]
    newweights[i] <- sum(weights[(ids[i]+1):ids[i+1]])
  }
  newres <- newres*newweights
  return(newres)
}
