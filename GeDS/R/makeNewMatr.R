makeNewMatr <- function(matrice, tab, by.row=F){
  recurr <-if(by.row) c(t(tab)) else c(tab)
  recurr <- recurr[recurr!=0]
  ids <- cumsum(recurr)
  ids<- c(0,ids)
  newres <- numeric((length(ids)-1))
  newX <- numeric((length(ids)-1))
  newY <- numeric((length(ids)-1))
  for(i in 1:(length(ids)-1)){
    newres[i] <- sum(matrice[(ids[i]+1):ids[i+1],3])
    newX[i] <- matrice[i,1]
    newY[i] <- matrice[i,2]
  }
  return(cbind(newX,newY,newres))
}
