TileBM <- function(BestMatches, Lines, Columns){
# TiledBestMatches=TileBM(BestMatches)
# Tiles a set of bestmatches, as given by ReadBM()
#
#Input:
#BestMatches: The bestmatches, as given by ReadBM(), a n by 3 Matrix with BMKey, x and y coordinates of the n bestmatches in the corresponding columns
#Lines, Columns: The dimensions of the SOM-Matrix to which the BestMatches correspond
#
#Output:
#TiledBestMatches: n*4 by 3 Matrix, containing the BestMatches in tiled format, meaning four times the BestMatches
# author Raphael Paebst (@Florian: achtung Funktion pruefen)
# 1.Editor: MT
  
# the first n rows contain the normal BestMatches, the second n rows contain the BestMatches, with 'x' increased by the number of columns, the third n rows contain BestMatches with 'Y' increased by the number of rows and 
# the fourth n rows contain BestMatches with both 'x' and 'y' increased, thus creating four tiles with BestMatches
  
  BestMatches=as.matrix(unname(BestMatches))
if(ncol(as.matrix(BestMatches))==3){
    BestMatches=cbind(as.numeric(BestMatches[,1]),as.numeric(BestMatches[,2]),as.numeric(BestMatches[,3]))
    
    TiledBestMatches <- rbind(BestMatches, cbind(BestMatches[,1], (BestMatches[,2]+Lines), BestMatches[,3]), cbind(BestMatches[,1:2],(BestMatches[,3]+Columns)), cbind(BestMatches[,1],(BestMatches[,2]+Lines),(BestMatches[,3]+Columns))) 
# adding Columns, Lines and Columns and Lines to 'x', 'y' and 'X' and 'Y'
}else{
  if (ncol(as.matrix(BestMatches))==2){
    #warning('Number of Rows is not 3, generating Key')
    BestMatches=cbind(1:length(BestMatches[,1]),BestMatches)
    
    BestMatches=cbind(as.numeric(BestMatches[,1]),as.numeric(BestMatches[,2]),as.numeric(BestMatches[,3]))
    
    TiledBestMatches <- rbind(BestMatches, cbind(BestMatches[,1], (BestMatches[,2]+Lines), BestMatches[,3]), cbind(BestMatches[,1:2],(BestMatches[,3]+Columns)), cbind(BestMatches[,1],(BestMatches[,2]+Lines),(BestMatches[,3]+Columns))) 
  }
  else{
  stop('Error: Number of Rows is not 2 or 3, nothing could be done')
  }
}
return(TiledBestMatches = TiledBestMatches)
}
