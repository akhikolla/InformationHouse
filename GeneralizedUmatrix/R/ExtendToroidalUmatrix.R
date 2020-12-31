ExtendToroidalUmatrix=function(Umatrix,Bestmatches,ExtendBorders){
  bmu=Bestmatches
  
  if(is.finite(ExtendBorders)){
  if(ExtendBorders>0){
    Umatrix<-cbind(Umatrix[,c((ncol(Umatrix)-ExtendBorders-1):ncol(Umatrix))],Umatrix,Umatrix[,1:ExtendBorders])
    Umatrix<-rbind(Umatrix[c((nrow(Umatrix)-ExtendBorders-1):nrow(Umatrix)),],Umatrix,Umatrix[1:ExtendBorders,])

    bmu[,1]=bmu[,1]+length(c((nrow(Umatrix)-ExtendBorders-1):nrow(Umatrix)))
    bmu[,2]=bmu[,2]+length(c((ncol(Umatrix)-ExtendBorders-1):ncol(Umatrix)))
  }
  
  if(ExtendBorders<0){
    warning("ExtendToroidalUmatrix assumes an extend matrix which can be shortened with an negative 'ExtendBorders' values.")
          Umatrix<-Umatrix[-c((nrow(Umatrix)-ExtendBorders+1):nrow(Umatrix)),]
          Umatrix<-Umatrix[-c(1:(ExtendBorders+1)),]
          Umatrix<-Umatrix[,-c((ncol(Umatrix)-ExtendBorders+1):ncol(Umatrix))]
          Umatrix<-Umatrix[,-c(1:(ExtendBorders+1))]
          
          bmu[,1]=bmu[,1]-length(c((nrow(Umatrix)-ExtendBorders+1):nrow(Umatrix)))
          bmu[,2]=bmu[,2]-length(c(ncol(Umatrix)-ExtendBorders+1):ncol(Umatrix))
  }
    
  }
  
  return(list(Umatrix=Umatrix,Bestmatches=bmu))
}