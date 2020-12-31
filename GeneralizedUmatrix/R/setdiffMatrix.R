setdiffMatrix <-
function(Matrix2Curt,Matrix2compare){
# CurtedMatrix<-setdiffMatrix(Matrix2Curt,Matrix2compare){
# setdiffMatrix shortens Matrix2Curt by those rows that are in both matrices.
# INPUT
# Matrix2Curt[n,k]            matrix, which will be shortened by x rows
# Matrix2compare[m,k]         matrix whose rows will be compared to those of Matrix2Curt
#                             x rows in Matrix2compare equal rows of Matrix2Curt (order of rows is irrelevant).
#							  Has the same number of columns as Matrix2Curt.
#
# OUTPUT List V
# V$CurtedMatrix[n-x,k]           Shortened Matrix2Curt
# V$IndCurted[x]				          Indizes, which Matrix2Curt was shortened with
# V$IndexInMatrix2comp[m]					Indizes, in Matrix2compare, which where found in  Matrix2Curt
# author: CL,MT 12/2014

 #require('fastmatch')
  #Fehlerabfang:
  
  if(is.null(Matrix2Curt)){
    return(list(CurtedMatrix=Matrix2Curt,IndCurted=NULL,IndexInMatrix2comp=NULL))
  }
  if(is.null(Matrix2compare)){
    return(list(CurtedMatrix=Matrix2Curt,IndCurted=NULL,IndexInMatrix2comp=NULL))
  }
  if(!is.matrix(Matrix2Curt)){
    return(list(CurtedMatrix=Matrix2Curt,IndCurted=NULL,IndexInMatrix2comp=NULL))
  }
  if(!is.matrix(Matrix2compare)){
    return(list(CurtedMatrix=Matrix2Curt,IndCurted=NULL,IndexInMatrix2comp=NULL))
  }
#   #RmRowIndex<-which(split(Matrix2Curt,row(Matrix2Curt)) %in% split(Matrix2compare,row(Matrix2compare)))
#   #split(matrix,row(matrix)) provides us the rows of "matrix" as list.
#   #fmatch ist schneller als match, sichtbar bei 2000 Databots: 50% Gesamtlaufzeitgewinn
#   
#   #RmRowIndex<-na.omit(fmatch(split(Matrix2compare,row(Matrix2compare)),split(Matrix2Curt,row(Matrix2Curt))))
#   a=lapply(split(Matrix2compare,row(Matrix2compare)),as.numeric)
#   b=lapply(split(Matrix2Curt,row(Matrix2Curt)),as.numeric)
#   #RmRowIndex<-as.vector(na.omit(fmatch(a,b)))  
#   #RmRowIndex<-na.omit(fmatch(a,b))
#   #DBIndex<-as.vector(na.omit(fmatch(b,a)))
#   if(nofmatch){
#     RmRowIndex=match(a,b)
#   }else{
#   RmRowIndex=fmatch(a,b)
#   }
#   IndexInMatrix2comp<-which(!is.na(RmRowIndex))  #Hoffe das ist schneller als nochmal fmatch zu machen
#   RmRowIndex=na.omit(RmRowIndex)
#   Matrix2Curt[RmRowIndex,] #Hier erst na.omit
#   #Matrix2compare[IndexInMatrix2comp,]
#   
# if(length(RmRowIndex)>0){
#     return(list(CurtedMatrix=Matrix2Curt[-RmRowIndex,],IndCurted=as.vector(RmRowIndex),IndexInMatrix2comp=IndexInMatrix2comp))
# }else{
#     return(list(CurtedMatrix=Matrix2Curt,IndCurted=NULL,IndexInMatrix2comp=IndexInMatrix2comp))
# }
# if(nofmatch){
 # library('dplyr')
 #schneller wie oben:
#requireRpackage('dplyr')
#   y=data.frame(Matrix2compare)
#   x=data.frame(Matrix2Curt)
#   FoundPositions=semi_join(x,y,by=c('row','col'))
#   OpenPositions=anti_join(x,y,by=c('row','col'))
#Komplexe zalen trick: 100mal schneller!!!

komp_gesucht     <- Matrix2Curt[ , 1] + 1i*Matrix2Curt[ , 2]
komp_durchsuchen <- Matrix2compare[ , 1] + 1i*Matrix2compare[ , 2]

Y <- komp_gesucht %in% komp_durchsuchen

FoundPositions <- Matrix2Curt[Y, ]
OpenPositions  <- Matrix2Curt[!Y, ]

  return(list(CurtedMatrix=OpenPositions,IndexInMatrix2comp=FoundPositions))
#   #return(RmRowIndex)
# }else{
# Dauert an laengesten
#     X     <- duplicated(M[  , 1:2])
#     Found <- M[X ,1:2]
#     Open  <- M[M[,3]==1 & !X,][  ,1:2]
#     return(list(CurtedMatrix=Open,IndexInMatrix2comp=Found))
# }

}
