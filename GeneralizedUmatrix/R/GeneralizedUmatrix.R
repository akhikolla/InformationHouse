GeneralizedUmatrix=function(Data,ProjectedPoints,PlotIt=FALSE,Cls=NULL,Toroid=TRUE,Tiled=FALSE,ComputeInR=FALSE){
#V=GeneralizedUmatrix(ProjectedPoints,Data)
#V=GeneralizedUmatrix(ProjectedPoints,Data,TRUE,Cls,TRUE)
# Generalisierte U-Matrix fuer Projektionsverfahren
# Erklaerung der Funktion in
#D:\Subversion\lehre\Vorlesungen\KnowledgeDiscovery\01Text\04ProjektionenUndVisualisierung\81DatenbionischeProjektionen\05Schwaerme\DataBionicSwarm.pptx
#Folie 45-51
# INPUT
# Data[1:n,1:d]                          array of data: n cases in rows, d variables in columns
# ProjectedPoints[1:n,OutputDimension]   n by OutputDimension matrix containing coordinates of the Projection: A matrix of the fitted configuration.
                                         # see Projections/R/..
# OPTIONAL
# PlotIt                                 bool, defaut=FALSE, if =TRUE: U-Marix of every current Position of Databots will be shown
# toroid
## ComputeInR                                  =T: Rcode, =F Cpp Code

# Output
# Umatrix[Lines,Columns                     umatrix (see ReadUMX() dbt.DataIO)
# EsomNeurons[1:Lines,1:Columns,1:weights]       3-dimensional numeric array (wide format), not wts (long format)
# Bmu[1:n,OutputDimension]                 GridConverted Projected Points information converted by convertProjectionProjectedPoints() to predefined Grid by Lines and Columns
# gplotres                                Ausgabe von ggplot
# unbesetztePositionen                    Umatrix[unbesetztePositionen] =NA
# author: MT 06/2015
#1.Editor; MT 12/2015
#2 Editor: MT 02/2020   switched to faster plotting
toroid=Toroid
# ErrorHandling
if(any(is.nan(Data))) stop('Data contains NaNs')
if(!is.matrix(ProjectedPoints)) stop('ProjectedPoints has to be a matrix')
if(!is.matrix(Data)){
  stop('Data has to be a matrix')
}
n=nrow(ProjectedPoints)
c=ncol(ProjectedPoints)
if(c>3 |c<2)
  stop(paste0('Wrong number of Columns of ProjectedPoints: ',c))
if(c==3){
  ProjectedPoints=ProjectedPoints[,2:3]
}
#end  ErrorHandling


##########################################################################
# calcUmatrixHLP function
##########################################################################
calcUmatrixHLP <- function(EsomNeurons,Toroid=T){
  # Umatrix=calcUmatrixHLP(wts)
  # Calculate the Umatrix for given EsomNeurons projection
  # INPUT
  # EsomNeurons[Lines,Columns,weights]		neuronen aus EsomNeurons
  # OPTIONAL
  # Toroid				planar=F
  # OUTPUT
  # Umatrix[Lines,Columns]
  #############################################
  ## Nachbarn()
  nachbarn <- function(k, EsomNeurons, Toroid=FALSE){
    # INPUT
    # k Gitterpunkt
    # EsomNeurons[Lines,Columns,weights]
    # Toroid
    # OUTPUT
    # nb
    M <- dim(EsomNeurons)[1]
    N <- dim(EsomNeurons)[2]
    if(Toroid){ 
      pos1 = c(k[1]-1,k[1]-1,k[1]-1,k[1],k[1],k[1]+1,k[1]+1,k[1]+1) %% M
      pos1[which(pos1==0)] = M
      pos2 = c(k[2]-1,k[2],k[2]+1,k[2]-1,k[2]+1,k[2]-1,k[2],k[2]+1) %% N
      pos2[which(pos2==0)] = N
      nb = cbind(pos1,pos2)
    }else{# planar
      if(k[1] == 1){
        if (k[2] == 1){
          nb = rbind(c(1,2), c(2,2), c(2,1))
        }else{ 
          if (k[2] == N){
            nb = rbind(c(1,(N-1)), c(2,(N-1)), c(2,N))
          }else{ 
            nb = rbind(c(1,(k[2]-1)), c(1,(k[2]+1)), c(2,(k[2]-1)), c(2,k[2]), c(2,(k[2]+1)))
          }
        }
      }#end fall 1 fuer planar
      if(k[1] == M){
        if(k[2] == 1){
          nb = rbind(c((M-1),1), c((M-1),2), c(M,2))
        }else{
          if(k[2] == N){
            nb = rbind(c((M-1),(N-1)), c((M-1),N), c(M,(N-1)))
          }else{
            nb = rbind(c((M-1),(k[2]-1)), c((M-1),k[2]), c((M-1),(k[2]+1)), c(M,(k[2]-1)), c(M,(k[2]+1)))
          }
        }
      }#end fall 2 fuer planar
      if(k[1] != 1 && k[1] != M){
        if(k[2] == 1){
          nb = rbind(c((k[1]-1),1), c((k[1]-1),2), c(k[1],2),c((k[1]+1),1), c((k[1]+1),2))
        }else{ 
          if(k[2] == N){
            nb = rbind(c((k[1]-1),(N-1)), c((k[1]-1),N), c(k[1],(N-1)), c((k[1]+1),(N-1)), c((k[1]+1),N))
          }else{
            nb = rbind(c((k[1]-1),(k[2]-1)), c((k[1]-1),k[2]), c((k[1]-1),(k[2]+1)), c(k[1],(k[2]-1)),c(k[1],(k[2]+1)), c((k[1]+1),(k[2]-1)), c((k[1]+1),k[2]), c((k[1]+1),(k[2]+1)))
          }
        }
      }#end fall 3 fuer planar
    }# end if Toroid
    return(nb)
  }
  ############################################
  k = dim(EsomNeurons)[1]
  m = dim(EsomNeurons)[2]
  Umatrix = matrix(0,k,m)
  d=dim(EsomNeurons)[3]
  if(is.null(d)){#wts als liste
    stop('EsomNeurons wts has to be an array[1:Lines,1:Columns,1:Weights], use ListAsEsomNeurons')
  }
  for(i in 1:k){
    for(j in 1:m){
      nbs=nachbarn(c(i,j),EsomNeurons,Toroid)
      wij=EsomNeurons[i,j,]
      n.nbs=dim(nbs)[1]
      for(l in 1:n.nbs){
        nij=EsomNeurons[nbs[l,1],nbs[l,2],]
        Umatrix[i,j]=Umatrix[i,j]+sqrt(sum((wij-nij)^2))
      }	
      Umatrix[i,j]=Umatrix[i,j]/n.nbs
    }
  }
  return(Umatrix)
}
##############End calcUmatrixHLP()
#MT: Shortcut not required anymore, Now Umatrix package availible on CRAN
# ##########################################################################
# # gridNeighbourhoodPattern function
# ##########################################################################
# 
# gridNeighbourhoodPattern <- function(Radius){
#   # gridNeighbourhoodPattern(Radius)
#   # returns index for neighbours relative to (0,0) and their distances
#   #
#   # INPUT
#   # Radius		      Radius in which neighbours should be searched
#   # 
#   # OUTPUT
#   # Neighbourhood(1:m,1:3)      positions of all weights in the Neighbourhood of c(0,0) together with their distance
#   # author: Florian Lerch
#   # gridNeighbourhoodPattern(3)
#   
#   # the pattern is starting from point 0,0
#   #author: FL
#   startPoint = c(0,0) 
#   
#   # get all neighbours
#   Neighbourhood <- c()
#   
#   # make sure that Radius is an integer
#   Radius <- floor(Radius)
#   
#   # calculate only a quarter of the sphere and get the rest through symmetry
#   for (y in (startPoint[1] - Radius):(startPoint[1])){
#     for (x in (startPoint[2] - Radius):(startPoint[2])){
#       if ((y - startPoint[1])^2 + (x - startPoint[2])^2 <= Radius^2) {
#         ySym <- startPoint[1] - (y - startPoint[1])
#         xSym <- startPoint[2] - (x - startPoint[2])
#         
#         # add the new found neighbours and the points symmetric to them into the Neighbourhood
#         nn1 <- c(y, x, sqrt(y^2+x^2))
#         nn2 <- c(y, xSym, sqrt(y^2+xSym^2))
#         nn3 <- c(ySym , x, sqrt(ySym^2+x^2))
#         nn4 <- c(ySym, xSym, sqrt(ySym^2+xSym^2))
#         
#         newNeighbours <-  matrix(c(nn1, nn2, nn3, nn4),ncol=3,byrow=TRUE)
#         Neighbourhood <- rbind(Neighbourhood,newNeighbours)
#       }
#     }
#   }
#   
#   # remove duplicated entries
#   Neighbourhood <- unique(Neighbourhood)
#   
#   # sort by distance, this makes it easier to remove the farther neighbour on a Toroid
#   Neighbourhood <- Neighbourhood[order(Neighbourhood[,3]),] 
#   
#   Neighbourhood
# }
# # end GridneighbourhoodThroughPattern function
# ##########################################################################
# # neighbourhoodThroughPattern function
# ##########################################################################
# neighbourhoodThroughPattern <- function(Index, Pattern, Columns, Lines, Toroid=TRUE, RadiusBiggerThanTorus=TRUE){
#   # neighbourhoodThroughPattern(Index, Pattern, Columns, Lines, Toroid)
#   # gets a Pattern for neighbourhood and returns all indices within that pattern starting from
#   # Index
#   #
#   # INPUT
#   # Index			position of a weightvector for which the neighbours should be returned. Index in form (lines,columns)
#   # Pattern			The neighbourhoodpattern. (see gridNeighbourhoodPattern)
#   # Columns			Width of the grid
#   # Lines			Height of the grid
#   # OPTIONAL
#   # Toroid			Is the Pattern build on a Toroid?
#   # RadiusBiggerThanTorus		If the currently used radius is bigger than min(Columns/2,Lines/2)
#   #				        there needs to be an expensive check for neighbours on two
#   #				        sides over the Toroid
#   # OUTPUT
#   # Neighbourhood(1:m,1:2)       indices of all weights in the Neighbourhood of Index together with their distance
#   # author: Florian Lerch
#   # neighbourhoodThroughPattern(c(3,5), Pattern, 80, 50)
#   
#   # move points according to difference
#   Neighbourhood <- addRowWiseC(Pattern,c(Index,0))
#   
#   if(Toroid==FALSE){ # just keep points within the borders 
#     Neighbourhood <- subset(Neighbourhood, Neighbourhood[,1] >= 1 & Neighbourhood[,2] >= 1 & 
#                               Neighbourhood[,1] <= Lines & Neighbourhood[,2] <= Columns)
#     rows = Neighbourhood[,1]
#     cols = Neighbourhood[,2]
#   }
#   else if(Toroid==TRUE){ 
#     # use modulo operator to get all entries into the borders
#     rows = (Neighbourhood[,1]) %% Lines 
#     cols = (Neighbourhood[,2]) %% Columns 
#     
#     rows[(rows == 0)] = Lines
#     cols[(cols == 0)] = Columns
#   }
#   
#   indices <- (rows-1)*Columns + cols
#   
#   Neighbourhood <- cbind(indices,Neighbourhood[,3])
#   
#   #RadiusBiggerThanTorus = F
#   # on a Toroid its possible that a value exists two times in the Neighbourhood
#   if(RadiusBiggerThanTorus){
#     Neighbourhood = Neighbourhood[!duplicated(Neighbourhood[,1]),]
#   }
#   
#   Neighbourhood
# }
# ##############End neighbourhoodThroughPattern()
##########################################################################
# Koordinatenbestimmung auf Gitter
##########################################################################
if(n>4096/8)
  minNeurons=n*8
else
  minNeurons=4096
coordsres=XYcoords2LinesColumns(X=ProjectedPoints[,1],Y=ProjectedPoints[,2],PlotIt=F,minNeurons = minNeurons)

Lines=coordsres$LC[1]
Columns=coordsres$LC[2]
BMUs=coordsres$GridConvertedPoints
BMUs[,1]=Lines-BMUs[,1]+1
d=ncol(Data) #NumberOfweights
k=Lines+5
m=Columns+5

#Hier ist nur relevant, das der Radius mindestens so gross ist, das in innerhalb des Radius um ein BMU
# noch jeweils ein anderer BMU liegt
# ansonsten fuert ein hoeherer Anfangsradius nur zu einer laengeren Berechnungsdauer
HeuristischerParameter=max(round(Columns/6,0),20)
#print(HeuristischerParameter)
#HeuristischerParameter=20
##########################################################################
# Bestimmung der Positionen, auf welchen sESOM keinen Einfluss haben kann
##########################################################################

###
#only for NoN CRAN intern usage:
###
# pattern=Umatrix::gridNeighbourhoodPattern(Radius)
# bmPos=lapply(1:nrow(coordsres$GridConvertedPoints),FUN=function(i,coordsres,pattern,Radius){
#   index=coordsres$GridConvertedPoints[i,]
#   res=Umatrix::neighbourhoodThroughPattern(index,pattern,coordsres$LC[2],coordsres$LC[1],T)
#   t(sapply(res[,1],FUN=function(ind,Columns){
#     row = ((ind-1) %/% Columns) + 1
#     col = ((ind-1) %% Columns) + 1
#     c(row,col)
#   },coordsres$LC[2]))
# },coordsres,pattern,Radius)
# 
# besetztPos=do.call(rbind,bmPos)

# Radius2=max(coordsres$LC[1],coordsres$LC[2])
# pattern2=Umatrix::gridNeighbourhoodPattern(Radius2)
# res2=Umatrix::neighbourhoodThroughPattern(c(1,1),pattern2,coordsres$LC[1],coordsres$LC[2],T)
# allePositionen=t(sapply(res2[,1],FUN=function(ind,Columns){
#   row = ((ind-1) %/% Columns) + 1
#   col = ((ind-1) %% Columns) + 1
#   c(row,col)
# },coordsres$LC[2]))

# unbesetztePositionen=setdiffMatrix(allePositionen,besetztPos)$CurtedMatrix

###
#only for NoN CRAN intern usage
###
##########################################################################
# Begin. simplified ESOM Algorithmus for given BestMatching units
##########################################################################

rnd=runif(n=d*k*m, min =min(Data), max = max(Data)) #besser als min(data) bis max(data)
wts<- array(rnd,c(k,m,d)) #[Lines,Columns,weights]
# BestMatches werden festgehalten
#for(i in c(1:nrow(BMUs))){
#  wts[BMUs[i,1],BMUs[i,2],] = Data[i,]

#}
#Jeder Radius sollte min. 1 Eppoche durchlaufen werden, mehr als eine Eppoche fuehrte nicht zu mehr Emergenz
# s. auch Experimente mit iUmatrix(), wo eine Umatrix als Video pro Eppoche bei diverser Parameterwahl gezeichnet wird
epochs=HeuristischerParameter
AnfangsRadius=HeuristischerParameter
vec=pmax(seq(from=AnfangsRadius-1,by=-1,length.out = HeuristischerParameter),1)
for (i in vec){
  CurrentRadius =  i#max(AnfangsRadius-i,1) #Endradius=1
  #Algorithmus
  wts=sESOM4BMUs(BMUs,Data, wts, toroid, CurrentRadius,ComputeInR=ComputeInR)
  print(paste0('Operator: getUmatrix4BMUs() at ',round(1-(i/HeuristischerParameter),2)*100,'%'))
} # end 1:epochs
##########################################################################
# End, simplified ESOM Algorithmus for given BestMatching units
##########################################################################
  print(paste0('Operator: calcUmatrix() generates U-Matrix now...'))

  Umap=calcUmatrixHLP(wts,Toroid=toroid)

  #Umap[unbesetztePositionen]=NA
  gplotres=NULL
if(isTRUE(PlotIt)){
  if(length(Cls)!=nrow(BMUs)) Cls=rep(1,nrow(BMUs))
  gplotres=TopviewTopographicMap(GeneralizedUmatrix = Umap,BestMatchingUnits = BMUs,Cls = Cls,Tiled =Tiled,BmSize =12) #Sondern Gebirge=Unbekannte Orte der U-Matrix
  print(gplotres)
}

return(list(Umatrix=Umap,EsomNeurons=wts,Bestmatches=BMUs,Lines=Lines,Columns=Columns,
       sESOMparamaters=list(Eppochs=HeuristischerParameter,Rmax=HeuristischerParameter,Rmin=1,CoolingStrategie='Linear, Lernratate ist const =1',Toroid=toroid),
       gplotres=gplotres))


}
