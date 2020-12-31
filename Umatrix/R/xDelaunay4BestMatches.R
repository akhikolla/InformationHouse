Delaunay4BestMatches <- function(BestMatches, MatrixOrSize, IsToroid = TRUE,PlotIt=FALSE,XYcoords=F){
# Delaunay4BestMatches(BestMatches, MatrixOrSize, IsToroid,PlotIt)
# Calculates the adjacency matrix of the delaunay graph for bestmatches in tiled form if BMs are located on a toroid grid
#
# INPUT
# BestMatches[1:n,1:3]            n by 3 matrix containing the BMKey, X and Y coordinates of the n BestMatches
#                                 BestMatches NEED NOT BE UNIQUE!
#                                 however, there is an edge in the Deaunay between duplicate points!  
#
# OPTIONAL
# MatrixOrSize[2]                 Defalut 50,82; A vector of length 2, containing the number of lines and columns of the SOM corresponding to the BestMatches
# IsToroid                        logical, indicating if BM's are on a toroid grid. Default is True
# PlotIt                          Set PlotIt=TRUE, if you want to see the Plots
# XYcoords                        default: FALSE, if TRUE: XYKoordinaten werden nicht zu Lines,Columns umgerechnet
# OUTPUT
# a list of:
# Delaunay[1:n,1:n]               adjacency matrix of the Delaunay-Graph
# TiledDelaunay                   if IsToroid,  this is the Tiled Delaunay Graph 
# X, Y                            x- and y-coordinates of the tiled bestmatches (or just bestmatches if IsToroid = FALSE).
  
# authors:  ALU / Michael Thrun    Aug.2014
# uses   DelaunayGraphMatrix()

  if(is.vector(MatrixOrSize)&length(MatrixOrSize)==2){ # MatrixOrSize=c(..,..) 
    Lines  =MatrixOrSize[1]
    Columns=MatrixOrSize[2]
  }else{ # MatrixOrSize=0 Matrix
    Columns=ncol(MatrixOrSize)
    Lines  =nrow(MatrixOrSize)
  }# end(is.vector(MatrixOrSize)&length(MatrixOrSize)==2)
  # jetzt gibts Lines und Columns

  if(IsToroid){      # behandlung eines pack-man (toroiden) MapSpaces
      TiledBestMatches=TileBM(BestMatches, Lines = Lines, Columns = Columns);

      if(XYcoords){
        TiledX=TiledBestMatches[,2]
        TiledY=TiledBestMatches[,3]
      }else{
        TiledX =TiledBestMatches[,3]; 
        TiledY = Lines*2+1-TiledBestMatches[,2]; # So rechnet man BM koordinaten in XY um
      }

      D <- DelaunayGraphMatrix(TiledX,TiledY,PlotIt=PlotIt);
      TiledDelaunay =  D$Delaunay;
      Key = BestMatches[,1];  
      AnzBestMatches = length(Key);
      Offset = c(0,1,2,3)*AnzBestMatches
      DelaunayAdjazenzMatrix = TiledDelaunay[c(1:AnzBestMatches),c(1:AnzBestMatches)]; # initialisieren mit dem oberen viertel
      for (i in Offset){
        for(j in Offset){
          DelaunayAdjazenzMatrix = DelaunayAdjazenzMatrix + TiledDelaunay[c(1:AnzBestMatches)+i,c(1:AnzBestMatches)+j];
        }# end for j
      }# end for i 
      # zurueckrechnen auf die ungekachelten werte
      #DelaunayAdjazenzMatrix <-UntileGraph4BM(TiledDelaunay, TiledBestMatches, Lines, Columns);
      DelaunayAdjazenzMatrix =(DelaunayAdjazenzMatrix>0)*1; # alles auf 0/1 reduzieren
      X=TiledX #MT: Berechnung fehl!
      Y=TiledY
  
  }
  else{# MapSpace ist planar
    if(XYcoords){
      X = BestMatches[,2];          #So rechnet man BM koordinaten in XY um
      Y = BestMatches[,3];  #So rechnet man BM koordinaten in XY um
    }else{
      X = BestMatches[,3];          #So rechnet man BM koordinaten in XY um
      Y = Lines+1-BestMatches[,2];  #So rechnet man BM koordinaten in XY um
    }
      D <- DelaunayGraphMatrix(X,Y,PlotIt);
      DelaunayAdjazenzMatrix  = D$Delaunay;
      TiledDelaunay = DelaunayAdjazenzMatrix;
  }# end if(IsToroid)

  # Ausgabe zusammenstellen
  return(list(Delaunay = DelaunayAdjazenzMatrix, TiledDelaunay = TiledDelaunay, X =X, Y =Y))
}# end function Delaunay4BestMatches
