gplot <- function(AdjacencyMatrix,Coordinates,Xlabel=' ',Ylabel=' ',xlim=c(min(VertexX),max(VertexX)),ylim=c(min(VertexY),max(VertexY)),LineSpec= 'black'){
# gplot(AdjacencyMatrix,Coordinates,Xlabel,YlabelLineSpec)
# as in MATLAB:
# plots a graph as in graph theory for an AdjacencyMatrix
# given the (2D or 3D )coordinates of the points
#
# INPUT
# AdjacencyMatrix[1:n, 1:n]  where d is2 or 3 Adjazenzmatrix des Graphen
# Coordinates[1:n, 1:d]      where d is2 or 3 coordinates of the vertices
# Xlabel,Ylabel              labels for x and y axis
# LineSpec                   line colour see plot(...col =LineSpec)
#
# Author ALU 2014

#   AdjacencyMatrix = Delaunay
#   Coordinates = cbind(x,y)
#   Xlabel =' '
#   Ylabel =' '

  
  Columns=ncol(AdjacencyMatrix)
  Lines  =nrow(AdjacencyMatrix)
  n = min(Columns,Lines)
  d = ncol(Coordinates)
  VertexX = Coordinates[,1];
  VertexY = Coordinates[,2];
  
  xlab = Xlabel; ylab =Ylabel;
  
  if (d==2){# planar graph
    plot(VertexX,VertexY,  type = 'p', xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab); # die knoten
    for(i in c(1:n)){ # die matrix durchgehen und plotten
      for(j in c(1:n)){ # die matrix durchgehen und plotten
        if (AdjacencyMatrix[i,j]>0){ # kante i-> j zeichnen 
          FromToX = c( VertexX[i],VertexX[j]);
          FromtoY = c( VertexY[i],VertexY[j]);
          par(new = TRUE)
          plot(FromToX,FromtoY,col = LineSpec,axes = FALSE, type = 'l', xlim = xlim, ylim = ylim, xlab =xlab, ylab = ylab);
        }# end if verbindung zeichnen
      }# end for j
    }# end for i
    
  }else{ # 3d graph
    print('gplot(): 3d graph plotting not jet implemented')
  }# end if (d==2)
  
} # end function gplot()

