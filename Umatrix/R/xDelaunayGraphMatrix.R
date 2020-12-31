DelaunayGraphMatrix <- function(X,Y,PlotIt=FALSE){
# D <- DelaunayGraphMatrix(X,Y,PlotIt);
# Calculates the Adjacency Martix for a Delaunay-Graph
# INPUT
# X[1:n],  Y[1:n]         point coordinates for which a planar Delaunay graph is constucted
#                         POINTS NEED NOT BE UNIQUE!
#                         however, there is an edge in the Delaunay between duplicate points!
# OPTIONAL
# PlotIt                  if TRUE Delaunay graph is plotted
#
# OUTPUT
# a list of
# Delaunay[1:n,1:n]       The adjacency matrix of the Delaunay-Graph.
#
# TRI                     NOT JET IMPLEMENTED:  for the Moment TRI == delaunayn(X,Y)
requireNamespace('geometry')
# uses packages deldir, geometry

 # if(!require(deldir)){ # ggf package deldir laden
    # Calculates the Delaunay triangulation and the
    # Voronoi tessellation (with respect to the entire plane) of a planar point set.
 #   install.packages('deldir')
 #   library(deldir)
 # }else{
#    library(deldir)
#  }# end if(!require(deldir))

#  if(!require(geometry)){ #to calculate the corners of the delaunay triangles
 #   install.packages('geometry')
#    library(geometry)
#  }else{
#    library(geometry)
 # }# end if(!require(geometry))

  # Punkte unique machen
  unique = uniquePoints(cbind(X,Y))
  UniqXY = unique$unique
  UniqueInd   = unique$sortind
  Uniq2DataInd = unique$mergeind
  IsDuplicate = unique$IsDuplicate
  UniqX =  UniqXY[,1];
  UniqY =  UniqXY[,2];

  # Delaunay ausrechnen mit deldir
  DeldirOutput = deldir(UniqX ,UniqY ); #
  PointsAndIndices = DeldirOutput$delsgs # dadrin stecken die indices des Delaunays von -> Nach
  FromInd =  PointsAndIndices$ind1;    #  indices der Ausgangspunkte des Delanays
  ToInd   =  PointsAndIndices$ind2     #  indices der Endpunkte des Delanays

  # Adjazenzmatrix befuellen
  UniqDelaunay = zeros(length(UniqX),length(UniqY)); # Adjazenzmatrix initialisieren
  for (i in c(1:length(FromInd))) {
    UniqDelaunay[FromInd[i],  ToInd[i]]  <- 1 ;#Only Direct neighbours A and B get an one from A to B
    UniqDelaunay[ToInd[i]  ,FromInd[i]]  <- 1 ;#Only Direct neighbours A and B get an one from A to B
  } # end for i neighbours

#gplot(UniqDelaunay,UniqXY); # Testplot

 # jetzt uniqe points wieder auf originale uebertragen
  Delaunay = zeros(length(X),length(Y));
  Delaunay = UniqDelaunay[Uniq2DataInd,Uniq2DataInd]
  Delaunay = Delaunay + IsDuplicate # noch je eine Verbindung zwischen den Doubletten eintragen


  if (PlotIt) { # Plot der gesamten Adjazenzmatrix mit doppelten punkte
    gplot(Delaunay,cbind(X,Y)); # Plot der gesamten Adjazenzmatrix mit doppelten punkten
  }# end if

 # default fuer TRI
  TRI = geometry::delaunayn(cbind(X,Y));

 # ausgabe zusammenstellen
  return(list(Delaunay = Delaunay, TRI = TRI))
}# end function  DelaunayGraphMatrix

  #  ALTER CODE !
#   if(!require(tripack)){ # checking for required library tripack and installing or loading it
#      install.packages("tripack")
#      library(tripack)
#   }else{
#     library(tripack)
#   }
#   if(!require(geometry)){ # checking for required library geometry and installing or loading it
#      install.packages("geometry")
#      library(geometry)
#   }else{
#     library(geometry)
#   }
#
#
# numberofbms <- length(X); #getting the number of points in X (and therefore Y as well).
# uniquified <- uniquePoints(cbind(X,Y)); # Removing all duplicate points.
# numberofuniquepoints <- nrow(uniquified$unique); # getting the number of unique points, needed for the removal of simulated toroid points
# uniqueX <- uniquified$unique[,1];
# uniqueY <- uniquified$unique[,2];
# uniqueindex <- uniquified$sortind;
# recoveryindex <- uniquified$mergeind; # the vector of indices to recover removed duplicates.
# Delaunay <- matrix(0, nrow = numberofuniquepoints, ncol = numberofuniquepoints); # initializing the adjacency matrix of dimensions length(x) by
# if (Tiling) { #Recommended case with a simulated toroid universe
# toroidify <- CornerPoints(uniqueX, uniqueY); #Adding points for the toroid universe
# toroidX <- toroidify$XPlusCorner;
# toroidY <- toroidify$YPlusCorner;
# #del <- tri.mesh(toroidX,toroidY, "remove") #calculating Delaunay-Tesselation
# voronoi <- voronoi.mosaic(toroidX, toroidY, "remove") # calculating the voronoi tesselation
# deladjs <- cells(voronoi)
# for (i in 1:numberofuniquepoints){
# neighbours <- which(deladjs[[i]]$neighbours <= numberofuniquepoints)
# Delaunay[i, neighbours] <- 1; #going through deladjs by row. all indices in a row are neighbours and are set to 1
# }
# TRI <- delaunayn(cbind(toroidX, toroidY))
# TRI <- TRI[- which(TRI[,] > numberofuniquepoints, arr.ind = TRUE),]; # removing triangles containing simulated points from the list
# }
# else { # calculating the graph with no simulated toroid universe, not recommended
# #del <- tri.mesh(uniqueX, uniqueY, "remove"); # like above, calculating the Delaunay graph.
# voronoi <- voronoi.mosaic(uniqueX, uniqueY, "remove")
# deladjs <- cells(voronoi) # like above, extracting the indices of adjacent vertices
# for (i in 1:length(deladjs)) { # like above, filling the adjacency matrix
# Delaunay[i, deladjs[[i]]$neighbours]  <- 1; #going through deladjs by row. The two indices in a row are adjacent points, so the corresponding value in the adjacency matrix is set to 1
# }
# TRI <- delaunayn(cbind(uniqueX, uniqueY))
# }
# Delaunay <- Delaunay[recoveryindex, recoveryindex]; # completing the adjacency matrix for the duplicate points.
#
# return(list(Delaunay = Delaunay, TRI = TRI))
# }
