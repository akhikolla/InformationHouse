XYcoords2LinesColumns=function(X,Y,minNeurons=4096,MaxDifferentPoints=F,PlotIt=T){
  #XYcoords2LinesColumns(X,Y)
  # Converts points given as x(i),y(i) coordinates to integer coordinates Columns(i),Lines(i)
  # INPUT
  # X(1:n), Y(1:n)                    coordinates: x(i),y(i) is the i-th point on a plane
  #
  # OPTIONAL
  # minNeurons                       minimal size of the corresponding grid i.e max(Lines)*max(Columns)>=MinGridSize , default MinGridSize = 4096
  #                                  defined by the numer of neurons
  # MaxDifferentPoints                TRUE: the discretization error is minimal
  #                                  FALSE: number of Lines and Columns is minimal
  # PlotIt                           Plots the result
  # OUTPUT
  # GridConvertedPoints[1:Columns,1:Lines,2]  IntegerPositions on a grid corresponding to x,y
  # LC vector of:
  # Lines                                   GridLines==max(GridConvertedPoints[,2])
  # Columns                                 GridColumns==max(GridConvertedPoints[,1])

  # author MT 08/2015
  multiplot <-
    function(...,
             plotlist = NULL,
             file,
             cols = 1,
             layout = NULL) {
      requireNamespace('grid')
      requireNamespace('ggplot2')
      # Make a list from the ... arguments and plotlist
      plots <- c(list(...), plotlist)
      
      numPlots = length(plots)
      
      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                         ncol = cols,
                         nrow = ceiling(numPlots / cols))
      }
      
      if (numPlots == 1) {
        print(plots[[1]])
        
      } else {
        # Set up the page
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
          # Get the i,j matrix positions of the regions that contain this subplot
          matchidx <-
            as.data.frame(which(layout == i, arr.ind = TRUE))
          
          print(plots[[i]],
                vp = grid::viewport(
                  layout.pos.row = matchidx$row,
                  layout.pos.col = matchidx$col
                ))
        }
      }
    }
  
  temp = cbind(X, Y)
  minX <- min(X)
  minY <- min(Y)
  maxX <- max(X)
  maxY <- max(Y)
  RangeXfirst = maxX - minX
  RangeYfirst = maxY - minY
  
  frac = RangeYfirst / RangeXfirst
  
  Lines = -frac - 1 + sqrt(frac ^ 2 + 1 + (minNeurons - 2) * frac)
  Columns = RangeXfirst / RangeYfirst * (Lines - 1) + 1
  
  i = 0
  while (round(Lines, 0) * round(Columns, 0) < minNeurons) {
    i = i + 1
    Lines = -frac - 1 + sqrt(frac ^ 2 + 1 + (minNeurons - 2) * frac) + i
    Columns = RangeXfirst / RangeYfirst * (Lines - 1) + 1
  }
  
  # Update min and  max.
  minX <- min(X)
  minY <- min(Y)
  maxX <- max(X)
  maxY <- max(Y)
  RangeX = maxX - minX
  RangeY = maxY - minY
  
  X <- floor((X - minX) / RangeX * (Columns - 1) + 1)
  Y <- floor((Y - minY) / RangeY * (Lines - 1) + 1)
  
  maxX <- max(X)
  maxY <- max(Y)
  if (maxX > Columns)
    warning('Columnsrange exceeded')
  if (maxY > Lines)
    warning('Linesrange exceeded')
  
  rCoord = cbind(Y, X)
  
  requireNamespace('mgcv')
  nvorher = nrow(mgcv::uniquecombs(temp))
  nnachher = nrow(mgcv::uniquecombs(rCoord))
  
  #nvorher = nrow(uniquePoints(temp)$unique)
  #nnachher = nrow(uniquePoints(rCoord)$unique)
  
  if (MaxDifferentPoints) {
    i = 0
    while (0.9 * nvorher > nnachher) {
      i = i + 1
      Lines = -frac - 1 + sqrt(frac ^ 2 + 1 + (minNeurons - 2) * frac) + i
      Columns = RangeXfirst / RangeYfirst * (Lines - 1) + 1
      X = temp[, 1]
      Y = temp[, 2]
      minX <- min(X)
      minY <- min(Y)
      maxX <- max(X)
      maxY <- max(Y)
      RangeX = maxX - minX
      RangeY = maxY - minY
      
      X <- floor((X - minX) / RangeX * (Columns - 1) + 1)
      Y <- floor((Y - minY) / RangeY * (Lines - 1) + 1)
      maxX <- max(X)
      maxY <- max(Y)
      if (maxX > Columns)
        warning('Columnsrange exceeded')
      if (maxY > Lines)
        warning('Linesrange exceeded')
      
      rCoord = cbind(Y, X)
      nvorher = nrow(mgcv::uniquecombs(temp))
      nnachher = nrow(mgcv::uniquecombs(rCoord))
      #nvorher = nrow(uniquePoints(temp)$unique)
      #nnachher = nrow(uniquePoints(rCoord)$unique)
      if (nvorher != nnachher) {
        print(paste(
          'Different Points are now on the same grid position:',
          nvorher,
          nnachher
        ))
      }
    } # end while
  }# end if MaxDifferentPoints
  if (PlotIt) {
    #  requireRpackage('ggplot2')
      requireNamespace('ggplot2')
    sp2 <-
      ggplot2::ggplot(as.data.frame(rCoord), ggplot2::aes(X, Y)) + ggplot2::geom_point() + ggplot2::ggtitle("Projected Points")
    sp1 <-
      ggplot2::ggplot(as.data.frame(temp), ggplot2::aes(X, Y)) + ggplot2::geom_point() + ggplot2::ggtitle("Grid Converted Points")
    multiplot(sp1, sp2, cols = 2)
    sp2 + ggplot2::coord_fixed(ratio = Columns / Lines)
    sp1 + ggplot2::coord_fixed(ratio = RangeXfirst / RangeYfirst)
    
    #   windows()
    #   par(mfrow = c(1,2))
    #   plot(sp1)
    #   plot(sp2)
    # PlotProjectedPoints(temp,NULL,FALSE,main='ProjectedPoints')
    #PlotProjectedPoints(rCoord,NULL,T,main='Corresponding BestMatches')
  }
  return(list(
    GridConvertedPoints = rCoord,
    LC = c(Lines = round(Lines, 0), Columns = round(Columns, 0))
  ))
}