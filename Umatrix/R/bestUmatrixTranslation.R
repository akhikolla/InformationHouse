bestUmatrixTranslation <- function(Umatrix, BestMatches=NULL){
# rotation = bestUmatrixTranslation(Umatrix, BestMatches)
# Get the ideal amount of lines and columns to rotate the Umatrix by, so that the borders
# are on positions with the greatest Heights and fewest BestMatches
# INPUT
# Umatrix(1:Lines, 1:Columns)	  Umatrix
# OPTIONAL
# BestMatches(1:n,1:2)		  BestMatches
# OUTPUT
# list with
# lines				  number of lines to rotate up
# cols				  number of columns to rotate left
# author: Florian Lerch, Michael Thrun

  ncols = ncol(Umatrix)
  nrows = nrow(Umatrix)

  # calculate Height for every line and column
  lineHeights = rowSums(Umatrix)
  colHeights = colSums(Umatrix) 

  if(!is.null(BestMatches)){
    # add bm keys if not given
    if(ncol(BestMatches) == 2) BestMatches <- cbind(1:nrow(BestMatches),BestMatches)

    # find intruding BestMatches on every line and column
    lineIntruders = rep(1,nrows)
    colIntruders = rep(1,ncols)
    for(i in 1:nrows) lineIntruders[i] = length(which(BestMatches[,2] == i))
    for(i in 1:ncols) colIntruders[i] = length(which(BestMatches[,3] == i))
    
    lineIntruders2 = rep(1,nrows)
    colIntruders2 = rep(1,ncols)
    for(i in 2:(nrows-1)) lineIntruders2[i] = lineIntruders[(i-1)]+lineIntruders[i]+lineIntruders[(i+1)]

    for(i in 2:(ncols-1)) colIntruders2[i] = colIntruders[(i-1)]+colIntruders[i]+colIntruders[(i+1)]

    lineIntruders2[1]=lineIntruders[nrows]+lineIntruders[1]+lineIntruders[2]
    lineIntruders2[nrows]=lineIntruders[nrows]+lineIntruders[1]+lineIntruders[nrows-1]
    
    colIntruders2[1]=colIntruders[ncols]+colIntruders[1]+colIntruders[2]
    colIntruders2[ncols]=colIntruders[ncols]+colIntruders[1]+colIntruders[ncols-1]
    
    colIntruders2[colIntruders2==0]=1
    lineIntruders2[lineIntruders2==0]=1
    
    # norm Heights with intruding BestMatches
    lineHeights = lineHeights / lineIntruders2
    colHeights = colHeights / colIntruders2
  }

  # calculate the upper left cutting point
  cutLine = which.max(lineHeights)
  cutCol = which.max(colHeights)

  list(lines = cutLine, cols = cutCol)
}
