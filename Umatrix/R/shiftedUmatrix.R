shiftedUmatrix <- function(rows, cols, Umatrix, BestMatches=NULL, Cls=NULL, bufferLines=0, bufferCols=0){
# Umatrix <- shiftedUmatrix(rows, cols, Umatrix)
# shifts a matrix and positions of BestMatches by rows and cols
# shifting direction is up left (positioning by lines and rows)
# INPUT
# rows				Rows to be shifted down
# cols				Cols to be shifted right
# Umatrix(1:Lines,1:Columns)	Umatrix to be shifted
# OPTIONAL
# BestMatches(1:n,1:2)		n positions of BestMatches to be shifted
# Cls(1:n)                      Classes of the BestMatches
# bufferLines			How many lines to add above and under the Umatrix
# bufferCols                    How many Columns to add left and right of the Umatrix
# OUTPUT
# list with
# Umatrix                       Shifted Umatrix
# BestMatches			Shifted (and possibly duplicated) BestMatches
# Cls				Possibly duplicated classes
# author: Florian Lerch

  nrows <- nrow(Umatrix)
  ncols <- ncol(Umatrix)

  rows = -rows
  cols = -cols

  # change format of Cls if it is given as a vector
  if(!is.null(Cls)){
    if(is.null(dim(Cls))) Cls <- matrix(Cls)
  }

  if((bufferLines > nrows) | (bufferCols > ncols)) stop("buffer is bigger than the Umatrix itself")

  # shift the bestmatches
  if(!is.null(BestMatches)){
    # add bestmatch keys if none given
    if(ncol(BestMatches) == 2) bm_keys <- 1:nrow(BestMatches)
    else{
      bm_keys = BestMatches[,1]
      BestMatches = BestMatches[,c(2,3)]
    }

    bmrows = BestMatches[,1] + rows
    bmcols = BestMatches[,2] + cols

    bmrows = bmrows %% nrows
    bmcols = bmcols %% ncols

    bmrows[(bmrows == 0)] = nrows
    bmcols[(bmcols == 0)] = ncols

    BestMatches = cbind(bm_keys,bmrows, bmcols)
  }

  rows <- rows%%nrows
  cols <- cols%%ncols

  if(rows != 0){
    upperPart <- Umatrix[1:(nrows-rows),]
    lowerPart <- Umatrix[(nrows-rows+1):nrows,]
    Umatrix <- rbind(lowerPart, upperPart)
  }

  if(cols != 0){
    leftPart <- Umatrix[,1:(ncols-cols)]
    rightPart <- Umatrix[,(ncols-cols+1):ncols]
    Umatrix <- cbind(rightPart, leftPart)
  }

  # add buffer (/border)
  if(bufferCols>0){
    leftPart = Umatrix[,1:bufferCols]
    rightPart = Umatrix[,(ncols-bufferCols+1):ncols]
    Umatrix = cbind(rightPart,Umatrix,leftPart)
  }
  if(bufferLines > 0){
    upperPart = Umatrix[1:bufferLines,]
    lowerPart = Umatrix[(nrows-bufferLines+1):nrows,]
    Umatrix = rbind(lowerPart,Umatrix,upperPart)
  }

  # shift BestMatches according to border
  if(!is.null(BestMatches)){
    BestMatches[,2] = BestMatches[,2] + bufferLines
    BestMatches[,3] = BestMatches[,3] + bufferCols
  }

  # add BestMatches to the border
  if((!is.null(BestMatches)) & ((bufferLines+bufferCols)!=0)){
    leftBestMatchesInd = BestMatches[,3] <= (bufferCols + bufferCols)
    rightBestMatchesInd = BestMatches[,3] > ncols # equals bufferCols + (ncols - bufferCols)

    leftBestMatches <- cbind(BestMatches[leftBestMatchesInd,c(1,2)], BestMatches[leftBestMatchesInd, 3] + ncols)
    rightBestMatches = cbind(BestMatches[rightBestMatchesInd,c(1,2)], BestMatches[rightBestMatchesInd,3]-ncols)

    BestMatches = rbind(BestMatches, leftBestMatches, rightBestMatches)

    if(!is.null(Cls))
      Cls <- rbind(Cls, matrix(Cls[leftBestMatchesInd]), matrix(Cls[rightBestMatchesInd]))

    topBestMatchesInd = BestMatches[,2] <= bufferLines + bufferLines
    bottomBestMatchesInd = BestMatches[,2] > nrows

    topBestMatches <- cbind(BestMatches[topBestMatchesInd,1], BestMatches[topBestMatchesInd,2] + nrows, BestMatches[topBestMatchesInd,3])
    bottomBestMatches <- cbind(BestMatches[bottomBestMatchesInd,1], BestMatches[bottomBestMatchesInd,2] - nrows, BestMatches[bottomBestMatchesInd,3])

    BestMatches = rbind(BestMatches, topBestMatches, bottomBestMatches)

    if(!is.null(Cls))
      Cls <- rbind(Cls, matrix(Cls[topBestMatchesInd]), matrix(Cls[bottomBestMatchesInd]))
  }

  list(Umatrix=Umatrix, BestMatches=BestMatches, Cls=Cls)
}
