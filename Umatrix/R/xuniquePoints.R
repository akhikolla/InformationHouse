uniquePoints <- function(data,mindist = 1e-10){
  # U <- uniquePoints(data)
  # return only the unique points in data
  #
  # INPUT
  # data[1:n,1:d]   The vector/matrix with n data points of dimension d
  #				          the points are in the  rows
	# mindist					everything with less distance is considered equal
  #
  # OUTPUT
  # a list U containg:
  # UniqueData  <- U$unique[1:u,1:d]       the data points  without duplicate points
  # UniqueInd   <- U$sortind[1:u]		       an index vector such that unique ==  data[sortind,]
  # Uniq2DataInd<- U$mergeind[1:n] 	       an index vector such that   data ==  unique[mergeind,]
  # IsDuplicate <- U$IsDuplicate[1:n,1:n]  for i!=j IsDuplicate[i,j]== 1  iff data[i,] == data[j,] ;  IsDuplicate[i,i]==0
  #
  # Complete Code rewrite by FP 03/17 
  # 1.Editor: MT 03/17: mindist als Parameter eingefuehrt
  
  
  # when data is a vector, convert to matrix
  if (class(data) == "numeric" || class(data) == "complex") {
    data <- matrix(data, ncol = 1)
  } else if (class(data) == "data.frame") {
    data <- as.matrix(data)
  } else if (class(data) != "matrix") {
    stop("uniquePoints input is neither a (numeric or complex) vector, matrix or data.frame.")
  }
  
  NumData = nrow(data)
  
  distsmat = as.matrix(dist(data))
  # * 1 is to "hack" TRUE to 1 and FALSE to 0
  IsDuplicate = ((distsmat) < mindist) * 1 - diag(NumData)
  
  # Hack: set distance of every node to itself to NA (we're not interested in these distances)
  # otherwise every element would be a duplicate of itself later.
  distsmat = distsmat + diag(rep(NA, nrow(distsmat)))
  
  # No duplicates found
  if (length(which((distsmat) < mindist)) == 0) {
    return(list(
      "unique" = unique(data),
      "sortind" = c(1:NumData),
      "mergeind" = c(1:NumData),
      IsDuplicate = IsDuplicate
    ))
  }
  
  # save rownames
  origrows = rownames(data)
  # use rownames to save original index
  rownames(data) = c(1:NumData)
  
  #  find out which distances are smaller than mindist (and as such duplicates)
  # (except the diagonal, so we artificially increase it)
  dups = which((distsmat) < mindist, arr.ind = T)
  
  # indicies of rows which will be deleted
  delinds = c()
  while (dim(dups)[1] > 0) {
    ind = dups[1, 1]
    delinds = c(delinds, ind)
    # remove every row in which the duplicate also occurs.
    # this is to prevent to also delete the "original"
    # if not save, encapuslate right hand side in matrix( ... , nrow = 2)
    dups = dups[-(which(dups == ind, arr.ind = T)[, 1]), ]
  }
  # delete duplicates, stay as matrix (important for vector input) and keep rownames
  uniquedata = data[-delinds, ]
  uniquedata = matrix(uniquedata, ncol = ncol(data))
  rownames(uniquedata) = rownames(data)[-delinds]
  
  # the rownames contain the indicies in the original data
  sortind = as.numeric(rownames(uniquedata))
  
  # calculate the mergind. At the end there should not be a NA in mergeind left
  mergeind = rep(NA, NumData)
  mergeind[sortind] = 1:nrow(uniquedata)
  # this works due to the fact that:
  #   union(sortind, delinds) == 1:NumData (if everything is sorted)
  for (i in delinds) {
    candidates = which(IsDuplicate[i, ] == 1)
    original = intersect(candidates, sortind)
    mergeind[i] = which(sortind == original)
  }
  
  # restore rownames of the original data
  rownames(uniquedata) = origrows[-delinds]
  
  return(
    list(
      "unique" = as.matrix(uniquedata),
      "sortind" = sortind,
      "mergeind" = mergeind,
      IsDuplicate = IsDuplicate
    )
  )
  
}
