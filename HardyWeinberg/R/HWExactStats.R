HWExactStats <- function(X,x.linked=FALSE,plinkcode=TRUE,midp=FALSE,...) {
  X <- as.matrix(X)
  nmarkers <- nrow(X)
  pvalvec <- numeric(nmarkers)
  if(x.linked) {
    if(plinkcode) {
      for (i in 1:nmarkers) {
        pvalvec[i] <- SNPHWEX(X[i,4],X[i,3],X[i,5],X[i,1],X[i,2],midp)  
      }
    } else {
      for (i in 1:nmarkers) {
        pvalvec[i] <- HWExact(X[i,],x.linked=x.linked,...)$pval
      }  
    }
  } else { # not x.linked
    if(plinkcode) {
      for (i in 1:nmarkers) {
        pvalvec[i] <- SNPHWE2(X[i,2],X[i,1],X[i,3],midp)
      }
    } else {
      for (i in 1:nmarkers) {
        pvalvec[i] <- HWExact(X[i,],x.linked=x.linked,...)$pval
      }  
    }
  }
  return(pvalvec)
}
