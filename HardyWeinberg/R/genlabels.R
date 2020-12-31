genlabels <- function(X) {
    if(is.vector(X)) {
        if(length(X)==3) names(X) <- c("AA","AB","BB")
        if(length(X)==5) names(X) <- c("A","B","AA","AB","BB")
    }
    if(is.matrix(X) | is.data.frame(X)) {
        if(ncol(X)==3) colnames(X) <- c("AA","AB","BB")
        if(ncol(X)==5) colnames(X) <- c("A","B","AA","AB","BB")
    }
  return(X)
}
