outcomes.3 <- function(n.a) {
  nAA <- 0
  Res <- NULL
  nAlim <- limit(n.a[1])
  for(i in 0:nAlim) {
    nAA <- i
    for(j in 0:(n.a[1]-2*i)) {
      nAB <- n.a[1] - 2*i - j
      nAC <- j
      nb.left <- n.a[2] - nAB
      nc.left <- n.a[3] - nAC
      if(nb.left <= nc.left) {
        for(l in 0:limit(nb.left)) {
          nBB <- l
          nBC <- nb.left - 2*l
          nCC <- (nc.left - nb.left)/2 + l
          r <- c(nAA,nAB,nAC,nBB,nBC,nCC)
          names(r) <- c("AA","AB","AC","BB","BC","CC")
          Res <- rbind(Res,r)
        }
      } else {
        for(l in 0:limit(nc.left)) {
          nCC <- l
          nBC <- nc.left - 2*l
          nBB <- (nb.left - nc.left)/2 + l
          r <- c(nAA,nAB,nAC,nBB,nBC,nCC)
          names(r) <- c("AA","AB","AC","BB","BC","CC")
          Res <- rbind(Res,r)
        }
      }
    }
  }
  #  rownames(Res) <- 1:nrow(Res)
  return(Res)
}