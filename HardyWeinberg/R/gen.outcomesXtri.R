gen.outcomesXtri <- function(na,nb,nc,nm,nf) {
  #
  # generates all outcomes of tri-allelelic exact test for the X chromosome
  #
  original.allele.string <- LETTERS[1:3]
  Res <- NULL
  ma.max <- min(nm,na)
  for(ma in ma.max:0) {
    mbleft <- min(nm - ma,nb)
    for(mb in mbleft:0) {
      mc <- min(nm-ma-mb,nc)
      aleft <- na-ma
      bleft <- nb-mb
      cleft <- nc-mc
      al.left <- c(aleft,bleft,cleft) # left for females
      #      print(al.left)
      ind <- order(al.left)
      al <- original.allele.string[ind]
      genotype.string <- c(paste(c(al[1],al[1]),collapse=""),
                           paste(sort(c(al[1],al[2])),collapse=""),
                           paste(sort(c(al[1],al[3])),collapse=""),
                           paste(c(al[2],al[2]),collapse=""),
                           paste(sort(c(al[2],al[3])),collapse=""),
                           paste(c(al[3],al[3]),collapse="")) # order of genotypes according to the original minor allele coding
      acounts <- al.left[ind] # sorted
      Out <- outcomes.3(acounts)
      colnames(Out) <- genotype.string
      Out <- Out[,c( "AA","AB","AC","BB","BC","CC")] # re-order female genotypes to correspond to original allele coding
      if(is.vector(Out)) {
        Out <- matrix(Out,nrow=1)
      }
      nout <- nrow(Out)
      Block <- cbind(rep(ma,nout),rep(mb,nout),rep(mc,nout),Out)
      Res <- rbind(Res,Block)
    }
  }
  colnames(Res) <- c("A","B","C","AA","AB","AC","BB","BC","CC")
  rownames(Res) <- 1:nrow(Res)
  return(Res)
}
