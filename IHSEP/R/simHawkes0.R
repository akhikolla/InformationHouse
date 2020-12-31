simHawkes0 <-
function(nu,g,cens=1,
                      nuM=max(optimize(nu,c(0,cens),maximum=TRUE)$obj,
                        nu(0),nu(cens))*1.1,
                      gM=max(optimize(g,c(0,cens),maximum=TRUE)$obj,
                        g(0),g(cens))*1.1
                      ){
  tms <- list()
  ng <- 0; #no. of gen initialized to 1
  tmp <- simPois(nu,cens=cens,int.M=nuM)#birth times of gen ng;
  gs <- length(tmp) # size of generation ng
  while(gs>0){
    ng <- ng+1
    tms[[ng]] <- tmp
    tmp <- vector('list',gs)
    for(i in 1:gs)
      tmp[[i]] <- tms[[ng]][i]+simPois(int=g,cens=cens-tms[[ng]][i],int.M=gM)
    tmp <- unlist(tmp)
    gs <- length(tmp)
  }
  tms
}
