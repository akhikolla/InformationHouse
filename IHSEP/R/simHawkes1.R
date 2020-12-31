simHawkes1 <-
  function(nu=NULL,g=NULL,cens=1,
           nuM=max(optimize(nu,c(0,cens),maximum=TRUE)$obj,
             nu(0),nu(cens))*1.1,
           gM=max(optimize(g,c(0,cens),maximum=TRUE)$obj,
             g(0),g(cens))*1.1,
           exp.g=FALSE,
           gp=c(1,2)
           ){
    tms <- list()
    ng <- 0; #no. of gen initialized to 1
    tmp <- simPois(nu,cens=cens,int.M=nuM)#birth times of gen ng;
    gs <- length(tmp) # size of generation ng
    pois.shift <- function(shift){
      shift+simPois(int=g,cens=cens-shift,int.M=gM,exp.int=exp.g,par=gp)
    }

    while(gs>0){
      ng <- ng+1
      tms[[ng]] <- tmp
      tmp <- unlist(sapply(tmp,pois.shift,simplify=FALSE))
      gs <- length(tmp)
    }
    tms
  }
