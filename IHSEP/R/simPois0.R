simPois0 <-
function(int,cens=1,int.M=optimize(int,
                                  c(0,cens),maximum=TRUE)$obj*1.1
                     ){
  N <- rpois(1,int.M*cens)
  if(N==0)return(numeric(0))
  tms <- runif(N,0,cens)
  tms <- tms[as.logical(mapply(rbinom,n=1,size=1,
                               prob=int(tms)/int.M) )]
  sort(tms)
}
