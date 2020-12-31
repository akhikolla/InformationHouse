lzstar <- function(
                    data, 
                    thetas, 
                    betas, 
                    lowerAs=NULL, 
                    slopes=NULL,
                    higherAs=NULL,
                    method=c("mle","wle","map"), 
                    mu=0, 
                    sigma=1 
                    ){ 

  if(is.null(slopes)){slopes <- rep(1,length(betas))}
  if(is.null(lowerAs)){lowerAs <- rep(0,length(betas))}
  if(is.null(higherAs)){higherAs <- rep(1,length(betas))}
  if(!all(apply(data,2,function(x) { all(na.omit(x) %in% 0:1) }))) stop("Please check the input, only 0/1/NA are allowed \n")

  ai <- slopes
  ci <- lowerAs
  di <- higherAs
  # processed Items 
  Xproc <- 1 * !is.na(data)  
  Nproc <- rowSums(Xproc) #instead on ncol(data)
  
  # calculate the propability of each pearson to NOT solve an item
  submatrix <-   t(outer(thetas,betas,'-'))*ai

  # probabilit.
  Pi <- t(ci + ( (( di-ci )*exp(submatrix)) / ( 1+exp(submatrix) ) ))

  Qi <- 1-Pi

  # first derivation of Pi, it is needed for `r`
  d1Pi <- t(((di-ci) * ai * exp(submatrix)) / (1 + exp(submatrix))^2)
  
  # different calculation of ri, in dependend of the choosen method
  # CHECK....
  ri <- d1Pi / (Pi * Qi)
    if( all( method!=c("mle","wle","map") ) ){
      stop("Pleas choose either 'mle', 'wle' or 'map' \n")
    }
      if(method=="mle"){
        # "MLE"
        r0 <- 0
      }
      if(method=="wle"){
        # "WLE"
        d2Pi <- t( ((di-ci) * (ai^2) * exp(submatrix) * (1-exp(submatrix))) / (1 + exp(submatrix))^3 )
        r0 <- rowSums( ri * d2Pi ) / (2 * rowSums( ri * d1Pi  ))       
      }
      if(method=="map"){
        # dependend on the choosen distribution
        r0 <-  (mu - thetas) / sigma^2
      }
  wi            <- log( Pi / Qi )
  Wn            <- rowSums( Xproc*(data - Pi) * wi ,na.rm=TRUE)
  sigmaNtheta2  <- rowSums(wi^2 * Pi  * Qi) / Nproc
  cn            <- rowSums(d1Pi * wi) / rowSums(d1Pi * ri)
  widach        <- wi - cn * ri
  tau2          <- rowSums(widach^2 * Pi * Qi) / Nproc
  
  lzstern       <- (Wn + cn*r0) / sqrt(Nproc * tau2)
  lzstern <- round(lzstern,6)
  out <- matrix(lzstern,ncol=1)
  colnames(out) <- "lzstar"
  return(out)
}
