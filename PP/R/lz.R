# LZ-Statistik
lz <- function( data, 
                thetas, 
                betas, 
                lowerAs=NULL, 
                slopes=NULL,
                higherAs=NULL,...){ 
  if(is.null(slopes)) slopes <- rep(1,length(betas))
  if(is.null(lowerAs))lowerAs <- rep(0,length(betas))
  if(is.null(higherAs))higherAs <- rep(1,length(betas))
  if(!all(apply(data,2,function(x) { all(na.omit(x) %in% 0:1) }))) stop("Please check the input, only 0/1/NA are allowed \n")

  
     ai <- slopes
     ci <- lowerAs
     di <- higherAs
     
     # processed Items 
     Xproc <- 1 * !is.na(data)  
     Nproc <- rowSums(Xproc)
     
  # calculate the propability of each person to solve an item ("-x" within the probabilit.)
  submatrix <- (matrix( thetas, ncol = nrow( data ), nrow = ncol( data ) ,byrow = TRUE) - betas) * ai

   Pi     <- t( ci + (di-ci)/(1+exp(-submatrix)) )
   Pi_1   <- 1-Pi

  l0      <- rowSums( ( (data * log(Pi)) + ((1-data) * log(Pi_1)) )*Xproc/Nproc,na.rm=TRUE)
  mean_l0 <- rowSums( Xproc*((Pi * log(Pi)) + ((Pi_1) * log(Pi_1))) ,na.rm=TRUE)
  var_l0  <- rowSums( Xproc*(Pi * (Pi_1)) * (log(Pi/Pi_1)^2) ,na.rm=TRUE)
  lz      <- (Nproc*l0 - mean_l0) / sqrt(var_l0) 
  lz <- round(lz,6)
  l0 <- round(l0,6)
  out <- cbind("lz"=lz,"lz_unst"=l0)
  return(out)
  
}
