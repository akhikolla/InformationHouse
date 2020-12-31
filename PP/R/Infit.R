Infit <- function( data, 
                   thetas, 
                   betas, 
                   lowerAs=NULL, 
                   slopes=NULL, 
                   higherAs=NULL,... ){
  if(is.null(slopes)) slopes <- rep(1,length(betas))
  if(is.null(lowerAs))lowerAs <- rep(0,length(betas))
  if(is.null(higherAs))higherAs <- rep(1,length(betas))
  
  if(!all(apply(data,2,function(x) { all(na.omit(x) %in% 0:1) }))) stop("Please check the input, only 0/1/NA are allowed \n")
  X <- data
  # processed Items 
  Xproc <- 1 * !is.na(X)  
  N <- as.numeric(apply(X,1,function(x) sum(!is.na(x))))
  
  ai <- slopes
  ci <- lowerAs
  di <- higherAs
  # calculate the propability of each pearson to not solve an item
  submatrix <- (matrix( thetas, ncol = nrow( data ), nrow = ncol( data ) ,byrow = TRUE) - betas) * ai
  
  # --------------------------------------------------------------------------------
  # probabilit. of solving an Item Expected Response (in der Literatur Pi nik)
  Pni     <- t( ci + (di-ci)/(1+exp(-submatrix)) )
  # --------------------------------------------------------------------------------
  
  # second: calculate for each the Prop
  # ------------------------------------------------  
  Qni <- Xproc * (Pni * ( 1-Pni ))
  
  # Variance of xni
  Wni <- Xproc * ( (((0-Pni)^2) * (1-Pni) ) + (((1-Pni)^2) * Pni) )
  
  # Kurtosis of xni
  Cni <- Xproc * ( (((0-Pni)^4) * (1-Pni)) + (((1-Pni)^4) * Pni) )
  
  # Score Residual: (Pni in literature Eni)
  Yni <- Xproc * (X - Pni)
  
  # Standardized Residual (Qni in literature Wni)
  Zni <- Yni / sqrt( Qni )
  
  # Fit Mean Square
  
  # INFIT  (weighted)
  Vn <- rowSums( (Yni)^2, na.rm=TRUE ) / rowSums( Qni, na.rm=TRUE )
  
  # standardized INFIT
  # Variance term
  qni2 <- rowSums(Cni - Wni^2,na.rm=TRUE ) / (rowSums(Wni,na.rm=TRUE))^2
  # infit.z
  ti <- ( (Vn^(1/3)) - 1)*(3/sqrt(qni2))+(sqrt(qni2)/3)
  

  
  # additional output
  chisq   <- rowSums( ((Zni)^2), na.rm=TRUE )
  Zni2 <- rowSums( ((Zni)^2), na.rm=TRUE )
  pvalue  <- 1 - pchisq(Zni2, N-1 )
  df      <- N - 1 
  
  out <- cbind(
    "infit"    = round(Vn,6),
    "in_t"     = round(ti,6),
    "in_chisq" = round(chisq,3),
    "in_df"    = df,
    "in_pv"    = round(pvalue,3)

  )
  return(out)
  
}