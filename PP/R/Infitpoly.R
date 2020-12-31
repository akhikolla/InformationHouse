Infitpoly <- function( data,
                       thetas,
                       thresholds, 
                       slopes=NULL
){
  # ------------------------------------------------------------------------------------------------
  betas <- as.vector(apply(thresholds[-1,],2,cumsum))
  betas <- betas[!is.na(betas)]
  if(is.null(slopes)) slope <- rep(1,length(betas))
  # ------------------------------------------------------------------------------------------------
  #  information
  X <- data
  N.mat <- apply(X, 1, function(x){ length(na.exclude(x)) })
  k <- apply(X,2,max,na.rm=TRUE)
  # processed Items 
  Xproc <- 1 * !is.na(X)  
  # ------------------------------------------------------------------------------------------------
  # parameter information
  ai <- slope
  # ------------------------------------------------------------------------------------------------
  # category information
  k_seq <- sequence(k)
  k_seq0 <- integer(0)
  for (i in k) k_seq0 <- c(k_seq0, 0:i)
  k_seq1 <- integer(0)
  for (i in k) k_seq1 <- c(k_seq1, c(1,1:i))
  # ------------------------------------------------------------------------------------------------
  # Itemspecific information
  k_item <- rep(1:length(k),times = k)
  k_item0 <- rep(1:length(k),times = k+1)
  theta_mat <- matrix( thetas, nrow = nrow( X ), ncol = sum( k ) ) 
  # ------------------------------------------------------------------------------------------------
  theta_mat_kat <- tcrossprod(theta_mat , diag(k_seq))
  # ---------------------------------------------
  submatrix_pijx <- (t(theta_mat_kat) - betas) #* ai 
  
  cat_pijx_0 <- tapply(1:length(betas),k_item,function(x) {
    submat <- rbind(rep(0,times = ncol(submatrix_pijx)),submatrix_pijx[x,])
      tmp <- apply(submat,2,function(xx){ 
        kat <- exp(xx)
        nom <- sum(kat,na.rm = TRUE)
        dev <- kat / nom
        return(dev)
      })
    return(tmp)   
  })
  
  cat_pijx <- lapply(cat_pijx_0,function(x){
    x <- x[-1,]
  })
  
  Pijx <- t(do.call(rbind,cat_pijx))
  Pijx <- Xproc[,k_item] * Pijx
  
  Emat.l <- tcrossprod(Pijx , t(diag(k_seq)))
  Emat <- tapply(1:length(betas),k_item,function(x) rowSums(Emat.l[,x],na.rm=TRUE))
  
  Eij <- t(do.call(rbind,Emat))
  Eij <- Xproc * Eij
  
  Pijx.0 <- t(do.call(rbind,cat_pijx_0))
  
  Pijx.0 <- Xproc[,k_item0] * Pijx.0
  
  Eij.0 <- t( apply(Eij[,k_item0],1,function(x) {k_seq0 - x}) )
  Eij.0 <- Xproc[,k_item0] * Eij.0
  # Qni <- Xproc * (Pni * ( 1-Pni ))
  
  # Variance
  Vmat.cat <- (Eij.0^2)*Pijx.0
  Vmat.l <- tapply(1:length(k_seq0),k_item0,function(x) rowSums(Vmat.cat[,x],na.rm=TRUE))
  Wni.mat <- Xproc * t(do.call(rbind,Vmat.l))
  # Wni.mat <- t(do.call(rbind,Vmat.l))
  
  Cmat.cat <- (Eij.0)^4*Pijx.0
  Cmat.l <- tapply(1:length(k_seq0),k_item0, function(x) {rowSums(Cmat.cat[,x],na.rm=TRUE)})
  Cni.mat <- Xproc * t(do.call(rbind,Cmat.l))
  # Cni.mat <- t(do.call(rbind,Cmat.l))
  
  Yni.mat <-  X - Eij
  #Yni.mat <-  X - Eij
  Zni.mat <- Yni.mat / sqrt( Wni.mat )
  Zni2.mat <- Zni.mat^2

  
  # INFIT MEANSQ
  Vn <- rowSums( (Yni.mat^2), na.rm=TRUE ) / rowSums( Wni.mat, na.rm=TRUE ) 
  
  # standardized INFIT
  # Variance term
  qni2 <- rowSums(Cni.mat - Wni.mat^2,na.rm=TRUE ) / (rowSums(Wni.mat,na.rm=TRUE))^2
  # infit
  ti <- ( (Vn^(1/3)) - 1)*(3/sqrt(qni2))+(sqrt(qni2)/3)
  
  # additional output
  chisq   <- rowSums( (Zni2.mat), na.rm=TRUE )
  pvalue  <- 1 - pchisq(rowSums( (Zni2.mat), na.rm=TRUE ), N.mat-1 )
  df      <- N.mat - 1 
  
  out <- cbind(
    "infit"    = round(Vn,6),
    "in_t"     = round(ti,6),
    "in_chisq" = round(chisq,3),
    "in_df"    = df,
    "in_pv"    = round(pvalue,3)
  )
  return(out)
}