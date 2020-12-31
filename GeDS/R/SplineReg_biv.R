SplineReg_biv <- function(X,Y,Z,W=NULL,weights=rep(1,length(X)),InterKnotsX,InterKnotsY,n,
                          Xextr=range(X),Yextr=range(Y),flag=TRUE,center=c(sum(Xextr)/2,sum(Yextr)/2)){
  matriceX <- splineDesign(knots=sort(c(InterKnotsX,rep(Xextr,n))),derivs=rep(0,length(X)),
                            x=X,ord=n,outer.ok = T)
  matriceY <- splineDesign(knots=sort(c(InterKnotsY,rep(Yextr,n))),derivs=rep(0,length(Y)),
                            x=Y,ord=n,outer.ok = T)
  matriceY_noint <- cut_int(matriceY)

  matricebiv <- tensorProd(matriceX,matriceY_noint)
  matricebiv2 <- cbind(matricebiv,W)
  tmp <- lm.wfit(matricebiv2, Z, weights)
  if(any(is.na(coef(tmp)))) {
    warning("NAs in regression coefficients")
    tmp <- lm(Z~-1+matricebiv2,weights =  weights)
  }

  resid <- residuals(tmp)
  theta <- as.numeric(coef(tmp))

  #  poly<-cbind(as.matrix(expand.grid(Yknots,Xknots)),theta[1:(NCOL(matricebiv))])
  out <- list("Theta"=theta,"Predicted"=matricebiv2%*%theta,
              "Residuals"=resid,"RSS"=t(resid)%*%resid,
              "XBasis"= matriceX, "YBasis"= matriceY_noint,#"Poly"=poly,
              "Xknots" = sort(c(InterKnotsX,rep(Xextr,n))),
              "Yknots" = sort(c(InterKnotsY,rep(Yextr,n))),
              "temporary"=tmp)
  #  }
  return(out)
}


cut_int <- function(mat){
  d<-dim(mat)
  #we delete the intercept
  #otherwise tensor product basis is rank deficient
  mat_star <- t(rep(1,d[1])%*%mat)
  Q_star<-(qr.Q(qr(mat_star),complete=T))[,-1]
  #  return(mat%*%Q_star)
  return(mat)
}

