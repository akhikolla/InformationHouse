SplineReg_fast_biv <- function(X,Y,Z,W=NULL,weights=rep(1,length(X)),InterKnotsX,InterKnotsY,n,
                               Xextr=range(X),Yextr=range(Y),flag=TRUE,center=c(sum(Xextr)/2,sum(Yextr)/2)){
  matriceX <- splineDesign(knots=sort(c(InterKnotsX,rep(Xextr,n))),derivs=rep(0,length(X)),
                            x=X,ord=n,outer.ok = T)
  matriceY <- splineDesign(knots=sort(c(InterKnotsY,rep(Yextr,n))),derivs=rep(0,length(Y)),
                            x=Y,ord=n,outer.ok = T)
  matriceY_noint <- cut_int(matriceY)
  matricebiv <- tensorProd(matriceX,matriceY_noint)
  matricebiv2 <- cbind(matricebiv,W)


  #fff <- !rankMatrix(matricebiv)==8
  #  Xknots<-makenewknots(sort(c(InterKnotsX,rep(Xextr,n-1))),deg=n)
  #  Yknots<-makenewknots(sort(c(InterKnotsY,rep(Yextr,n-1))),deg=n)
  if(all(weights==1)){
    tmp <- .lm.fit(matricebiv2, Z)
    if (tmp$rank<ncol(matricebiv2)){
      tmp <- lm.fit(matricebiv2, Z)
    }
    resid <- residuals(tmp)
  } else {
    tmp <- lm.wfit.light(matricebiv2, Z, weights) #ccc<-lm(Z ~ -1+matricebiv) #
    resid <- tmp$residuals


    if (tmp$rank<ncol(matricebiv2)){
      tmp <- lm.wfit(matricebiv2, Z, weights)
      resid <- residuals(tmp)

    }
  }
  theta <- as.numeric(coef(tmp))
#  theta[is.na(theta)] <- 0
  out <- list("Theta"=theta,"Predicted"=matricebiv2%*%theta,
              "Residuals"=resid,"RSS"=t(resid)%*%resid,
              "XBasis"= matriceX, "YBasis"= matriceY_noint,#"Poly"=poly,
              "temporary"=tmp)
  return(out)
}
