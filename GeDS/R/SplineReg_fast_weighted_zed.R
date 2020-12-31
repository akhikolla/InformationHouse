SplineReg_fast_weighted_zed <- function(X,Y,Z,offset,weights=rep(1,length(X)),InterKnots,n,extr=range(X)){
  matrice <- splineDesign(knots=sort(c(InterKnots,rep(extr,n))),derivs=rep(0,length(X)),x=X,ord=n,outer.ok = T)
  matrice2 <- cbind(matrice,Z)

  Y0 <- Y-offset
  tmp <-  if(all(weights==1)) .lm.fit(matrice2, Y0) else lm.wfit.light(matrice2, Y0,weights)
  theta <- coef(tmp)
  predicted <- matrice2%*%theta+offset
  resid <- Y - predicted
  out <- list("Theta"=theta,"Predicted"=predicted,
              "Residuals"=resid,"RSS"=t(resid)%*%resid,
              "Basis"= matrice2,
              "temporary"=tmp)
  return(out)
}
