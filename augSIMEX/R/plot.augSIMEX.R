plot.augSIMEX<-function(x, ...){
  extrapolation<-"both"
  xlim=c(-1.2,2.2)
  p.names <- names(x$coefficients)
  variable = x$err.true
  nvar<-length(variable)
  ### if the user give the input of the who variable directly
  if (nvar>length(p.names)) {
    variable<-substitute(variable)
    nvar<-1}
  ### if the user give a vector of the name of variable or the index of variable
  if (mode(variable)=="character") {
    if (nvar==1){
    index<-which(p.names==variable)
    } else index<-unlist(lapply(1:nvar,FUN = function(i){which(p.names[i]==variable)}))
  } else index<-variable 
  
  for (ii in index){
      betas<-x$coefmatrix[,ii]
      lambda<-x$lambda
      plot(lambda,betas,
           xlab="lambda", ylab="Estimated Coefficient",xlim=xlim,...)
      if (extrapolation=="linear"|extrapolation=="both"){
        model.linear <- lm(betas~lambda)
        extrapoint.linear<-data.frame(lambda=-1)
        extravalue.linear<-predict(model.linear,extrapoint.linear)
        points(-1,extravalue.linear,col=2,pch = 16)
        points.linear<-data.frame(lambda=seq(-1.5, 2.2,0.01))
        values.linear<-predict(model.linear,points.linear)
        lines(points.linear$lambda,values.linear,col=2,lty=1,lwd=2)
        if (extrapolation=="linear") {
          legend("topright",legend=c("linear"),col=c(2),lty=c(1),lwd=2)
        }
      }
      if (extrapolation=="quadratic"|extrapolation=="both"){
        lambda2<-lambda^2
        model.quadratic <- lm(betas~lambda+lambda2)
        extrapoint.quadratic<-data.frame(lambda=-1,lambda2=1)
        extravalue.quadratic<-predict(model.quadratic,extrapoint.quadratic)
        points(-1,extravalue.quadratic,col=4,pch = 16)
        points.quadratic<-data.frame(lambda=seq(-1.5, 2.2,0.01),
                                     lambda2=seq(-1.5, 2.2,0.01)^2)
        values.quadratic<-predict(model.quadratic,points.quadratic)
        lines(points.quadratic$lambda,values.quadratic,col=4,lty=8,lwd=2)
      }
      if (extrapolation=="quadratic") {
        legend("topright",legend=c("quadratic"),col=c(4),lty=c(8),lwd=2)
      } else legend("topright",legend=c("linear","quadratic"),col=c(2,4),lty=c(1,8),lwd=2)
  }
}



