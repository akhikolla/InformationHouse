jackknife <- function(data,pp,fit){
  out <- matrix(NA,ncol=ncol(data),nrow=nrow(data))
  data <- as.matrix(data)
  beta <- pp$ipar$thres[2,]
  if(fit=="lz"){
    for(i in 1:ncol(data)){
  capture.output(pp.lz <- PP_4pl(respm=data[,-i],thres=pp$ipar$thres[2,][-i], lowerA=pp$ipar$lowerA[-i], slopes=pp$ipar$slopes[-i], upperA=pp$ipar$upperA[-i],type=pp$type, mu=pp$ipar$mu[-i], sigma2=sqrt(pp$ipar$sigma2[-i])))
      out[,i] <- lz(data=data[,-i], thetas=pp.lz$resPP$resPP[,"estimate"], betas=pp.lz$ipar$thres[2,], lowerAs=pp.lz$ipar$lowerA, slopes=pp.lz$ipar$slopes, higherAs=pp.lz$ipar$upperA)[,fit]
    }
  }
  if(fit=="lzstar"){
    for(i in 1:ncol(data)){
      capture.output(pp.ls <- PP_4pl(respm=data[,-i],thres=pp$ipar$thres[2,][-i], lowerA=pp$ipar$lowerA[-i], slopes=pp$ipar$slopes[-i], upperA=pp$ipar$upperA[-i],type=pp$type, mu=pp$ipar$mu[-i], sigma2=sqrt(pp$ipar$sigma2[-i])))
      out[,i] <- lzstar(data=data[,-i], thetas=pp.ls$resPP$resPP[,"estimate"], betas=pp.ls$ipar$thres[2,], lowerAs=pp.ls$ipar$lowerA, slopes=pp.ls$ipar$slopes, higherAs=pp.ls$ipar$upperA,method=pp.ls$type, mu=pp.ls$ipar$mu, sigma=sqrt(pp.ls$ipar$sigma2))
    }
  }
    
  out.se <- sqrt( ((ncol(data)-1)/ncol(data)) * rowSums((out - rowMeans(out, na.rm=TRUE))^2, na.rm=TRUE) )
  return(out.se)
}
