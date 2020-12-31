# #' Fitter function for GeD Spline Regression for bivariate data
# #' @param X vector of the first independent variable.
# #' @param Y vector of the second independent variable.
# #' @param Z vector of response variable.
# #' @param W design matrix for other variables.
# #' @param weights vector of weights to be used in the fitting process.
# #' @param Indicator contingency table of \code{X} and \code{Y}.
# #' @param beta parameter for the GeDS fitting algorith. See documentation.
# #' @param alpha parameter for the stop criterium. See documentation.
# #' @param min.intknots minimum number of internal knots required.
# #' @param max.intknots maximum number of internal knots required.
# #' @param q number of iterations considered for the stop criterium. See documentation.
# #' @param Xextr boundary knots in the \code{X} direction, default to the range \code{X}
# #' @param Yextr boundary knots in the \code{Y} direction, default to the range \code{Y}
# #' @param show.iters logical indicating whether to print or not step by step informations.
# #' @param tol tolerance to be used in checking whether two knots should be considered different.
# #' @return an object of class \link{GeDS-Class}. See \link{GeDS} for further details.
# #' @description Function called from \link{GeDS},  intended to be  used   directly.
# #' @seealso \link{GeDS}
BivariateFitter <- function(X, Y, Z, W, weights=rep(1,length(X)), Indicator, beta=.5, alpha = 0.5, min.intknots = 0,
                                     max.intknots = 300, q = 2, Xextr=range(X), Yextr=range(Y), show.iters=TRUE, tol = as.double(1e-12)
) {
  save <- match.call()
  RSSnew <- numeric()
  args <- list("X" = X, "Y" = Y, "Z" = Z,"W" = W, "weights"=weights, "beta" = beta, "alpha" = alpha,"min.intknots" = min.intknots,
               "max.intknots" = max.intknots, "q" = q, "Xextr" = Xextr, "Yextr" = Yextr, "tol" = tol)
  previousX <- matrix(nrow=max.intknots+1, ncol=max.intknots+4)
  previousY <- matrix(nrow=max.intknots+1, ncol=max.intknots+4)
  nw <- NCOL(W)
  oldcoef <- matrix(nrow=max.intknots+1, ncol=(2+max.intknots/2)*(2+max.intknots/2)+nw)
  Xnodi = NULL
  Ynodi = NULL
  ordX <- order(X,Y)
  ordY <- order(Y,X)
  matrX <- matrix(ncol=3,nrow=length(X))
  nintX <- 10
  upperX <- seq(from=Xextr[1],to=Xextr[2],length= nintX+1)[-1]
  nintY <- 10
  upperY <- seq(from=Yextr[1],to=Yextr[2],length= nintY+1)[-1]
  Xctrl <- FALSE
  Yctrl <- FALSE

  for(j in 1:(max.intknots+1)){#4){#54){#
    if(j >1)  {
      if(Xctrl)  Xnodi <- sort(Xnodi)
      if(Yctrl)  Ynodi <- sort(Ynodi)
    }
    # usual regression
    first.deg<-SplineReg_fast_biv(X=X, Y=Y, Z=Z, W=W, weights=weights, InterKnotsX=Xnodi,
                                  InterKnotsY=Ynodi, Xextr=Xextr, Yextr=Yextr, n=2)#first regression #fino a 1.7 sec
    previousX[j,1:(length(Xnodi)+4)] <- sort(c(Xnodi,rep(Xextr,2)))
    previousY[j,1:(length(Ynodi)+4)] <- sort(c(Ynodi,rep(Yextr,2)))
    oldcoef[j,1:length(first.deg$Theta)] <- first.deg$Theta
    matr <- cbind(X,Y,first.deg$Residuals*weights)
    RSS.tmp <- first.deg$RSS
    RSSnew <- c(RSSnew,RSS.tmp)
      if(j>q && (j-q > min.intknots)){
        if(RSSnew[j]/RSSnew[j-q] >= alpha) {      # stop rule
          break
        }
    }



    # divide the X domain in nintX intervals
    dX <- numeric(nintX)
    upperX <- upperX+1e-15
    dX[1] <- sum(unique(X)<=upperX[1])
    for(i in 2:nintX){
      dX[i] <- sum(unique(X)<=upperX[i])-sum(dX)
    }
    zeroes <- dX==0
    matrX <- matr[ordX,]
    matrX <- makeNewMatr(matrX, Indicator, by.row=FALSE)#che va pulita qui per le osservazioni multiple
    dcumX <- cumsum(dX)


    Ymean <- numeric()
    Ywidth <- numeric()
    Ynum <- numeric()
    kk <- 1
    dXY <- numeric()

    #formazione dei cluster per ogni slice
    for(i in 1:nintX){
      if(!zeroes[i]) {
        if (i == 1){
          tmpY <- matrX[1:dcumX[i],2]
          if(length(tmpY)!=1){
            matrX[1:dcumX[i],] <- matrX[1:dcumX[i],][order(tmpY),]
          }
          tmpR <- matrX[1:dcumX[i],3]

        } else {
          tmpY <- matrX[(dcumX[i-1]+1):dcumX[i],2]
          if(length(tmpY)!=1){
            matrX[(dcumX[i-1]+1):dcumX[i],] <- matrX[(dcumX[i-1]+1):dcumX[i],][order(tmpY),]
          }
          tmpR <- matrX[(dcumX[i-1]+1):dcumX[i],3]
        }
        segni <- sign(tmpR)
        for(jj in 1:length(tmpR)){
          if(all(segni==segni[1])){
            dXY[kk]<-length(segni)
            kk = kk+1
            break
          } else {
            dXY[kk]<-min(which(segni!=segni[1])-1)
            segni <- segni[-(1:dXY[kk])]
            kk = kk+1
          }
        }
      }
    }

    #computations of residual mean and width of each cluster
    dcumXY <- cumsum(dXY)
    Ymean <- numeric(length(dXY))
    Ywidth <- numeric(length(dXY))
    Ymean[1] <- abs(mean(matrX[1:dcumXY[1],3]))
    Ywidth[1] <- diff(range(matrX[1:dcumXY[1],2]))
    for(i in 2:length(dXY)){
      Ymean[i] <- abs(mean(matrX[(dcumXY[i-1]+1):dcumXY[i],3]))
      Ywidth[i] <- diff(range(matrX[(dcumXY[i-1]+1):dcumXY[i],2]))
    }
    Ymean <- Ymean/max(Ymean)
    Ywidth <- Ywidth/max(Ywidth)

    Yweights <- beta*Ymean+(1-beta)*Ywidth
    u <- length(dXY)

    #avoids problems
    flagY <- F

    # Y knot placement
      for(i in 1:u){
        if(all(Yweights<=0)) {

          flagY <- TRUE
          break
        }
        indice <- which.max(Yweights)
        if (indice == 1) {dcumInf = 1} else {dcumInf = dcumXY[indice-1]+1}
        dcumSup <- dcumXY[indice]
        sup <- matrX[dcumSup,2]
        inf <- matrX[dcumInf,2]
        Ynewknot <- matrX[dcumSup:dcumInf,3]%*%matrX[dcumSup:dcumInf,2]/
          sum(matrX[dcumSup:dcumInf,3])
        #        tmp <- sort(c(rep(Yextr,2+1),Ynewknot,Ynodi))
        #        print(Ynewknot)
        #        print(tmp)
        #        tmp0 <- sum(tmp<=as.numeric(Ynewknot))
        #        tmp <- tmp[(tmp0-3):(tmp0+3)]
        #        tmp2 <- cbind(tmp[1:4],tmp[4:7])
        #        tmp3 <- logical(4)
        #        for(kk in 1:4){
        #          tmp3[kk]<-any(((tmp2[kk,1]<Y#matrX[dcumSup:dcumInf,2]
        #                          ) * (tmp2[kk,2]>Y#matrX[dcumSup:dcumInf,2]
        #                               )))
        #        }
        #check both conditions - no knots between Xs and Xs between knots
        if(#all(tmp3) &&
          (((dcumSup-dcumInf)!=0) && (!any((Ynodi>=inf) * (Ynodi<=sup))))  ||
          #                   (((dcumsup-dcuminf)==0) && (!any( abs(Ynodi-inf)< tol ))) )
          ((dcumSup-dcumInf)==0 && (dcumInf==1 || dcumSup==length(X)))
        ) {
          break
        } else {
          Yweights[indice] <- -Inf
        }
      }
    weightY <- Yweights[indice]


    ### qua inizia la X del nodo

    dY <- numeric(nintY)
    upperY <- upperY+1e-15
    dY[1] <- sum(unique(Y)<=upperY[1])
    for(i in 2:nintY){
      dY[i] <- sum(unique(Y)<=upperY[i])-sum(dY)
    }
    zeroes <- dY==0
    dcumY <- cumsum(dY)
    matrY <- matr[ordY,]
    matrY <- makeNewMatr(matrY, Indicator, by.row=TRUE)#che va pulita qui per le osservazioni multiple
    Xmean <- numeric()
    Xwidth <- numeric()
    Xnum <- numeric()
    kk <- 1
    dYX <- numeric()

    # andrebbe aggiunta questione indicator
    for(i in 1:nintY){
      if(!zeroes[i]) {
        if (i == 1){
          tmpX <- matrY[1:dcumY[i],1]
          if(length(tmpX)!=1){
            matrY[1:dcumY[i],] <- matrY[1:dcumY[i],][order(tmpX),]}
          tmpR <- matrY[1:dcumY[i],3]
        } else {
          tmpX <- matrY[(dcumY[i-1]+1):dcumY[i],1]
          if(length(tmpX)!=1){
            matrY[(dcumY[i-1]+1):dcumY[i],] <- matrY[(dcumY[i-1]+1):dcumY[i],][order(tmpX),]
          }
          tmpR <- matrY[(dcumY[i-1]+1):dcumY[i],3]
        }
        segni <- sign(tmpR)
        for(jj in 1:length(tmpR)){
          if(all(segni==segni[1])){
            dYX[kk]<-length(segni)
            kk = kk+1
            break
          } else {
            dYX[kk]<-min(which(segni!=segni[1])-1)
            segni <- segni[-(1:dYX[kk])]
            kk = kk+1
          }
        }
      }
    }
    dcumYX <- cumsum(dYX)
    Xmean <- numeric(length(dYX))
    Xwidth <- numeric(length(dYX))
    Xmean[1] <- abs(mean(matrY[1:dcumYX[1],3]))
    Xwidth[1] <- diff(range(matrY[1:dcumYX[1],1]))
    for(i in 2:length(dYX)){
      Xmean[i] <- abs(mean(matrY[(dcumYX[i-1]+1):dcumYX[i],3]))
      Xwidth[i] <- diff(range(matrY[(dcumYX[i-1]+1):dcumYX[i],1]))
    }
    Xmean <- Xmean/max(Xmean)
    Xwidth <- Xwidth/max(Xwidth)
    Xweights <- beta*Xmean+(1-beta)*Xwidth
    u <- length(dYX)

    flagX <- F
    for(i in 1:u){
        if(all(Xweights<0)) {
          flagX <- TRUE
          break
        }
        indice <- which.max(Xweights)
        if (indice == 1) {dcumInf = 1} else {dcumInf = dcumYX[indice-1]+1}
        dcumSup <- dcumYX[indice]
        sup <- matrY[dcumSup,1]
        inf <- matrY[dcumInf,1]
        Xnewknot<-matrY[dcumSup:dcumInf,3]%*%matrY[dcumSup:dcumInf,1]/
          sum(matrY[dcumSup:dcumInf,3])
        #    tmp <- sort(c(rep(Xextr,2+1),Xnewknot,Xnodi))
        #    tmp0 <- sum(tmp<=as.numeric(Xnewknot))
        #    tmp <- tmp[(tmp0-3):(tmp0+3)]
        #    tmp2 <- cbind(tmp[1:4],tmp[4:7])
        #    tmp3 <- logical(4)
        #    for(kk in 1:4){
        #      tmp3[kk]<-any(((tmp2[kk,1]<X#matrY[dcumSup:dcumInf,1]
        #                      ) * (tmp2[kk,2]>X#matrY[dcumSup:dcumInf,1]
        #                           )))
        #    }
        #check both conditions - no knots between Ys and Ys between knots
        if(#all(tmp3) &&
          (((dcumSup-dcumInf)!=0) && (!any((Xnodi>=inf) * (Xnodi<=sup))))  ||
          #                   (((dcumsup-dcuminf)==0) && (!any( abs(Xnodi-inf)< tol ))) )
          ((dcumSup-dcumInf)==0 && (dcumInf==1 || dcumSup==length(Y)))
        ) {break } else {
          Xweights[indice] <- -Inf
        }
      }
    weightX <- Xweights[indice]

    if(flagX && flagY) {
      print("Unable to find other knots satisfying required conditions")
      break
    } else {
      if(flagX) {
        weightX <- 0
        weightY <- 1
      } else {
        if(flagY) {
          weightY <- 0
          weightX <- 1
        }
      }
    }

    if(weightX>weightY){
      Ynewknot <- NULL
      Xctrl <- TRUE
      if(show.iters) {
        if(j>q) {toprint<-paste0("Iteration ",j,": New X Knot = ", round(Xnewknot,3)  ,
                                 ", RSS = " , round(RSSnew[j],3), ", phi = ", round(RSSnew[j]/RSSnew[j-q],3))} else {
                                   toprint<-paste0("Iteration ",j,": New X Knot = ", round(Xnewknot,3)  ,", RSS = " , round(RSSnew[j],3))}
        print(toprint)}

    } else {
      Xnewknot <- NULL
      Yctrl <- TRUE
      if(show.iters) {
        if(j>q) {toprint<-paste0("Iteration ",j,": New Y Knot = ", round(Ynewknot,3)  ,
                                 ", RSS = " , round(RSSnew[j],3), ", phi = ", round(RSSnew[j]/RSSnew[j-q],3))} else {
                                   toprint<-paste0("Iteration ",j,": New Y Knot = ", round(Ynewknot,3)  ,", RSS = " , round(RSSnew[j],3))}
        print(toprint)}

    }

    Ynodi<-c(Ynodi,Ynewknot)
    Xnodi<-c(Xnodi,Xnewknot)


    if((length(Ynodi)+3)*(length(Xnodi)+3)>=length(Z)) {
      warning("Exiting stage A: Too many knots found")
      break
      }
  }

  #STAGE B
  if (j == max.intknots + 1) {
    warning("Maximum number of iterations exceeded")
    toBeSaved <- sum(!is.na(previousX[j,]))
    previousX <- previousX[,-((toBeSaved+1):(max.intknots+4))]
    toBeSaved <- sum(!is.na(previousY[j,]))
    previousY <- previousY[,-((toBeSaved+1):(max.intknots+4))]
    lastXknots <- sum(!is.na(previousX[j,]))
    lastYknots <- sum(!is.na(previousY[j,]))
    oldcoef <- oldcoef[-((j+1):(max.intknots+1)),]
    oldcoef <- oldcoef[,-((j+2):(max.intknots+2))]
    if(j<2){
      warning("Too few internal knots found: Linear spline will not be computed. Try to set a different value for 'q' or a different treshold")
      llX <- NULL
      llY <- NULL
      lin <- NULL} else {
        ikX <- previousX[j,3:(lastXknots-2)]
        ikY <- previousY[j,3:(lastYknots-2)]
        llX <- makenewknots(ikX,2)#lin #l'approssimazione
        llY <- makenewknots(ikY,2)#lin #l'approssimazione
        lin <- SplineReg_biv(X=X, Y=Y, Z=Z, InterKnotsX=ikX,InterKnotsY=ikY, Xextr=Xextr, Yextr=Yextr, n=2)
      }
    if(j<3){
      warning("Too few internal knots found: Quadratic spline will not be computed. Try to set a different value for 'q' or a different treshold")
      qqX <- NULL
      qqY <- NULL
      squ <- NULL} else {
        qqX <- makenewknots(ikX,3)#quad
        qqY <- makenewknots(ikY,3)#quad
        squ <- SplineReg_biv(X=X, Y=Y, Z=Z, InterKnotsX=qqX,InterKnotsY=qqY, Xextr=Xextr, Yextr=Yextr, n=3)
      }
    if(j< 4){
      warning("Too few internal knots found: Cubic spline will not be computed. Try to set a different value for 'q' or a different treshold")
      ccX <- NULL
      ccY <- NULL
      cub <- NULL} else {
        ccX <- makenewknots(ikY,4)#cub
        ccY <- makenewknots(ikY,4)#cub
        cub <- SplineReg_biv(X=X, Y=Y, Z=Z, InterKnotsX=ikX,InterKnotsY=ikY, Xextr=Xextr, Yextr=Yextr, n=4)
      }

  } else {
    previousX <- previousX[-((j+1):(max.intknots+1)),] #delete also last two
    toBeSaved <- sum(!is.na(previousX[j,]))
    lastXknots <- sum(!is.na(previousX[j-q,]))
    lastYknots <- sum(!is.na(previousY[j-q,]))
    previousX <- previousX[,-((lastXknots+1):(max.intknots+4))]
    previousY <- previousY[-((j+1):(max.intknots+1)),] #delete also last two
    toBeSaved <- sum(!is.na(previousY[j,]))
    previousY <- previousY[,-((toBeSaved+1):(max.intknots+4))]
    oldcoef <- oldcoef[-((j+1):(max.intknots+1)),]
    oldcoef <- oldcoef[,-((j+2):(max.intknots+2))]
    if(j-q<2){
      warning("Too few internal knots found: Linear spline will not be computed. Try to set a different value for 'q' or a different treshold")
      llX <- NULL
      llY <- NULL
      lin <- NULL} else {
        ikX <- previousX[j-q,3:(lastXknots-2)]
        ikY <- previousY[j-q,3:(lastYknots-2)]
        llX <- makenewknots(ikX,2)#lin #l'approssimazione
        llY <- makenewknots(ikY,2)#lin #l'approssimazione
        lin <- SplineReg_biv(X=X, Y=Y, Z=Z, InterKnotsX=ikX,InterKnotsY=ikY, Xextr=Xextr, Yextr=Yextr, n=2)
      }
    if(j-q<3){
      warning("Too few internal knots found: Quadratic spline will not be computed. Try to set a different value for 'q' or a different treshold")
      qqX <- NULL
      qqY <- NULL
      squ <- NULL} else {
        qqX <- makenewknots(ikX,3)#quad
        qqY <- makenewknots(ikY,3)#quad
        squ <- SplineReg_biv(X=X, Y=Y, Z=Z, InterKnotsX=qqX,InterKnotsY=qqY, Xextr=Xextr, Yextr=Yextr, n=3)
      }
    if(j-q < 4){
      warning("Too few internal knots found: Cubic spline will not be computed. Try to set a different value for 'q' or a different treshold")
      ccX <- NULL
      ccY <- NULL
      cub <- NULL} else {
        ccX <- makenewknots(ikX,4)#cub
        ccY <- makenewknots(ikY,4)#cub
        cub <- SplineReg_biv(X=X, Y=Y, Z=Z, InterKnotsX=ccX,InterKnotsY=ccY, Xextr=Xextr, Yextr=Yextr, n=4)
      }

  }

  out <- list("Type" = "LM - biv","Linear.Knots"=list("Xk" = llX,"Yk" = llY),"Quadratic.Knots"=list("Xk" = qqX,"Yk" = qqY),
              "Cubic.Knots"=list("Xk" = ccX,"Yk" = ccY),"Dev.Linear" = lin$RSS,
              "Dev.Quadratic" = squ$RSS,"Dev.Cubic" = cub$RSS,
              "RSS" = RSSnew, "Linear" = lin, "Quadratic" = squ, "Cubic" = cub, "Stored" = previousX,
              "Args"= args,"Call"= save, "Nintknots"= list("X"= length(llX), "Y"= length(llY)),"iters" = j, "Guesses" = NULL,
              "Coefficients" = oldcoef)
  class(out) <- "GeDS"
  return(out)
}
