CGGP_internal_imputesomegrid <- function(CGGP,y,theta,...,ystart = NULL) {
  Q  = max(CGGP$uo[1:CGGP$uoCOUNT,]) # Max value of all blocks
  
  if(!is.matrix(y)){
    numoutputs = 1
    yimputed <- y
  }else{
    numoutputs = dim(y)[2]
    yimputed <- y
  }
  for(oplcv in 1:numoutputs){
    if(!is.matrix(y)){
      y.thisloop = y
    }else{
      y.thisloop = as.vector(y[,oplcv])
    }
    
    if(any(is.na(y.thisloop))){
      Is = sort(which(is.na(y.thisloop)))
      
      if(!is.matrix(y)){
        thetaMAP.thisloop = theta
        
        cholS = list(matrix(1,1,1),Q*CGGP$d) # To store choleskys
        # Loop over each dimension
        for (dimlcv in 1:CGGP$d) {
          # Loop over each possible needed correlation matrix
          for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
            Xbrn = CGGP$xb[1:CGGP$sizest[levellcv]]
            Xbrn = Xbrn[order(Xbrn)]
            Sstuff = CGGP$CorrMat(Xbrn, Xbrn ,  thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],return_dCdtheta = FALSE)
            S = Sstuff
            # When theta is large (> about 5), the matrix is essentially all 1's, can't be inverted
            solvetry <- try({
              cS = chol(S)
              cholS[[(dimlcv-1)*Q+levellcv]]= cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
            }, silent = TRUE)
          }
        }
        cholS.thisloop =cholS
      }else{
        if(is.matrix(theta)){
          thetaMAP.thisloop = theta[,oplcv]
          cholS = list(matrix(1,1,1),Q*CGGP$d) # To store choleskys
          # Loop over each dimension
          for (dimlcv in 1:CGGP$d) {
            # Loop over each possible needed correlation matrix
            for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
              Xbrn = CGGP$xb[1:CGGP$sizest[levellcv]]
              Xbrn = Xbrn[order(Xbrn)]
              Sstuff = CGGP$CorrMat(Xbrn, Xbrn ,  thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],return_dCdtheta = FALSE)
              S = Sstuff
              # When theta is large (> about 5), the matrix is essentially all 1's, can't be inverted
              solvetry <- try({
                cS = chol(S)
                cholS[[(dimlcv-1)*Q+levellcv]]= cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
              }, silent = TRUE)
            }
          }
          cholS.thisloop =cholS
        }else{
          thetaMAP.thisloop = theta
          
          if(oplcv <1.5){
            cholS = list(matrix(1,1,1),Q*CGGP$d) # To store choleskys
            # Loop over each dimension
            for (dimlcv in 1:CGGP$d) {
              # Loop over each possible needed correlation matrix
              for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
                Xbrn = CGGP$xb[1:CGGP$sizest[levellcv]]
                Xbrn = Xbrn[order(Xbrn)]
                Sstuff = CGGP$CorrMat(Xbrn, Xbrn ,  thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],return_dCdtheta = FALSE)
                S = Sstuff
                # When theta is large (> about 5), the matrix is essentially all 1's, can't be inverted
                solvetry <- try({
                  cS = chol(S)
                  cholS[[(dimlcv-1)*Q+levellcv]]= cS+t(cS)-diag(diag(cS)) #store the symmetric version for C code
                }, silent = TRUE)
              }
            }
            cholS.thisloop =cholS
          }
        }
      }
      
      xp = as.matrix(CGGP$design[Is,])
      if(dim(xp)[2]<CGGP$d){
        xp = t(xp)
      }
      n2pred = dim(xp)[1]
      
      warmstart = 1
      if(is.null(ystart)){
        warmstart = 0
      }else{
        if(!is.matrix(y)){
          if(is.matrix(ystart)){
            warmstart = 0
          }else{
            if(length(ystart)!=length(y)){
              warmstart =0
            }else{
              warmstart = 1
              yn0 = ystart
              yhat0 = ystart[Is]
              w = rep(1,n2pred)
            }
          }
        }else{
          if(!is.matrix(ystart)){
            warmstart = 0
          }else{
            if(dim(ystart)[1] !=dim(y)[1] || dim(ystart)[2] !=dim(y)[2] ){
              warmstart =0
            }else{
              warmstart = 1
              yn0 = ystart[,oplcv]
              yhat0 = ystart[Is,oplcv]
              w = rep(1,n2pred)
            }
          }
        }
      }
      
      if(warmstart<0.5){
      brokenblocks = unique(CGGP$blockassign[Is])
      possblocks = 1:max(CGGP$uoCOUNT,20)
      
      
      for (lcv in 1:n2pred){
        Bs = CGGP$blockassign[Is[lcv]]
        possblocks = unique(c(possblocks,CGGP$uala[Bs,1:sum(CGGP$uo[Bs,]>1.5)]))
      }
      possblocks = sort(possblocks)
      keepthisone = rep(1,length(possblocks))
      
      for (lcv in 1:length(possblocks)){
        if (any(Is %in% CGGP$dit[possblocks[lcv], 1:CGGP$gridsizet[possblocks[lcv]]])){
          keepthisone[lcv] = 0;
        }
      }
      possblocks = possblocks[which(keepthisone>0.5)]
      
      # Cp is sigma(x_0) in paper, correlation vector between design points and xp
      Cp = matrix(0,n2pred,CGGP$ss)
      GGGG = list(matrix(1,n2pred,length(CGGP$xb)),CGGP$d)
      for (dimlcv in 1:CGGP$d) { # Loop over dimensions
        V = CGGP$CorrMat(xp[,dimlcv], CGGP$xb[1:CGGP$sizest[max(CGGP$uo[,dimlcv])]],
                         thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
                         returnlogs=TRUE)
        GGGG[[dimlcv]] = exp(V)
        Cp = Cp+V[,CGGP$designindex[,dimlcv]]
      }
      Cp = exp(Cp)
      
      ME_t = matrix(1,n2pred,1)
      MSE_v = list(matrix(0,n2pred,2),(CGGP$d+1)*(CGGP$maxlevel+1)) 
      Q  = max(CGGP$uo[1:CGGP$uoCOUNT,])
      for (dimlcv in 1:CGGP$d) {
        for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
          gg = (dimlcv-1)*Q
          INDSN = 1:CGGP$sizest[levellcv]
          INDSN = INDSN[sort(CGGP$xb[1:CGGP$sizest[levellcv]],
                             index.return = TRUE)$ix]
          MSE_v[[(dimlcv)*CGGP$maxlevel+levellcv]] =
            CGGP_internal_postvarmatcalc_fromGMat(GGGG[[dimlcv]],
                                                  c(),
                                                  as.matrix(
                                                    cholS.thisloop[[gg+levellcv]]
                                                  ),
                                                  c(),
                                                  INDSN,
                                                  CGGP$numpara,
                                                  returndiag=TRUE)
        }
      }
      
      ME_t = matrix(1,n2pred,length(possblocks))
      for (blocklcv in 1:length(possblocks)) {
        ME_s = matrix(1,n2pred,1)
        for (dimlcv in 1:CGGP$d) {
          levelnow = CGGP$uo[possblocks[blocklcv],dimlcv]
          ME_s = ME_s*as.matrix(MSE_v[[(dimlcv)*CGGP$maxlevel+levelnow]])
        }
        ME_t[,blocklcv] = as.vector(ME_s)
      }
      
      Jstar = possblocks[apply(ME_t, 1, which.max)]
      w = rep(0,n2pred)
      w = 1-apply(ME_t, 1, max)
      
      #find x0
      yn0 = y.thisloop
      yhat0 = rep(0,n2pred)
      gg = (1:CGGP$d-1)*Q
      for(lcv in 1:length(Jstar)){
        Q  = max(CGGP$uo[1:CGGP$uoCOUNT,]) # Max value of all blocks
        IS = CGGP$dit[Jstar[lcv], 1:CGGP$gridsizet[Jstar[lcv]]];
        B = y.thisloop[IS]
        rcpp_kronDBS(unlist(cholS.thisloop[gg+CGGP$uo[Jstar[lcv],]]), B, CGGP$gridsizest[Jstar[lcv],])
        yhat0[lcv] = Cp[lcv,IS]%*%B
        yn0[Is[lcv]] = yhat0[lcv]
      }
      }
      gg = (1:CGGP$d-1)*Q
      #find g at x0
      pwforg = rep(0,length(y.thisloop))
      for (blocklcv in 1:CGGP$uoCOUNT) {
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]]
        if(any(Is %in% IS) ){
          B = yn0[IS]
          rcpp_kronDBS(unlist(cholS.thisloop[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
          pwforg[IS] = pwforg[IS]+CGGP$w[blocklcv] * B
        }
      }
      pw0 = pwforg
      pwforg = pwforg[Is]
      dir0 = pwforg
      
      wgst1 = rep(0,length(y.thisloop))
      wgst2 = rep(0,length(y.thisloop))
      wgst1[Is] = w*dir0
      dotheseblocks = rep(0,CGGP$uoCOUNT)
      for (blocklcv in 1:CGGP$uoCOUNT) {
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]]
        if(any(Is %in% IS) ){
          B = wgst1[IS]
          rcpp_kronDBS(unlist(cholS.thisloop[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
          wgst2[IS] = wgst2[IS]+CGGP$w[blocklcv] * B
        }
      }
      
      lambdas = pmax(pmin(sum(wgst2*yn0)/sum(wgst1[Is]*wgst2[Is]),1.1),-0.1)
      yhat1 = yhat0-lambdas*w*dir0;
      yn1 = yn0
      yn1[Is] = yhat1
      
      pwforg = rep(0,length(y.thisloop))
      for (blocklcv in 1:CGGP$uoCOUNT) {
        IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]]
        if(any(Is %in% IS) ){
          B = yn1[IS]
          rcpp_kronDBS(unlist(cholS.thisloop[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
          pwforg[IS] = pwforg[IS]+CGGP$w[blocklcv] * B
        }
      }
      pw1 = pwforg
      pwforg = pwforg[Is]
      dir1 = pwforg
      
      M = 20
      s = matrix(0,length(dir1),M)
      ny = matrix(0,length(dir1),M)
      dirsave = matrix(0,length(dir1),M)
      xsave = matrix(0,length(dir1),M)
      rho = rep(0,M)
      gamma = rep(0,M)
      
      #x0 = yhat0, x1 = yhat1
      #g(x0)= dir0, x1 = dir1
      lcv = 1
      s[,lcv] = yhat1-yhat0
      ny[,lcv] = dir1-dir0
      xsave[,lcv] = yhat1
      dirsave[,lcv] = dir1
      
      L = rep(0,400)
      for(lcv in 1:400){
        
        q = dir1
        if(lcv > 1.5){
          for(k in 2:min(lcv,M)){#q(min(lcv,M),1,by=-1)
            rho[k] = 1/sum(s[,k]*ny[,k])
            gamma[k] = rho[k]*sum(s[,k]*q)
            q = q-gamma[k]*ny[,k]
          }
        }
        
        r = q*mean(s[,1]*ny[,1])/mean(ny[,1]*ny[,1])
        if(lcv > 1.5){
          for(k in min(lcv,M):2){#q(min(lcv,M),1,by=-1)
            beta = rho[k]*sum(ny[,k]*r)
            r =r+(gamma[k]-beta)*s[,k]
          }
        }
        dir1u = -r
        
        wgst1 = rep(0,length(y.thisloop))
        wgst2 = rep(0,length(y.thisloop))
        wgst1[Is] = dir1u
        dotheseblocks = rep(0,CGGP$uoCOUNT)
        for (blocklcv in 1:CGGP$uoCOUNT) {
          IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]]
          if(any(Is %in% IS) ){
            B = wgst1[IS]
            rcpp_kronDBS(unlist(cholS.thisloop[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
            wgst2[IS] = wgst2[IS]+CGGP$w[blocklcv] * B
          }
        }
        
        lambdas =-sum(wgst2*yn1)/sum(wgst1[Is]*wgst2[Is])
        
        for(k in  (min(lcv,M-1)):1){
          s[,k+1] = s[,k]
          ny[,k+1] = ny[,k]
          xsave[,k+1] = xsave[,k]
          dirsave[,k+1] = dirsave[,k]
        }
        
        yhat2 = yhat1+lambdas*dir1u
        
        yn2 = yn1
        yn2[Is] = yhat2
        
        pwforg = rep(0,length(y.thisloop))
        for (blocklcv in 1:CGGP$uoCOUNT) {
          IS = CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]]
          if(any(Is %in% IS) ){
            B = yn2[IS]
            rcpp_kronDBS(unlist(cholS.thisloop[gg+CGGP$uo[blocklcv,]]), B, CGGP$gridsizest[blocklcv,])
            pwforg[IS] = pwforg[IS]+CGGP$w[blocklcv] * B
          }
        }
        pw2 = pwforg
        pwforg = pwforg[Is]
        dir2 = pwforg
        
        s[,1] = yhat2-xsave[,2]
        ny[,1] = dir2-dirsave[,2]
        xsave[,1] = yhat2
        dirsave[,1] = dir2
        if(any(is.na(yn2))){
          break
        }
        
        dir1 = dir2
        yhat1 = yhat2
        yn1 = yn2
        pw1 = pw2
        
        L[lcv] = sum( pw2*yn2)
        if(lcv > (M+1)){
          if(max(abs(L[lcv:(lcv-3)]-L[(lcv-1):(lcv-4)])) < 10^(-3)*L[lcv]){
            break
          }
        }
        if(lcv > 10 && lcv <= M+1){
          if(max(abs(L[lcv:(lcv-3)]-L[(lcv-1):(lcv-4)])) < 10^(-3)*L[lcv]){
            break
          }
        }
      }
      
      
      if(!is.matrix(y)){
        yimputed <- yn1
      }else{
        yimputed[,oplcv] <- yn1
      }
    }
  }
  return(yimputed)
}