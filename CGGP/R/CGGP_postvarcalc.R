#' Calculate posterior variance, faster version
#'
#' @param GMat Matrix
#' @param dGMat Derivative of matrix
#' @param cholS Cholesky factorization of S
#' @param dSMat Deriv of SMat
#' @param INDSN Indices, maybe
#' @param numpara Number of parameters for correlation function
#' @param returnlogs Should log scale be returned
#' @param returnderiratio Should derivative ratio be returned?
#' @param returndG Should dG be returned
#' @param returndiag Should diag be returned
#' @param ... Placeholder
#'
#' @return Variance posterior
# @export
#' @noRd
CGGP_internal_postvarmatcalc_fromGMat <- function(GMat, dGMat,cholS,dSMat,INDSN,numpara,...,
                                         returnlogs=FALSE, returnderiratio =FALSE,
                                         returndG = FALSE,returndiag = FALSE) {
  # Next line was giving error with single value, so I changed it
  CoinvC1o = backsolve(cholS,backsolve(cholS, t(GMat[,INDSN, drop=F]), transpose = TRUE))
  # backsolve1 <- backsolve(cholS,if (is.matrix(GMat[,INDSN]) && ncol(GMat[,INDSN])>1) t(GMat[,INDSN]) else GMat[,INDSN], transpose = TRUE)
  # CoinvC1o = backsolve(cholS,backsolve1)
  if(returndiag){
    if(!returnlogs){
      Sigma_mat = rowSums(t((CoinvC1o))*((GMat[,INDSN])))
    }else{
      nlSm = rowSums(t((CoinvC1o))*((GMat[,INDSN])))
      Sigma_mat =log(nlSm)
    }
    np = length(Sigma_mat)
  }else{
    if(!returnlogs){
      Sigma_mat = t(CoinvC1o)%*%(t(GMat[,INDSN]))
    }else{
      nlSm =(t(CoinvC1o))%*%(t(GMat[,INDSN]))
      Sigma_mat =log(nlSm)
    }
    np = dim(Sigma_mat)[1]
  }
  
  
  if(returndG){
    nb = dim(GMat)[2]
    nc = dim(dSMat)[1]
    if(returndiag){
      dSigma_mat = matrix(0,np,numpara)
    }else{
      dSigma_mat = matrix(0,np,np*numpara)
    }
    
    for(k in 1:numpara){
      dS = dSMat[,(k-1)*nc+(1:nc)]
      CoinvC1oE = ((as.matrix(dS))%*%t(GMat[,INDSN]))
      if(returndiag){
        dCoinvC1o = backsolve(cholS,backsolve(cholS,CoinvC1oE, transpose = TRUE))
        dCoinvC1o = dCoinvC1o + backsolve(cholS,backsolve(cholS,t(dGMat[,(k-1)*nb+INDSN]), transpose = TRUE))
        
        if(!returnlogs && !returnderiratio){
          dSigma_mat[,k] = rowSums(t(dCoinvC1o)*(GMat[,INDSN]))+rowSums(t(CoinvC1o)*(dGMat[,(k-1)*nb+INDSN]))
        }else if(returnderiratio){
          dSigma_mat[,k] = rowSums(t(dCoinvC1o)*(GMat[,INDSN]))+rowSums(t(CoinvC1o)*(dGMat[,(k-1)*nb+INDSN]))/Sigma_mat
        }else{
          dSigma_mat[,k] = (rowSums(t(dCoinvC1o)*(GMat[,INDSN]))+rowSums(t(CoinvC1o)*(dGMat[,(k-1)*nb+INDSN])))/nlSm
        }
      }else{
        dCoinvC1_part1 = t(CoinvC1o)%*%(CoinvC1oE)
        dCoinvC1_part2 = t(CoinvC1o)%*%t(dGMat[,(k-1)*nb+INDSN])
        dCoinvC1_part2 = dCoinvC1_part2+t(dCoinvC1_part2)
        if(!returnlogs && !returnderiratio){
          dSigma_mat[,(k-1)*np + 1:np] =dCoinvC1_part1+dCoinvC1_part2
        }else if(returnderiratio){
          dSigma_mat[,(k-1)*np + 1:np] =( dCoinvC1_part1+dCoinvC1_part2)/Sigma_mat
        }else{
          dSigma_mat[,(k-1)*np + 1:np] =( dCoinvC1_part1+dCoinvC1_part2)/nlSm
        }
      }
    }
    return(list("Sigma_mat"= Sigma_mat,"dSigma_mat" = dSigma_mat))
  }else{
    return(Sigma_mat)
  }
}