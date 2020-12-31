####==========================================================================####
## The R code for estimating asymptotic covariance matrix of KNN estimators     ##
##                                                                              ##
##                                                                              ##
####==========================================================================####

psknn <- function(x, y, type = c("eucli", "manha", "canber", "lagran", "mahala")){
  n <- nrow(x)
  S.inv <- solve(cov(x))
  type <- match.arg(type)
  dista <- switch(type,
                  eucli = function(xi, xx, Sx.inv){
                    ss <- colSums((xi - t(xx))^2)
                    return(sqrt(ss))
                  },
                  manha = function(xi, xx, Sx.inv){
                    ss <- colSums(abs(xi - t(xx)))
                    return(ss)
                  },
                  canber = function(xi, xx, Sx.inv){
                    ss <- colSums( abs(xi - t(xx))/( abs(xi) + t(abs(xx)) ) )
                    return(ss)
                  },
                  lagran = function(xi, xx, Sx.inv){
                    ss <- apply(abs(xi - t(xx)), 2, max)
                    return(ss)
                  },
                  mahala = function(xi, xx, Sx.inv){
                    ss1 <- (xi - t(xx))
                    ss <- diag(t(ss1)%*%Sx.inv%*%ss1)
                    return(sqrt(ss))
                  }
  )
  res <- numeric(n)
  K1 <- rep(0, n) # rep(1,n)
  for(i in 1:n){
    dist.tem <- dista(x[i,], x, S.inv)
    id.tem <- order(dist.tem)[-1]
    y.temp <- y[id.tem]
    if(y[i] == 0){
      if(y.temp[1] == 0){
        K1[i] <- which((1 - y.temp) == 0)[1]
      }
      else{
        K1[i] <- which((1 - y.temp) == 1)[1]
      }
    }
    else{
      if(y.temp[1] == 1){
        K1[i] <- which((1 - y.temp) == 1)[1]
      }
      else{
        K1[i] <- which((1 - y.temp) == 0)[1]
      }
    }
    res[i] <- mean(y.temp[1:K1[i]])
  }
  return(res)
}

# ### Computing the variances
asy.Cov.KNN <- function(X, thet, bet, cc, pi.hat, rho.hat, k){
  n <- length(X)
  res <- matrix(0, 3, 3)
  bet[5] <- bet[2] - bet[3]
  ### Computing the variances
  omega1 <- mean((rho.hat[,1]*(1 - rho.hat[,1])*(1 - pi.hat))*
                  ((k + 1)/k + (1 - pi.hat)/pi.hat))
  omega2 <- mean((rho.hat[,2]*(1 - rho.hat[,2])*(1 - pi.hat))*
                  ((k + 1)/k + (1 - pi.hat)/pi.hat))
  omega3 <- mean((rho.hat[,3]*(1 - rho.hat[,3])*(1 - pi.hat))*
                  ((k + 1)/k + (1 - pi.hat)/pi.hat))
  sig.sq <- (thet*(1 - thet) + c(omega1, omega2, omega3))
  omega11 <- mean(((X >= cc[1])*rho.hat[,1]*(1 - rho.hat[,1])*
                    (1 - pi.hat))*((k + 1)/k + (1 - pi.hat)/pi.hat))
  omega12 <- mean((X >= cc[1])*(rho.hat[,2]*(1 - rho.hat[,2])*
                                 (1 - pi.hat))*((k + 1)/k + (1 - pi.hat)/pi.hat))
  omega22 <- mean((X >= cc[2])*(rho.hat[,2]*(1 - rho.hat[,2])*
                                 (1 - pi.hat))*((k + 1)/k + (1 - pi.hat)/pi.hat))
  omega23 <- mean((X >= cc[2])*(rho.hat[,3]*(1 - rho.hat[,3])*
                                 (1 - pi.hat))*((k + 1)/k + (1 - pi.hat)/pi.hat))
  omega <- omega12 - omega22
  sig.sq.bet <- (bet*(1 - bet) + c(omega11, omega12, omega22, omega23, omega))
  gam <- sig.sq.gam <- numeric(4)
  gam[1] <- thet[1] - bet[1]
  gam[2] <- thet[2] - bet[2]
  gam[3] <- thet[2] - bet[3]
  gam[4] <- thet[3] - bet[4]
  sig.sq.gam[1] <- (gam[1]*(1 - gam[1]) + omega1 - omega11)
  sig.sq.gam[2] <- (gam[2]*(1 - gam[2]) + omega2 - omega12)
  sig.sq.gam[3] <- (gam[3]*(1 - gam[3]) + omega2 - omega22)
  sig.sq.gam[4] <- (gam[4]*(1 - gam[4]) + omega3 - omega23)
  sig111 <- (sig.sq[1] + sig.sq.bet[1] - sig.sq.gam[1])/2
  sig212 <- (sig.sq[2] + sig.sq.bet[2] - sig.sq.gam[2])/2
  sig222 <- (sig.sq[2] + sig.sq.bet[3] - sig.sq.gam[3])/2
  sig323 <- (sig.sq[3] + sig.sq.bet[4] - sig.sq.gam[4])/2
  res[1,1] <- bet[1]^2*sig.sq[1]/thet[1]^4 + sig.sq.bet[1]/thet[1]^2 -
    2*bet[1]*sig111/thet[1]^3
  res[2,2] <- bet[5]^2*sig.sq[2]/thet[2]^4 + sig.sq.bet[5]/thet[2]^2 -
    2*bet[5]*(sig212 - sig222)/thet[2]^3
  res[3,3] <- bet[4]^2*sig.sq[3]/thet[3]^4 + sig.sq.bet[4]/thet[3]^2 -
    2*bet[4]*sig323/thet[3]^3
  ### Estimating xi_12
  omeg_sig12 <- mean((rho.hat[,1]*rho.hat[,2]*(1 - pi.hat))*
                      ((k + 1)/k + (1 - pi.hat)/pi.hat))
  omeg1112 <- mean((X >= cc[1])*(rho.hat[,1]*rho.hat[,2]*
                            (1 - pi.hat))*((k + 1)/k + (1 - pi.hat)/pi.hat))
  omeg1122 <- mean((X >= cc[2])*(rho.hat[,1]*rho.hat[,2]*
                            (1 - pi.hat))*((k + 1)/k + (1 - pi.hat)/pi.hat))
  omeg11.12_22 <- omeg1112 - omeg1122
  sig.11.12_22 <- -(bet[1]*bet[5] + omeg11.12_22)
  sig.1.12_22 <- -(thet[1]*bet[5] + omeg11.12_22)
  sig.2.11 <- -(thet[2]*bet[1] + omeg1112)
  sig.thet.12 <- -(thet[1]*thet[2] + omeg_sig12)
  res[1,2] <- -sig.11.12_22/(thet[1]*thet[2]) + bet[1]*sig.1.12_22/
                  (thet[1]^2 *thet[2]) - bet[5]*(bet[1]*sig.thet.12/thet[1]^2 -
                  sig.2.11/thet[1])/thet[2]^2
  ###Estimating xi_13
  omeg1123 <- mean((X >= cc[2])*(rho.hat[,1]*rho.hat[,3]*
                            (1 - pi.hat))*((k + 1)/k + (1 - pi.hat)/pi.hat))
  omeg311 <- mean((X >= cc[1])*(rho.hat[,1]*rho.hat[,3]*
                            (1 - pi.hat))*((k + 1)/k + (1 - pi.hat)/pi.hat))
  sig.1.23 <- -(thet[1]*bet[4] + omeg1123)
  sig.11.23 <- -(bet[1]*bet[4] + omeg1123)
  sig.3.11 <- -(thet[3]*bet[1] + omeg311)
  res[1,3] <- (bet[1]*sig.1.23/thet[1] - sig.11.23)/(thet[1]*thet[3]) +
              bet[4]*(bet[1]*(sig.sq[1] + sig.thet.12)/thet[1] + sig.3.11)/
              (thet[1]*thet[3]^2)
  ###Estimating xi_23
  omeg312 <- mean((X >= cc[1])*(rho.hat[,2]*rho.hat[,3]*
                          (1 - pi.hat))*((k + 1)/k + (1 - pi.hat)/pi.hat))
  omeg322 <- mean((X >= cc[2])*(rho.hat[,2]*rho.hat[,3]*
                          (1 - pi.hat))*((k + 1)/k + (1 - pi.hat)/pi.hat))
  omeg3.12_22 <- omeg312 - omeg322
  sig.23.12_22 <- -bet[4]*bet[5]
  sig.2.23 <- -(thet[2]*bet[4] + omeg322)
  sig.3.12_22 <- -(thet[3]*bet[5] + omeg3.12_22)
  res[2,3] <- (sig.23.12_22 - bet[5]*sig.2.23/thet[2])/(thet[3]*thet[2]) +
              bet[4]*(-sig.3.12_22 - bet[5]*(sig.sq[2] + sig.thet.12)/
              thet[2])/(thet[2]*thet[3]^2)
  res[lower.tri(res)] <- res[upper.tri(res)]
  return(res/n)
}
