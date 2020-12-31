####==========================================================================####
## The R code for some specials functions, which are using in the computation   ##
## for the asymptotic covariances of TCF's at a specified cut point.            ##
##                                                                              ##
####==========================================================================####

## the derivative of function h(beta, theta, tau)
h_deriv <- function(thet, bet, n.nuis){
  res <- matrix(0, nrow = (6 + n.nuis), ncol = 3)
  res[, 1] <- c(bet[1]/thet[1]^2, 0, -1/thet[1], 0, 0, 0, rep(0,n.nuis))
  res[, 2] <- c(0, -(bet[2] - bet[3])/thet[2]^2, 0, 1/thet[2], -1/thet[2], 0,
                rep(0,n.nuis))
  res[, 3] <- c(bet[4]/(1 - thet[1]- thet[2])^2, bet[4]/(1 - thet[1]- thet[2])^2,
                0, 0, 0, 1/(1-thet[1]-thet[2]), rep(0,n.nuis))
  return(t(res))
}

## the first derivatives of rho.
rho_deriv <- function(u_matrix, rho_hat, ref_level = c("1", "2")){
  # u_matrix: the design matrix is used to fit the disease model.
  del_2 <- -u_matrix*rho_hat[,1]*rho_hat[,2]
  ref_level <- match.arg(ref_level)
  if(ref_level == "1"){
    del_1 <- u_matrix*rho_hat[,1]*(1 - rho_hat[,1])
    res <- cbind(del_1, del_2)
  }else if(ref_level == "2"){
    del_1 <- u_matrix*rho_hat[,2]*(1 - rho_hat[,2])
    res <- cbind(del_2, del_1)
  }else stop("Choose reference level in set {`1`, `2`}.")
  return(res)
}

## the first derivatives of 1/pi.
pi_inv_deriv <- function(u_matrix, pi_hat, pi_cov,
                       model = c("logit", "probit", "threshold")){
  # u_matrix: the design matrix is used to fit the vefication model.
  model <- match.arg(model)
  res <- switch(EXPR = model,
                logit = -u_matrix*(1 - pi_hat)/pi_hat,
                probit = -u_matrix*as.numeric(dnorm(u_matrix %*% pi_cov))/pi_hat,
                threshold = -u_matrix/pi_hat^2)
  return(res)
}

## Estimating function for rho_k = Pr(D_k = 1 | T, A).
mlogEstFunc <- function(dd, vv, u_matrix, rho_hat){
  temp <- u_matrix*vv
  B1 <- (dd[, 1] - rho_hat[, 1])*temp
  B2 <- (dd[, 2] - rho_hat[, 2])*temp
  return(cbind(B1, B2))
}

## Estimating function for pi = Pr(V = 1| T, A).
piEstFunc <- function(vv, u_matrix, pi_hat, pi_cov,
                      model = c("logit", "probit", "threshold")){
  model <- match.arg(model)
  res <- switch(EXPR = model,
                logit = (vv - pi_hat)*u_matrix,
                probit = u_matrix*as.numeric(dnorm(u_matrix %*% pi_cov))*
                  (vv - pi_hat)/(pi_hat*(1 - pi_hat)),
                threshold = (vv/pi_hat - 1)*u_matrix/(1 - pi_hat))
  return(res)
}

## Estimating function of Imputation Estimator (FI and MSI).
EstFuncIE <- function(dd, vv, tt, u_matrix, cc, cov.rho, rho.hat, thet, bet, m){
  n.par <- 6 + length(cov.rho)
  res <- matrix(0, nrow = nrow(u_matrix), ncol = n.par)
  res[,1] <- vv*(m*dd[,1] - thet[1] + (1 - m)*rho.hat[,1]) +
    (1 - vv)*(rho.hat[,1] - thet[1])
  res[,2] <- vv*(m*dd[,2] - thet[2] + (1 - m)*rho.hat[,2]) +
    (1 - vv)*(rho.hat[,2] - thet[2])
  res[,3] <- vv*((tt >= cc[1])*(m*dd[,1] + (1 - m)*rho.hat[,1]) - bet[1]) +
    (1 - vv)*((tt >= cc[1])*rho.hat[,1] - bet[1])
  res[,4] <- vv*((tt >= cc[1])*(m*dd[,2] + (1 - m)*rho.hat[,2]) - bet[2]) +
    (1 - vv)*((tt >= cc[1])*rho.hat[,2] - bet[2])
  res[,5] <- vv*((tt >= cc[2])*(m*dd[,2] + (1 - m)*rho.hat[,2]) - bet[3]) +
    (1 - vv)*((tt >= cc[2])*rho.hat[,2] - bet[3])
  res[,6] <- vv*((tt >= cc[2])*(m*dd[,3] + (1 - m)*rho.hat[,3]) - bet[4]) +
    (1 - vv)*((tt >= cc[2])*rho.hat[,3] - bet[4])
  res[,7:n.par] <- mlogEstFunc(dd, vv, u_matrix, rho.hat)
  return(res)
}

## Derivative of Estimating function of Imputation Estimator (FI and MSI).
EstFuncIE_deriv <- function(vv, tt, u_matrix, cc, cov.rho, rho.hat, rho.hes,
                            thet, bet, m){
  n.par <- 6 + length(cov.rho)
  der_rho_1 <- rho_deriv(u_matrix, rho.hat, ref_level = "1")
  der_rho_2 <- rho_deriv(u_matrix, rho.hat, ref_level = "2")
  der_rho_3 <- -(der_rho_1 + der_rho_2)
  temp1 <- (1 - m*vv)*der_rho_1
  temp2 <- (1 - m*vv)*der_rho_2
  A <- rbind(colSums(temp1), colSums(temp2))
  B1 <- colSums((tt >= cc[1])*temp1)
  B2 <- colSums((tt >= cc[1])*temp2)
  B3 <- colSums((tt >= cc[2])*temp2)
  B4 <- colSums((tt >= cc[2])*(1 - m*vv)*der_rho_3)
  B <- rbind(B1, B2, B3, B4)
  res.left <- rbind(diag(-length(tt), 6, 6),
                    matrix(0, nrow = length(cov.rho), ncol = 6))
  res.right <- rbind(A, B, -rho.hes)
  res <- cbind(res.left, res.right)
  return(res)
}

## Estimating function of IPW Estimator
EstFuncIPW <- function(dd, vv, tt, u_matrix, cc, cov.pi, pi.hat, thet, bet,
                       model){
  n.par <- 6 + length(cov.pi)
  res <- matrix(0, nrow = length(tt), ncol = n.par)
  res[,1] <- vv*(dd[,1] - thet[1])/pi.hat
  res[,2] <- vv*(dd[,2] - thet[2])/pi.hat
  res[,3] <- vv*((tt >= cc[1])*dd[,1] - bet[1])/pi.hat
  res[,4] <- vv*((tt >= cc[1])*dd[,2] - bet[2])/pi.hat
  res[,5] <- vv*((tt >= cc[2])*dd[,2] - bet[3])/pi.hat
  res[,6] <- vv*((tt >= cc[2])*dd[,3] - bet[4])/pi.hat
  res[,7:n.par] <- piEstFunc(vv, u_matrix, pi.hat, cov.pi, model)
  return(res)
}

## Derivative of Estimating function of IPW Estimator
EstFuncIPW_deriv <- function(dd, vv, tt, u_matrix, cc, cov.pi, pi.hat, pi.hes,
                             thet, bet, model){
  n.par <- 6 + length(cov.pi)
  der_pi_inv <- pi_inv_deriv(u_matrix, pi.hat, cov.pi, model)
  term <- vv*der_pi_inv
  A1 <- colSums(term*(dd[,1] - thet[1]))
  A2 <- colSums(term*(dd[,2] - thet[2]))
  B11 <- colSums(term*((tt >= cc[1])*dd[,1] - bet[1]))
  B12 <- colSums(term*((tt >= cc[1])*dd[,2] - bet[2]))
  B22 <- colSums(term*((tt >= cc[2])*dd[,2] - bet[3]))
  B23 <- colSums(term*((tt >= cc[2])*dd[,3] - bet[4]))
  res.left <- rbind(diag(-sum(vv/pi.hat), 6, 6),
                    matrix(0, nrow = length(cov.pi), ncol = 6))
  res.right <- rbind(A1, A2, B11, B12, B22, B23, -pi.hes)
  res <- cbind(res.left, res.right)
  return(res)
}

## Estimating function of SPE Estimator
EstFuncSPE <- function(dd, vv, tt, u_matrix_rho, u_matrix_pi, cc, cov.rho,
                       rho.hat, cov.pi, pi.hat, thet, bet, model){
  n.par <- 6 + length(cov.rho) + length(cov.pi)
  res <- matrix(0, nrow = length(tt), ncol = n.par)
  res[,1] <- vv*(dd[,1] - thet[1])/pi.hat -
    (rho.hat[,1] - thet[1])*(vv - pi.hat)/pi.hat
  res[,2] <- vv*(dd[,2] - thet[2])/pi.hat -
    (rho.hat[,2] - thet[2])*(vv - pi.hat)/pi.hat
  res[,3] <- vv*((tt >= cc[1])*dd[,1] - bet[1])/pi.hat -
    (vv - pi.hat)*((tt >= cc[1])*rho.hat[,1] - bet[1])/pi.hat
  res[,4] <- vv*((tt >= cc[1])*dd[,2] - bet[2])/pi.hat -
    (vv - pi.hat)*((tt >= cc[1])*rho.hat[,2] - bet[2])/pi.hat
  res[,5] <- vv*((tt >= cc[2])*dd[,2] - bet[3])/pi.hat -
    (vv - pi.hat)*((tt >= cc[2])*rho.hat[,2] - bet[3])/pi.hat
  res[,6] <- vv*((tt >= cc[2])*dd[,3] - bet[4])/pi.hat -
    (vv - pi.hat)*((tt >= cc[2])*rho.hat[,3] - bet[4])/pi.hat
  res[,7:(6 + length(cov.rho))] <- mlogEstFunc(dd, vv, u_matrix_rho, rho.hat)
  res[,(7 + length(cov.rho)):n.par] <- piEstFunc(vv, u_matrix_pi, pi.hat, cov.pi,
                                                model)
  return(res)
}

## Derivative of Estimating function of SPE Estimator
EstFuncSPE_deriv <- function(dd, vv, tt, u_matrix_rho, u_matrix_pi, cc, cov.rho,
                             rho.hat, rho.hes, cov.pi, pi.hat, pi.hes, thet, bet,
                             model){
  n.par <- 6 + length(cov.rho) + length(cov.pi)
  der_rho_1 <- rho_deriv(u_matrix_rho, rho.hat, ref_level = "1")
  der_rho_2 <- rho_deriv(u_matrix_rho, rho.hat, ref_level = "2")
  der_rho_3 <- -(der_rho_1 + der_rho_2)
  temp.v <- (1- vv/pi.hat)
  temp1 <- temp.v*der_rho_1
  temp2 <- temp.v*der_rho_2
  der_pi_inv <- pi_inv_deriv(u_matrix_pi, pi.hat, cov.pi, model)
  temp3 <- vv*der_pi_inv
  H <- rbind(colSums(temp1), colSums(temp2))
  G1 <- colSums((tt >= cc[1])*temp1)
  G2 <- colSums((tt >= cc[1])*temp2)
  G3 <- colSums((tt >= cc[2])*temp2)
  G4 <- colSums((tt >= cc[2])*temp.v*der_rho_3)
  G <- rbind(G1, G2, G3, G4)
  D1 <- colSums(temp3*(rho.hat[,1] - dd[,1]))
  D2 <- colSums(temp3*(rho.hat[,2] - dd[,2]))
  E11 <- colSums(temp3*((tt >= cc[1])*(rho.hat[,1] - dd[,1])))
  E12 <- colSums(temp3*((tt >= cc[1])*(rho.hat[,2] - dd[,2])))
  E22 <- colSums(temp3*((tt >= cc[2])*(rho.hat[,2] - dd[,2])))
  E23 <- colSums(temp3*((tt >= cc[2])*(rho.hat[,3] - dd[,3])))
  res.left <- rbind(diag(-length(tt), 6, 6),
                    matrix(0, nrow = length(cov.rho) + length(cov.pi), ncol = 6))
  res.right.rho <- rbind(H, G, -rho.hes,
                         diag(0, nrow = length(cov.pi), ncol = length(cov.rho)))
  res.right.pi <- rbind(D1, D2, E11, E12, E22, E23,
                        diag(0, nrow = length(cov.rho), ncol = length(cov.pi)),
                        -pi.hes)
  res <- cbind(res.left, res.right.rho, res.right.pi)
  return(res)
}
