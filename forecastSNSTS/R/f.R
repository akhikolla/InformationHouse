################################################################################
#' Compute \eqn{f(\delta)} for a tvAR(p) process
#'
#' This functions computes the quantity \eqn{f(\delta)} defined in (24) of
#' Kley et al. (2019) when the underlying process follows an tvAR(p) process.
#' Recall that, to apply Theorem 3.1 in Kley et al. (2019), the function
#' \eqn{f(\delta)} is required to be positive, which can be verified with the
#' numbers returned from this function.
#' The function returns a vector with elements \eqn{f(\delta)} 
#' for each \eqn{\delta} in \code{which.deltas}, with \eqn{f(\delta)}
#' defined as  
#' \deqn{f(\delta) := \min_{p_1,p_2 = 0, \ldots, p_{\max}} \min_{N \in \mathcal{N}} \Big| {\rm MSPE}_{s_1/T,m/T}^{(p_1,h)}(\frac{s_1}{T}) - (1+\delta) \cdot {\rm MSPE}_{N/T,m/T}^{(p_2,h)}(\frac{s_1}{T}) \Big|, \quad \delta \geq 0}
#' where \eqn{T, m, p_{\max}, h} are positive integers,
#' \eqn{\mathcal{N} \subset \{p_{\max}+1, \ldots, T-m-h\}}, and \eqn{s_1 := T-m-h+1}.
#' 
#' The function \eqn{{\rm MSPE}_{\Delta_1, \Delta_2}^{(p,h)}(u)} is defined, for real-valued \eqn{u} and
#' \eqn{\Delta_1, \Delta_2 \geq 0}, in terms of the second order properties of the process:
#' \deqn{{\rm MSPE}_{\Delta_1, \Delta_2}^{(p,h)}(u) := \int_0^1 g^{(p,h)}_{\Delta_1}\Big( u + \Delta_2 (1-x) \Big) {\rm d}x,}
#' with \eqn{g^{(0,h)}_{\Delta}(u) := \gamma_0(u)} and, for \eqn{p = 1, 2, \ldots},
#' \deqn{g^{(p,h)}_{\Delta}(u) := \gamma_0(u) - 2 \big( v_{\Delta}^{(p,h)}(u) \big)' \gamma_0^{(p,h)}(u) + \big( v_{\Delta}^{(p,h)}(u) \big)' \Gamma_0^{(p)}(u) v_{\Delta}^{(p,h)}(u)}
#' \deqn{\gamma_0^{(p,h)}(u) := \big( \gamma_h(u), \ldots, \gamma_{h+p-1}(u) \big)',}
#' where
#' \deqn{v^{(p,h)}_{\Delta}(u) := e'_1 \big( e_1 \big( a_{\Delta}^{(p)}(t) \big)' + H \big)^h,}
#' with \eqn{e_1} and \eqn{H} defined in the documentation of \code{\link{predCoef}} and,
#' for every real-valued \eqn{u} and \eqn{\Delta \geq 0},
#' \deqn{a^{(p)}_{\Delta}(u) := \Gamma^{(p)}_{\Delta}(u)^{-1} \gamma^{(p)}_{\Delta}(u),}
#' where
#' \deqn{\gamma^{(p)}_{\Delta}(u) := \int_0^1 \gamma^{(p)}(u+\Delta (x-1)) {\rm d}x, \quad \gamma^{(p)}(u) := [\gamma_1(u)\;\ldots\;\gamma_p(u)]',}
#' \deqn{\Gamma^{(p)}_{\Delta}(u) := \int_0^1 \Gamma^{(p)}(u+\Delta (x-1)) {\rm d}x, \quad \Gamma^{(p)}(u) := (\gamma_{i-j}(u);\,i,j=1,\ldots,p).}
#'
#' The local autocovariances \eqn{\gamma_k(u)} are defined as the lag-\eqn{k}
#' autocovariances of an AR(p) process which has coefficients
#' \eqn{a_1(u), \ldots, a_p(u)} and innovations with variance \eqn{\sigma(u)^2},
#' because the underlying model is assumed to be tvAR(p)  
#' \deqn{Y_{t,T} = \sum_{j=1}^p a_j(t/T) Y_{t-j,T} + \sigma(t/T) \varepsilon_{t},}
#' where \eqn{a_1, \ldots, a_p} are real valued functions (defined on \eqn{[0,1]}) and \eqn{\sigma} is a
#' positive function (defined on \eqn{[0,1]}).
#' 
#' @name f
#' @export
#' 
#' @importFrom stats integrate
#' 
#' @param which.deltas vector containing the \eqn{\delta}'s for which to
#'                     to compute \eqn{f(\delta)},
#' @param p_max        parameter \eqn{p_{\max}},
#' @param h            parameter \eqn{h},
#' @param T            parameter \eqn{T},
#' @param Ns           a vector containing the elements of the set
#'                     \eqn{\mathcal{N}},
#' @param m            parameter \eqn{m},
#' @param a            a list of real-valued functions, specifying the
#'                     coefficients of the tvAR(p) process,
#' @param sigma        a positive-valued function, specifying the variance
#'                     of the innovations of the tvAR(p) process,
#' 
#' @return Returns     a vector with the values \eqn{f(\delta)}, as defined in
#'                     (24) of Kley et al. (2019), where it is now denoted by \eqn{q(\delta)}, for each \eqn{\delta} in
#' 										 \code{which.delta}.
#'
#' @examples
#' \dontrun{
#' ## because computation is quite time-consuming.
#' n <- 100
#' a <- list( function(u) {return(0.8+0.19*sin(4*pi*u))} )
#' sigma <- function (u) {return(1)}
#' 
#' Ns <- seq( floor((n/2)^(4/5)), floor(n^(4/5)),
#'            ceiling((floor(n^(4/5)) - floor((n/2)^(4/5)))/25) )
#' which.deltas <- c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.4, 0.6)
#' P_max <- 7
#' H <- 1
#' m <- floor(n^(.85)/4)
#' 
#' # now replicate some results from Table 4 in Kley et al. (2019)
#' f( which.deltas, P_max, h = 1, n - m, Ns, m, a, sigma )
#' f( which.deltas, P_max, h = 5, n - m, Ns, m, a, sigma )
#' }
################################################################################

f <- function(which.deltas, p_max, h, T, Ns, m, a, sigma) {
  
  gamma <- function(k, u) {
    ## create a(u)
    au <- c()
    for (a_spec in a) {
      au <- c(au, a_spec(u))
    }
    acfARp(au, sigma(u), k )
  }
  
  ## Define the function g_{D1}^{(p,h)} from (49) of Kley et al (2017)
  g <- function( D1, p, h, u) {
    
    res <- gamma(0, u)
    
    if (p > 0) {
      
      gamma_D1 <- function(k,u) {
        aux <- Vectorize(function(x) { return( gamma(k, u + D1 * (x-1) ) ) })
        return (integrate(aux, 0, 1)$value)
      }
      
      gamma0_all <- Vectorize(function(k) {gamma(k, u)})(0:p)
      gamma_all  <- Vectorize(function(k) {gamma_D1(k, u)})(0:p)
      
      gammaP <- gamma_all[(1+1):(p+1)]
      
      Gamma0P <- matrix(gamma0_all[0+1], ncol = 1, nrow = 1)
      GammaP <- matrix(gamma_all[0+1], ncol = 1, nrow = 1)
      
      if (p >= 2) {
        for (i in 2:p) {
          
          right  <- matrix(gamma0_all[(i-1+1):(1+1)], ncol=1)
          bottom <- matrix(gamma0_all[(i-1+1):(0+1)], nrow=1)
          Gamma0P <- rbind(cbind( Gamma0P, right), bottom)
          
          right  <- matrix(gamma_all[(i-1+1):(1+1)], ncol=1)
          bottom <- matrix(gamma_all[(i-1+1):(0+1)], nrow=1)
          GammaP <- rbind(cbind( GammaP, right), bottom)
          
        }
      }
      
      YW <- solve(GammaP) %*% gammaP
      
      A_p <- rbind(t(YW), cbind(diag(rep(1, p-1)), rep(0,p-1)) )
      
      V_p <- diag(rep(1,p))
      for (i in 1:h) {
        V_p <- V_p %*% A_p
      }
      
      vP <- t(c(1,rep(0,p-1)) %*% V_p)
      
      res <- res - 2 * t(vP) %*% Vectorize(function(k) {gamma(k, u)})(h:(h+p-1))
      res <- res + t(vP) %*% Gamma0P %*% vP
      
    }
    
    return(res)
    
  } ## end definition of g
  
  ## Define the function  MSPE_{D1, D2}^{(p,h)}(u)
  MSPE <- function(D1, D2, p, h, u) {
    aux <- Vectorize(function(x) { return( g(D1, p, h, u + D2 * (1-x) ) ) })
    return (integrate(aux, 0, 1)$value)
  }
  
  M_s  <- array(0, dim = c(p_max + 1) )
  M_ls <- array(0, dim = c(p_max + 1, length(Ns)) ) 
  
  for (p1 in 0:p_max) {
    s1 <- T - m - h
    M_s[p1+1] <- MSPE( s1/T, m/T, p1, h, s1/T)
  }
  
  for (p2 in 0:p_max) {
    for (N.i in 1:length(Ns)) {
      M_ls[p2+1, N.i] <- MSPE( Ns[N.i]/T, m/T, p2, h, s1/T)
    }
  }
  
  fdelta <- array(10^10, dim = c( length(which.deltas), p_max ) )
  
  for (p in 1:p_max) {
    for (p1 in 0:p) {
      for (p2 in 0:p) {
        for (N.i in 1:length(Ns)) {
          for (delta.i in 1:length(which.deltas)) {
            dif <- abs(M_s[p1+1] - (1+which.deltas[delta.i]) * M_ls[p2+1, N.i])
            if (dif < fdelta[delta.i, p]) {
              fdelta[delta.i, ] <- dif
            }
          }
        }
      }
    }
  }
  
  return(fdelta[,p_max])

}