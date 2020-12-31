#' This function is to compute the centralized instrument given the covariates L
#' @keywords internal
#' @export
IV_center <- function (Z, L = NULL) {
  
  # Now we obtain the iid-decomp of \hat\theta
  expit <- function (x) {
    exp(x)/(1 + exp(x))
  }
  
  IV_type <- all(Z == 0 | Z == 1)

  # fitting model for expectation of instrument

  if (IV_type) {
    if (is.null(L)) {
      zmod <- glm(Z ~ 1, family = "binomial")
      E.dot <- matrix(rep(1, length(Z)), ncol = 1) * c(fitted(zmod) * (1 - fitted(zmod))) 
    }
    else {
      zmod <- glm(Z ~ L, family = "binomial")
      E.dot <- cbind(1, L) * c(fitted(zmod) * (1 - fitted(zmod))) 
    }
    
  }
  else{
    if (is.null(L)) {
      zmod <- lm(Z ~ 1)
      E.dot <- matrix(rep(1, length(Z)), ncol = 1)
    }
    else {
      zmod <- lm(Z ~ L)
      E.dot <- cbind(1, L)
    }
    
  }
  
  eps.theta <- as.matrix(iid(zmod))
  p.dim <- dim(eps.theta)[2]
  Z.c <- Z - fitted(zmod)
  
  iv_center <- list(Zc = Z.c,
                    epstheta = eps.theta,
                    Edot = E.dot,
                    pdim = p.dim)
  return(iv_center)
}

#' This function is to generate a grid of treatment status
#' @keywords internal
#' @export
treatment_status <- function (n,
                              K,
                              stime, 
                              treatment_init, 
                              treatment_shift_time, 
                              max.time) {

  D_status <- matrix(0, ncol = K, nrow = n)
  for (i in 1:n) {
    if (treatment_shift_time[i] <= 0 || treatment_shift_time[i] >= max.time) {
      D_status[i, ] = treatment_init[i]
      next
    }
    D_status[i, stime < treatment_shift_time[i]] = treatment_init[i]
    D_status[i, stime >= treatment_shift_time[i]] = 1 - treatment_init[i]
  }
  return(D_status)
}
