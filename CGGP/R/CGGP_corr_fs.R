#' Cauchy correlation function
#' 
#' Calculate correlation matrix for two sets of points in one dimension.
#' Note that this is not the correlation between two vectors.
#'
#' @param x1 Vector of coordinates from same dimension
#' @param x2 Vector of coordinates from same dimension
# ' @param LS Log of parameter that controls lengthscale
# ' @param FD Logit of 0.5*parameter  that controls the fractal dimension
# ' @param HE Log of parameter that controls the hurst effect
#' @param theta Correlation parameters:
#' \itemize{
#'     \item LS Log of parameter that controls lengthscale
#'     \item FD Logit of 0.5*parameter  that controls the fractal dimension
#'     \item HE Log of parameter that controls the hurst effect
#' }
#' @param return_dCdtheta Should dCdtheta be returned?
#' @param return_numpara Should it just return the number of parameters?
#' @param returnlogs Should log of correlation be returned?
#'
#' @return Matrix of correlation values between x1 and x2
#' @export
#' @family correlation functions
#'
#' @examples
#' CGGP_internal_CorrMatCauchy(c(0,.2,.4),c(.1,.3,.5), theta=c(-1,.9,.1))
CGGP_internal_CorrMatCauchy <- function(x1, x2, theta, return_dCdtheta=FALSE,
                                        return_numpara=FALSE,
                                        returnlogs=FALSE) {
  if(return_numpara){
    return(3)
  }else{ 
    if (length(theta) != 3) {stop("CorrMatCauchy theta should be length 3")}
    diffmat =abs(outer(x1,x2,'-')); 
    
    expLS = exp(3*(theta[1]))
    expHE = exp(3*(theta[2]))
    h = diffmat/expLS
    alpha = 2*exp(3*theta[3]+2)/(1+exp(3*theta[3]+2))
    halpha = h^alpha
    pow = -expHE/alpha
    
    if (!returnlogs) {
      C = (1+halpha)^pow
    } else {
      C = pow * log(1+halpha)
    }
    if(return_dCdtheta){
      if (!returnlogs) {
        dCdtheta = cbind(3*expHE*((1+halpha)^(pow-1))*(halpha),
                         3*C*pow*log(1+halpha),
                         (C*(expHE*log(1+halpha)/alpha^2 - 
                               expHE*halpha*log(h)/alpha/(1+halpha))) * 
                           6*exp(3*theta[3]+2)/(1+exp(3*theta[3]+2))^2)
      } else {
        dCdtheta = cbind(3*expHE*halpha/(1+halpha),
                         3*pow*log(1+halpha),
                         ((expHE*log(1+halpha)/alpha^2 - 
                             expHE*halpha*log(h)/alpha/(1+halpha))) * 
                           6*exp(3*theta[3]+2)/(1+exp(3*theta[3]+2))^2)
      }
      dCdtheta[is.na(dCdtheta)] = 0
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }
}

#' CauchySQT correlation function
#' 
#' Calculate correlation matrix for two sets of points in one dimension.
#' Note that this is not the correlation between two vectors.
#'
#' @inheritParams CGGP_internal_CorrMatCauchy
#'
#' @return Matrix of correlation values between x1 and x2
#' @export
#' @family correlation functions
#'
#' @examples
#' CGGP_internal_CorrMatCauchySQT(c(0,.2,.4),c(.1,.3,.5), theta=c(-.1,.3,-.7))
CGGP_internal_CorrMatCauchySQT <- function(x1, x2,theta, return_dCdtheta = FALSE,
                                           return_numpara=FALSE,
                                           returnlogs=FALSE) {
  if(return_numpara){
    return(3);
  } else{
    if (length(theta) != 3) {stop("CorrMatCauchySQT theta should be length 3")}
    expTILT = exp((theta[3]))
    expLS = exp(3*(theta[1]))
    x1t = (x1+10^(-2))^expTILT
    x2t = (x2+10^(-2))^expTILT
    x1ts = x1t/expLS
    x2ts = x2t/expLS
    
    diffmat =abs(outer(x1ts,x2ts,'-')); 
    expHE = exp(3*(theta[2]))
    h = diffmat
    alpha = 2*exp(5)/(1+exp(5))
    halpha = h^alpha
    pow = -expHE/alpha
    
    if (!returnlogs) {
      C = (1+halpha)^pow
    } else {
      C = pow * log(1+halpha)
    }
    
    if(return_dCdtheta){
      Q = ((1+halpha)^(pow-1))
      gt1 = x1t*log(x1+10^(-2))
      gt2 = x2t*log(x2+10^(-2))
      lh =outer(gt1,gt2,'-')
      hnabs = outer(x1ts,x2ts,'-')
      LO = alpha*expTILT*(pow/expLS)*(abs(h)^(alpha-1)*lh*sign(hnabs))
      if (!returnlogs) {
        dCdtheta = cbind(3*expHE*((1+halpha)^(pow-1))*(halpha),3*C*pow*log(1+halpha),LO*Q)
      } else {
        dCdtheta = cbind(3*expHE*halpha/(1+halpha), 3*pow*log(1+halpha), LO/(1+halpha))
      }
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }
}



#' CauchySQ correlation function
#' 
#' Calculate correlation matrix for two sets of points in one dimension
#' Note that this is not the correlation between two vectors.
#'
#' @inheritParams CGGP_internal_CorrMatCauchy
#'
#' @return Matrix of correlation values between x1 and x2
#' @export
#' @family correlation functions
#'
#' @examples
#' CGGP_internal_CorrMatCauchySQ(c(0,.2,.4),c(.1,.3,.5), theta=c(-.7,-.5))
CGGP_internal_CorrMatCauchySQ <- function(x1, x2,theta, return_dCdtheta = FALSE,
                                          return_numpara =FALSE,returnlogs = FALSE) {
  if(return_numpara){
    return(2);
  }else{ 
    if (length(theta) != 2) {stop("CorrMatCauchySQ theta should be length 2")}
    diffmat =abs(outer(x1,x2,'-')); 
    
    expLS = exp(3*theta[1])
    expHE = exp(3*theta[2])
    h = diffmat/expLS
    alpha = 2*exp(0+6)/(1+exp(0+6))
    halpha = h^alpha
    pow = -expHE/alpha
    
    if(!returnlogs){
      C = (1+halpha)^pow
    }else{
      C = pow*log(1+halpha)
    }
    if(return_dCdtheta){
      if(!returnlogs){
        dCdtheta = cbind(3*expHE*((1+halpha)^(pow-1))*(halpha),3*C*pow*log(1+halpha))
      }else{
        dCdtheta = cbind(3*expHE*halpha/(1+halpha),3*C)
      }
      dCdtheta[is.na(dCdtheta)] = 0
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }
}



#' Gaussian correlation function
#' 
#' Calculate correlation matrix for two sets of points in one dimension
#' Note that this is not the correlation between two vectors.
#' 
#' WE HIGHLY ADVISE NOT USING THIS CORRELATION FUNCTION.
#' Try Power Exponential, CauchySQT, Cauchy, or Matern 3/2 instead.
#'
#' @inheritParams CGGP_internal_CorrMatCauchy
#'
#' @return Matrix of correlation values between x1 and x2
#' @export
#' @family correlation functions
#'
#' @examples
#' CGGP_internal_CorrMatGaussian(c(0,.2,.4),c(.1,.3,.5), theta=c(-.7))
CGGP_internal_CorrMatGaussian <- function(x1, x2,theta, return_dCdtheta = FALSE,
                                          return_numpara=FALSE,
                                          returnlogs=FALSE) {
  if(return_numpara){
    return(1);
  }else{ 
    if (length(theta) != 1) {stop("CorrMatGaussian theta should be length 1")}
    diffmat =abs(outer(x1,x2,'-'))
    diffmat2 <- diffmat^2
    
    expLS = exp(3*theta[1])
    h = diffmat2/expLS
    
    # Gaussian corr is awful, always needs a nugget
    nug <- 1e-10
    if (!returnlogs) {
      C = (1-nug)*exp(-h) + nug*(diffmat<10^(-4))
      # C = exp(-h)
    } else {
      # C = -h
      C = (1-nug)*exp(-h) + nug*(diffmat<10^(-4))
      C <- log(C)
    }
    if(return_dCdtheta){
      if (!returnlogs) {
        dCdtheta <- 3*C*diffmat2 / expLS
      } else {
        dCdtheta <- 3*diffmat2 / expLS
      }
      dCdtheta[is.na(dCdtheta)] = 0
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }
}




#' Matern 3/2 correlation function
#' 
#' Calculate correlation matrix for two sets of points in one dimension.
#' Note that this is not the correlation between two vectors.
#'
#' @inheritParams CGGP_internal_CorrMatCauchy
#'
#' @return Matrix of correlation values between x1 and x2
#' @export
#' @family correlation functions
#'
#' @examples
#' CGGP_internal_CorrMatMatern32(c(0,.2,.4),c(.1,.3,.5), theta=c(-.7))
CGGP_internal_CorrMatMatern32 <- function(x1, x2,theta, return_dCdtheta=FALSE,
                                          return_numpara=FALSE,
                                          returnlogs=FALSE) {
  if(return_numpara){
    return(1);
  }else{ 
    if (length(theta) != 1) {stop("CorrMatMatern32 theta should be length 1")}
    diffmat =abs(outer(x1,x2,'-'))
    
    expLS = exp(3*theta[1])
    h = diffmat/expLS
    
    if (!returnlogs) {
      # C = (1-10^(-10))*(1+sqrt(3)*h)*exp(-sqrt(3)*h) + 10^(-10)*(diffmat<10^(-4))
      C = (1+sqrt(3)*h)*exp(-sqrt(3)*h)
    } else {
      C <- log(1+sqrt(3)*h) - sqrt(3)*h
    }
    if(return_dCdtheta){
      if (!returnlogs) {
        dCdtheta <- (sqrt(3)*diffmat*exp(-sqrt(3)*h) - sqrt(3)*C*diffmat) * (-3/expLS)
      } else {
        dCdtheta <- (sqrt(3)*diffmat/(1+sqrt(3)*h) - sqrt(3)*diffmat) * (-3/expLS)
      }
      dCdtheta[is.na(dCdtheta)] = 0
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }
}



#' Matern 5/2 correlation function
#' 
#' Calculate correlation matrix for two sets of points in one dimension.
#' Note that this is not the correlation between two vectors.
#'
#' @inheritParams CGGP_internal_CorrMatCauchy
#'
#' @return Matrix of correlation values between x1 and x2
#' @export
#' @family correlation functions
#'
#' @examples
#' CGGP_internal_CorrMatMatern52(c(0,.2,.4),c(.1,.3,.5), theta=c(-.7))
CGGP_internal_CorrMatMatern52 <- function(x1, x2,theta, return_dCdtheta=FALSE,
                                          return_numpara=FALSE,
                                          returnlogs=FALSE) {
  if(return_numpara){
    return(1);
  }else{ 
    if (length(theta) != 1) {stop("CorrMatMatern52 theta should be length 1")}
    diffmat =abs(outer(x1,x2,'-'))
    expLS = exp(3*theta[1])
    h = diffmat/expLS
    if (!returnlogs) {
      # C = (1-10^(-10))*(1+sqrt(5)*h+5/3*h^2)*exp(-sqrt(5)*h) + 10^(-10)*(diffmat<10^(-4))
      C = (1+sqrt(5)*h+5/3*h^2)*exp(-sqrt(5)*h)
    } else {
      C = log(1+sqrt(5)*h+5/3*h^2) - sqrt(5)*h
    }
    if(return_dCdtheta){
      if (!returnlogs) {
        dCdtheta <- ((sqrt(5)*diffmat+10/3*diffmat*h)*exp(-sqrt(5)*h) -
                       sqrt(5)*C*diffmat) * (-3/expLS)
      } else {
        dCdtheta <- ((sqrt(5)*diffmat+10/3*diffmat*h)/(1+sqrt(5)*h+5/3*h^2) -
                       sqrt(5)*diffmat) * (-3/expLS)
      }
      dCdtheta[is.na(dCdtheta)] = 0
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }
}



#' Power exponential correlation function
#' 
#' Calculate correlation matrix for two sets of points in one dimension.
#' Note that this is not the correlation between two vectors.
#'
#' @inheritParams CGGP_internal_CorrMatCauchy
#'
#' @return Matrix of correlation values between x1 and x2
# @rdname CGGP_internal_CorrMatCauchy
#' @export
#' @family correlation functions
#'
#' @examples
#' CGGP_internal_CorrMatPowerExp(c(0,.2,.4),c(.1,.3,.5), theta=c(-.7,.2))
CGGP_internal_CorrMatPowerExp <- function(x1, x2,theta,
                                          return_dCdtheta = FALSE,
                                          return_numpara=FALSE,
                                          returnlogs=FALSE) {
  if(return_numpara){
    return(2);
  }else{ 
    if (length(theta) != 2) {stop("CorrMatPowerExp theta should be length 2")}
    diffmat =abs(outer(x1,x2,'-'))
    tmax <- 3
    expLS = exp(tmax*theta[1])
    minpower <- 1
    maxpower <- 1.95
    alpha <- minpower + (theta[2]+1)/2 * (maxpower - minpower)
    h = diffmat/expLS
    if (!returnlogs) {
      # C = (1-nug)*exp(-(h)^alpha) + nug*(diffmat<10^(-4))
      C = exp(-(h)^alpha)
    } else {
      C = -(h^alpha)
    }
    if(return_dCdtheta){
      if (!returnlogs) {
        dCdtheta <- cbind(tmax*alpha*C*diffmat^alpha/expLS^alpha,
                          -C*h^alpha*log(h)/2 * (maxpower - minpower))
      } else {
        dCdtheta <- cbind(tmax*alpha*diffmat^alpha/expLS^alpha,
                          -h^alpha*log(h)/2 * (maxpower - minpower))
      }
      dCdtheta[is.na(dCdtheta)] = 0
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }
}



#' Wendland0 (Triangle) correlation function
#' 
#' Calculate correlation matrix for two sets of points in one dimension.
#' Note that this is not the correlation between two vectors.
#'
#' @inheritParams CGGP_internal_CorrMatCauchy
#'
#' @return Matrix of correlation values between x1 and x2
# @rdname CGGP_internal_CorrMatCauchy
#' @export
#' @family correlation functions
#'
#' @examples
#' CGGP_internal_CorrMatWendland0(c(0,.2,.4),c(.1,.3,.5), theta=-.7)
CGGP_internal_CorrMatWendland0 <- function(x1, x2,theta,
                                           return_dCdtheta = FALSE,
                                           return_numpara=FALSE,
                                           returnlogs=FALSE) {
  if(return_numpara){
    return(1)
  }else{ 
    if (length(theta) != 1) {stop("CorrMatWendland0 theta should be length 1")}
    diffmat =abs(outer(x1,x2,'-'))
    tmax <- 3
    expLS = exp(tmax*theta[1])
    h = diffmat/expLS
    if (!returnlogs) {
      C <- pmax(1 - h, 0)
    } else {
      C = log(pmax(1 - h, 0))
    }
    if(return_dCdtheta){
      if (!returnlogs) {
        dCdtheta <- ifelse(1-h > 0, tmax*diffmat/expLS, 0)
      } else {
        dCdtheta <- ifelse(1-h > 0, tmax*diffmat/expLS/pmax(1-h,0), 0)
      }
      dCdtheta[is.na(dCdtheta)] = 0
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }
}


#' Wendland1 1 correlation function
#' 
#' Calculate correlation matrix for two sets of points in one dimension.
#' Note that this is not the correlation between two vectors.
#'
#' @inheritParams CGGP_internal_CorrMatCauchy
#'
#' @return Matrix of correlation values between x1 and x2
# @rdname CGGP_internal_CorrMatCauchy
#' @export
#' @family correlation functions
#'
#' @examples
#' CGGP_internal_CorrMatWendland1(c(0,.2,.4),c(.1,.3,.5), theta=-.7)
CGGP_internal_CorrMatWendland1 <- function(x1, x2,theta,
                                           return_dCdtheta = FALSE,
                                           return_numpara=FALSE,
                                           returnlogs=FALSE) {
  if(return_numpara){
    return(1)
  }else{ 
    if (length(theta) != 1) {stop("CorrMatWendland1 theta should be length 1")}
    diffmat =abs(outer(x1,x2,'-'))
    tmax <- 3
    expLS = exp(tmax*theta[1])
    h = diffmat/expLS
    if (!returnlogs) {
      C <- pmax(1 - h, 0)^3 * (3*h + 1)
    } else {
      C = log(pmax(1 - h, 0)^3 * (3*h + 1))
    }
    if(return_dCdtheta){
      if (!returnlogs) {
        dCdtheta <- ifelse(1-h > 0,
                           tmax*expLS * (3*(1-h)^2*(h/expLS)*(3*h+1) + (1-h)^3*(-3*h/expLS)),
                           0)
      } else {
        dCdtheta <- ifelse(1-h > 0,
                           tmax*expLS * (3*(1-h)^2*(h/expLS)*(3*h+1) + (1-h)^3*(-3*h/expLS)) /
                             (pmax(1 - h, 0)^3 * (3*h + 1)),
                           0)
      }
      dCdtheta[is.na(dCdtheta)] = 0
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }
}



#' Wendland2 2 correlation function
#' 
#' Calculate correlation matrix for two sets of points in one dimension.
#' Note that this is not the correlation between two vectors.
#'
#' @inheritParams CGGP_internal_CorrMatCauchy
#'
#' @return Matrix of correlation values between x1 and x2
# @rdname CGGP_internal_CorrMatCauchy
#' @export
#' @family correlation functions
#'
#' @examples
#' CGGP_internal_CorrMatWendland2(c(0,.2,.4),c(.1,.3,.5), theta=-.7)
CGGP_internal_CorrMatWendland2 <- function(x1, x2,theta,
                                           return_dCdtheta = FALSE,
                                           return_numpara=FALSE,
                                           returnlogs=FALSE) {
  if(return_numpara){
    return(1)
  }else{ 
    if (length(theta) != 1) {stop("CorrMatWendland2 theta should be length 1")}
    diffmat =abs(outer(x1,x2,'-'))
    tmax <- 3
    expLS = exp(tmax*theta[1])
    h = diffmat/expLS
    if (!returnlogs) {
      C <- pmax(1 - h, 0)^5 * (8*h^2 + 5*h + 1)
    } else {
      C = log(pmax(1 - h, 0)^5 * (8*h^2 + 5*h + 1))
    }
    if(return_dCdtheta){
      if (!returnlogs) {
        dCdtheta <- ifelse(1-h > 0,
                           tmax*expLS * (5*(1-h)^4*(h/expLS)*(8*h^2 + 5*h + 1) + (1-h)^5*(-16*h^2/expLS-5*h/expLS)),
                           # tmax*expLS * (-h/expLS) * (-5*(1-h)^4*(8*h^2 + 5*h + 1) + (1-h)^5*(16*h+5)),
                           0)
      } else {
        dCdtheta <- ifelse(1-h > 0,
                           tmax*expLS * (5*(1-h)^4*(h/expLS)*(8*h^2 + 5*h + 1) + (1-h)^5*(-16*h^2/expLS-5*h/expLS)) /
                             (pmax(1 - h, 0)^5 * (8*h^2 + 5*h + 1)),
                           0)
      }
      dCdtheta[is.na(dCdtheta)] = 0
      out <- list(C=C,dCdtheta=dCdtheta)
      return(out)
    }else{
      return(C)
    }
  }
}
