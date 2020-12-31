##' Epidemiology problem, initial and rescaled to [0,1]^2 versions.
##' @title SIR test problem
##' @param x vector of size two
##' @rdname SIR
##' @references 
##' R. Hu, M. Ludkovski (2017), Sequential Design for Ranking Response Surfaces, SIAM/ASA Journal on Uncertainty Quantification, 5(1), 212-239.
##' @export
##' @examples 
##' ## SIR test problem illustration
##' ngrid <- 10 # increase
##' xgrid <- seq(0, 1, length.out = ngrid)
##' Xgrid <- as.matrix(expand.grid(xgrid, xgrid))
##' 
##' nrep <- 5 # increase
##' X <- Xgrid[rep(1:nrow(Xgrid), nrep),]
##' Y <- apply(X, 1, sirEval)
##' dataSIR <- find_reps(X, Y)
##' filled.contour(xgrid, xgrid, matrix(lapply(dataSIR$Zlist, sd), ngrid),
##'                xlab = "Susceptibles", ylab = "Infecteds", color.palette = terrain.colors)
##'
sirEval <- function(x){
  return(sirSimulate(S0 = x[1] * 600 + 1200, I0 = 200 * x[2], M = 2000, beta = 0.5, imm = 0)$totI/800)
}


##' @param S0 initial nunber of susceptibles
##' @param I0 initial number of infected
##' @param M total population
##' @param beta,gamma,imm control rates
##' @rdname SIR
##' @importFrom stats runif
##' @export
sirSimulate <- function(S0 = 1990, I0 = 10, M = S0 + I0, beta = 0.75, gamma = 0.5, imm = 0)
{
  curS <- rep(S0,M); curI <- rep(I0,M); curT <- rep(0,M)
  curS[1] <- S0
  curI[1] <- I0
  curT[1] <- 0
  count <- 1
  maxI <- I0
  
  # continue until no more infecteds
  while ( (curI[count] >0 & imm==0) | (curT[count] < 100 & imm > 0) ) {
    
    # gillespie SSA algorithm: 2 reactions possible
    infRate <- beta*curS[count]/(M)*(curI[count]+imm)
    recRate <- gamma*curI[count]
    infTime <- -1/infRate*log( runif(1))
    recTime <- -1/recRate*log( runif(1))
    
    if (infTime < recTime) {   # infection
      curS[count+1] <- curS[count] - 1
      curI[count+1] <- curI[count] + 1
      maxI <- max(maxI, curI[count+1])
    }     else {    # recovery
      curI[count+1] <- curI[count] - 1
      curS[count+1] <- curS[count]
    }
    curT[count+1] <- curT[count] + min( infTime,recTime)
    count <- count+1
  }
  return(list(maxI = maxI, totT = curT[count], totI = S0-curS[count],S=curS,I=curI,R=M-curS-curI,T=curT))
}

#' 1d test function (1)
#' @param x scalar or matrix (size n x 1) in [0,1]
#' @export
#' @references 
#' A. Forrester, A. Sobester, A. Keane (2008), Engineering design via surrogate modelling: a practical guide, John Wiley & Sons
#' @examples 
#' plot(f1d)
f1d <- function(x){
  if(is.null(dim(x))) x <- matrix(x, ncol = 1)
  return(((x*6-2)^2)*sin((x*6-2)*2))
}

#' Noisy 1d test function (1)
#' Add Gaussian noise with variance r(x) = scale * (1.1 + sin(2 pi x))^2 to \code{\link[hetGP]{f1d}}
#' @param x scalar or matrix (size n x 1) in [0,1]
#' @param scale scalar in [0, Inf] to control the signal to noise ratio
#' @export
#' @examples 
#' X <- matrix(seq(0, 1, length.out = 101), ncol = 1)
#' Xr <- X[sort(sample(x = 1:101, size = 500, replace = TRUE)),, drop = FALSE]
#' plot(Xr, f1d_n(Xr))
#' lines(X, f1d(X), col = "red", lwd = 2)
f1d_n <- function(x, scale = 1){
  if(is.null(dim(x))) x <- matrix(x, ncol = 1)
  return(rnorm(n = nrow(x), mean = f1d(x), sd = scale *  (1.1 + sin(2 *pi* x))^2))
}

#' 1d test function (2)
#' @param x scalar or matrix (size n x 1) in [0,1]
#' @export
#' @references 
#' A. Boukouvalas, and D. Cornford (2009), Learning heteroscedastic Gaussian processes for complex datasets, Technical report. \cr \cr
#' M. Yuan, and  G. Wahba (2004), Doubly penalized likelihood estimator in heteroscedastic regression, Statistics and Probability Letters 69, 11-20.
#' @examples 
#' plot(f1d2)
f1d2 <- function(x){
  if(is.null(dim(x))) x <- matrix(x, ncol = 1)
  return(2 * (exp(-30*(x-0.25)^2) + sin(pi * x^2)) - 2)
}

#' Noisy 1d test function (2)
#' Add Gaussian noise with variance r(x) = scale * (exp(sin(2 pi x)))^2 to \code{\link[hetGP]{f1d2}}
#' @param x scalar or matrix (size n x 1) in [0,1]
#' @param scale scalar in [0, Inf] to control the signal to noise ratio
#' @export
#' @examples 
#' X <- matrix(seq(0, 1, length.out = 101), ncol = 1)
#' Xr <- X[sort(sample(x = 1:101, size = 500, replace = TRUE)),, drop = FALSE]
#' plot(Xr, f1d2_n(Xr))
#' lines(X, f1d2(X), col = "red", lwd = 2)
f1d2_n <- function(x, scale = 1){
  if(is.null(dim(x))) x <- matrix(x, ncol = 1)
  return(rnorm(n = nrow(x), mean = f1d2(x), sd = scale *  (exp(sin(2*pi*x)))))
}

##' Portfolio value at risk test problem
##' @title Portfolio simulation
##' @param x two dimensional vector of inputs in [?, ?]^2
##' @param S1_0,S2_0,K1,K2,T1,T2,sigma1,sigma2,r,rho parameters
## ' @references 
##' @noRd
##' @examples 
##' # outer scenarios
##' nscen <- 1000; T0 <- 1; nvar <- 2
##' Xscen <- array(0, dim=c(nscen, 2))
##' z1 <- rnorm(nscen)
##' S1_0 <- 50; S2_0 <- 80; r <- 0.04
##' sigma1 <- 0.25; sigma2 <- 0.35; rho <- 0.3
##' Xscen[,1] <- S1_0 * exp((r-sigma1^2/2)*T0 + sigma1*sqrt(T0)*z1)
##' Xscen[,2] <- S2_0 * exp((r-sigma2^2/2)*T0 + sigma2*sqrt(T0)*(rho*z1 + sqrt(1-rho^2)*rnorm( nscen)))
##' 
##' ## Unique (randomly chosen) design locations
##' n <- 100
##' Xu <- Xscen[ sample(1:nscen, n),]
##' X <- Xu[sample(1:n, 100*n, replace = TRUE),]
##' 
##' ## obtain training data response at design locations X
##' Z <- bsVar(X)
##' 
##' ## Formating of data for model creation (find replicated observations) 
##' prdata <- find_reps(X, Z)
##' 
##' ## Model fitting
##' model <- mleHetGP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult), Z = prdata$Z,
##'                   lower = rep(0.1, nvar), upper = rep(200, nvar),
##'                                     covtype = "Matern5_2")
##'
##'## a quick view into the data stored in the "hetGP"-class object
##'summary(model)                  
##'
##'## prediction from the fit on the grid  
##'ngrid <- 51
##'xgrid <- matrix(seq(0, 1, length.out = ngrid), ncol = 1) 
##'Xgrid <- as.matrix(expand.grid(30+100*xgrid, 30+170*xgrid))
##'
##'predictions <- predict(x = Xgrid, object =  model)
##'
##'## Visualization of the predictive surface
##'par(mfrow = c(2, 1))
##'contour(x = 30+100*xgrid,  y = 30+170*xgrid, z = matrix(predictions$mean, ngrid), 
##'        main = "Predicted mean", nlevels = 20)
##'        points(X, col = 'blue', pch = 20)
##'        contour(x = 30+100*xgrid,  y= 30+170*xgrid, z = matrix(sqrt(predictions$nugs), ngrid), 
##'                main = "Predicted noise values", nlevels = 20)
##'                points(X, col = 'blue', pch = 20)
##'                par(mfrow = c(1, 1))
##'                
bsVar <- function(x, S1_0 = 50, S2_0 = 80, K1 = 40, K2 = 85, T1 = 1, T2 = 2, sigma1 = 0.25, sigma2 = 0.35, r = 0.04, rho = 0.3){
  if(is.null(nrow(x)))
    x <- matrix(x, nrow = 1)
  z1 <- rnorm(nrow(x))
  logs1 <- (r - sigma1^2/2)*T1 + sigma1*sqrt(T1)*z1
  logs2 <- (r - sigma2^2/2)*T2 + sigma2*sqrt(T2)*(rho*z1 + sqrt(1-rho^2)*rnorm(nrow(x)))
  portf <- exp(-r*T1)*100*pmax( x[,1]*exp(logs1) - K1,0)  - 
    exp(-r*T2)*50*pmax( x[,2]*exp(logs2) - K2, 0)
  return(portf)
}

