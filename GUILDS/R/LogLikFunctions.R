logLikguilds <- function(theta_x, theta_y,
                         alpha_x, alpha_y,
                         J, Sx, Sy, Nx, Ny,
                         KDA_X, KDA_Y,
                         prefactor1, prefactor2,
                         verbose=TRUE) {
 thrs <- 10
 
 f <- function(x) {
   return( -1 * evaluateLogLik(x, theta_x, theta_y, alpha_x, alpha_y,
                               J, Sx, Sy, Nx, Ny, KDA_X, KDA_Y) )
 }

 maxes <- pracma::fminbnd(f, 0, 1,
                          maxiter = 500, tol = 1e-4)

 ymax <- -1 * maxes$fmin
 xmax <- maxes$xmin
 xlft <- 0
 xrgt <- 1
 eps <- .Machine$double.eps

 check_left_x  <- evaluateLogLik(eps, theta_x, theta_y, alpha_x, alpha_y,
                              J, Sx, Sy, Nx, Ny, KDA_X, KDA_Y) - ymax + thrs
 check_right_x <- evaluateLogLik(1 - eps, theta_x, theta_y, alpha_x, alpha_y,
                              J, Sx, Sy, Nx, Ny, KDA_X, KDA_Y) - ymax + thrs

 if(check_left_x < 0) {
     g <- function(x) {
        return(evaluateLogLik(x, theta_x, theta_y, alpha_x, alpha_y,
                              J, Sx, Sy, Nx, Ny, KDA_X, KDA_Y) - ymax + thrs)
     }
     xlft <- (uniroot(g, c(eps, xmax)))$root
 }

 if(check_right_x < 0) {
     h <- function(x) {
        return(evaluateLogLik(x, theta_x, theta_y, alpha_x, alpha_y,
                              J, Sx, Sy, Nx, Ny, KDA_X, KDA_Y) - ymax + thrs)
     }
     xrgt <- (uniroot(h, c(xmax, 1 - eps)))$root
 }

 calc_ll_exp <- function(x) {
   out <- exp( evaluateLogLik( x, theta_x, theta_y, alpha_x, alpha_y,
                               J, Sx, Sy, Nx, Ny, KDA_X, KDA_Y) - ymax)
   return(out)
 }

 aux <- integrate(f = calc_ll_exp, lower = xlft, upper = xrgt, abs.tol=1e-9)

 y <- ymax + log(aux$value) + prefactor1 + Sx * log(theta_x / 2) + 
   prefactor2 + Sy * log(theta_y / 2) + lgamma((theta_x / 2) + (theta_y / 2)) - 
   (lgamma(theta_x/2) + lgamma(theta_y / 2))
 
 #if(verbose==TRUE) 
 # cat(sprintf("%4.3f       %4.4f       %4.4f       %4.4f        %4.3f\n",
 #  theta_x, theta_y, alpha_x, alpha_y, y)); flush.console(); 
 
 return(y)
}


maxLikelihood.Guilds <- function(init_vals, model, 
                                 method = "simplex", sadx, sady, 
                                 verbose = TRUE) {
  incorrectlength <- 0
  if (model == "D0" && length(init_vals) != 2) incorrectlength <- 1
  if (model == "D1" && length(init_vals) != 3) incorrectlength <- 1
  
  if (incorrectlength == 1) { 
    stop("maxLikelihood.Guilds: ",
         "Input vector is of incorrect length\n")
  }
  
  if(init_vals[1] < 1) {
    stop("maxLikelihood.Guilds: ",
         "initial theta can not be below one")
  }
  if(init_vals[2] < 0) {
    stop("maxLikelihood.Guilds: ",
         "initial alpha can not be below zero")
  }
  if(init_vals[2] > 1) {
    stop("maxLikelihood.Guilds: ",
         "initial alpha can not be above 1")
  }
  if(model == "D1") {
    if(init_vals[3] < 0 ) {
      stop("maxLikelihood.Guilds: ",
           "initial alpha_y can not be below 0")
    }
    if(init_vals[3] > 1 ) {
      stop("maxLikelihood.Guilds: ",
           "initial alpha_y can not be above 1")
    }
  }
  
  kda_x <- calcKDA(sadx)
  kda_y <- calcKDA(sady)
    
  x <- c(table(sadx))
  freq_x <- c()
  for (i in seq_along(x) ) freq_x[i] <- x[[i]]

  prefactor1 <- -( sum(log(sadx)) + sum(lgamma(1 + freq_x)) )
  
  x2 <- c(table(sady))
  freq_y <- c()
  for(i in seq_along(x2)) freq_y[i] <- x2[[i]]
  prefactor2 <- -( sum(log(sady)) + sum(lgamma(1 + freq_y)) )

  Sx <- length(sadx)
  Sy <- length(sady)
  Nx <- sum(sadx)
  Ny <- sum(sady)
  J <- Nx + Ny
  
  loglikverbose <- verbose
  if(method == "simplex") loglikverbose <- FALSE

  g <- function(v) {  
    theta_x <- v[1] * 2
    theta_y <- v[1] * 2
    alpha_x <- v[2]
    alpha_y <- v[2]
    
    if(model == "D1") {
      alpha_y <- v[3]
    }

    y <- 0
  	if (alpha_x < 0 || 
  	    alpha_y < 0 || 
  	    theta_x < 1 || 
  	    theta_y < 1  ||
  	    alpha_x > (1 - ( 1e-8)) ||
  	    alpha_y > (1 - ( 1e-8))
  	    ) {
  	  y <- -Inf
  	}
  	if(!is.infinite(y)) {
  	  y <- logLikguilds(theta_x, theta_y, alpha_x, alpha_y,
  	                 J, Sx, Sy, Nx, Ny, kda_x, kda_y,
  	                 prefactor1, prefactor2, verbose)
  	}
  	return(-y)
  }

  x
  if (method == "simplex") {
    	x <- simplex(init_vals,g,verbose)
  }
#  if (method == "subplex") {
#	  x <- subplex::subplex(init_vals,g)
#  }

  return(x)
}

logLikelihood.Guilds <- function(parameters, model, 
                                 sadx, sady, 
                                 verbose = TRUE) {
  kda_x <- calcKDA(sadx)
  kda_y <- calcKDA(sady)
  
  x <- c(table(sadx))
  freq_x <- c()
  for(i in seq_along(x)) freq_x[i] <- x[[i]]
  prefactor1 <- -1 * ( sum(log(sadx)) + sum(lgamma(1 + freq_x)) )
  
  x2 <- c(table(sady))
  freq_y <- c()
  for(i in seq_along(x2)) freq_y[i] <- x2[[i]]
  prefactor2 <- -1 * ( sum(log(sady)) + sum(lgamma(1 + freq_y)) )

  Sx <- length(sadx)
  Sy <- length(sady)
  Nx <- sum(sadx)
  Ny <- sum(sady)
  J <- Nx + Ny
  
  
  if (model == "D0") { 
    #because theta_x = theta_y = theta/2
    theta_x <- parameters[1] * 2
    theta_y <- parameters[1] * 2
    alpha_x <- parameters[2] 
    alpha_y <- parameters[2]
  }
  if(model == "D1") {
    theta_x <- parameters[1] * 2
    theta_y <- parameters[1] * 2
    alpha_x <- parameters[2]
    alpha_y <- parameters[3]
  }

  ll <- logLikguilds(theta_x, theta_y, alpha_x, alpha_y, J,
                     Sx, Sy, Nx, Ny, kda_x, kda_y, prefactor1,
                     prefactor2, verbose)
  
  if ( verbose == TRUE ) {
    cat("Likelihood is ", ll, "\n")
  }
  
  return(ll)
}