evaluate_cond_lik <- function(v, theta_x, theta_y, alpha_x, alpha_y, Nx, Ny) {
  nx <- v
  ny <- 1 - nx
  J <- Nx + Ny
  I_X <- alpha_x * nx * (J - 1) / (1 - alpha_x * nx - alpha_y * ny)
  I_Y <- alpha_y * ny * (J - 1) / (1 - alpha_x * nx - alpha_y * ny)

  a <- lgamma(J + 1)
  b <- rep(0, length(I_X))
  poch_X <- rep(0, length(I_X))
  poch_Y <- rep(0, length(I_X))

  for (cnt in seq_along(I_X)) {
    b[cnt] <- lgamma(I_X[cnt] + I_Y[cnt] + J) - lgamma(I_X[cnt] + I_Y[cnt])
    poch_X[cnt] <- lgamma(I_X[cnt] + Nx) - lgamma(I_X[cnt])
    poch_Y[cnt] <- lgamma(I_Y[cnt] + Ny) - lgamma(I_Y[cnt])
  }

  c <- a - b

  h <- poch_X + poch_Y - (lgamma(Nx + 1) + lgamma(Ny + 1))

  k  <-  lgamma( (theta_x / 2) + (theta_y / 2)) -
        (lgamma(theta_x / 2) + lgamma(theta_y / 2))
  l  <- ((theta_x / 2) - 1) * log(nx) + ((theta_y / 2) - 1) * log(ny)

  result <- c + h + k + l
  return(result)
}

calc_conditional <- function(v, model, Nx, Ny) {
  incorrectlength <- 0
  if (model == "D0" && length(v) != 2) incorrectlength <- 1
  if (model == "D1" && length(v) != 3) incorrectlength <- 1

  if (incorrectlength == 1) {
    stop("calcConditional:",
         "Input vector is of incorrect length\n")
  }

  theta_x <- v[1] * 2
  theta_y <- v[1] * 2
  alpha_x <- v[2]
  alpha_y <- v[2]

  if (model == "D1") {
    alpha_y <- v[3]
  }

  if (alpha_x < 0 ||
      alpha_y < 0 ||
      theta_x < 1 ||
      theta_y < 1 ||
      alpha_x > (1 - (1e-8)) ||
      alpha_y > (1 - (1e-8))
     ) {
    return(-Inf)
  }

  f <- function(x) {
    return( -1 * evaluate_cond_lik(x, theta_x, theta_y,
                                 alpha_x, alpha_y,
                                 Nx, Ny) )
  }

 maxes <- pracma::fminbnd(f, 0, 1, maxiter = 500, tol = 1e-4 )

 ymax <- -1 * maxes$fmin
 xmax <- maxes$xmin
 xlft <- 0
 xrgt <- 1
 eps <- .Machine$double.eps
 thrs <- 10

 check_left_x  <- evaluate_cond_lik(eps, theta_x, theta_y,
                               alpha_x, alpha_y, Nx, Ny)  - ymax + thrs
 check_right_x <- evaluate_cond_lik(1 - eps, theta_x, theta_y,
                               alpha_x, alpha_y, Nx, Ny)  - ymax + thrs

 if (check_left_x < 0) {
     g <- function(x) {
        return(evaluate_cond_lik(x, theta_x, theta_y,
                               alpha_x, alpha_y, Nx, Ny) - ymax + thrs)
     }
     xlft <- (uniroot(g, c(eps, xmax)))$root
 }

 if (check_right_x < 0) {
     h <- function(x) {
        return(evaluate_cond_lik(x, theta_x, theta_y,
                               alpha_x, alpha_y,
                               Nx, Ny) - ymax + thrs)
     }
     xrgt <- (uniroot( h, c( xmax, 1 - eps)))$root
 }

 calc_ll_exp <- function(x) {
   out <- exp(evaluate_cond_lik(x, theta_x, theta_y,
                              alpha_x, alpha_y,
                              Nx, Ny) - ymax)
   return(out)
 }

  aux <- integrate(f = calc_ll_exp,
                   lower = xlft,
                   upper = xrgt,
                   abs.tol = 1e-9)

  LL <- log(aux$value) + ymax

  return(LL)
}


logLikelihood.Guilds.Conditional <- function(parameters, model,
                                             sadx, sady,
                                             verbose = TRUE) {
  Nx <- sum(sadx)
  Ny <- sum(sady)

  if (verbose == TRUE) {
   cat("Chosen model: ", model, "\n")
   cat("Now starting to calculate likelihood of: \n")
   x2 <- parameters
   if (model == "D0") cat("Theta X =", x2[1],
                         " Theta Y =", "Theta X",
                         "\t Alpha X =", x2[2],
                         " Alpha Y =", "Alpha X", "\n")
   if (model == "D1") cat("Theta X =", x2[1],
                         " Theta Y =", "Theta X",
                         "\t Alpha X =", x2[2],
                         " Alpha Y =", x2[3], "\n")

   flush.console()
  }

  ll <- logLikelihood.Guilds(parameters, model, sadx, sady, verbose)

  #conditional part:
  conditional_part <- calc_conditional(parameters, model, Nx, Ny)

  output <- ll - conditional_part
  return(output)
}

conditional.LogLik <- function(v, model, J, Sx, Sy, Nx, Ny, kda_x, kda_y,
                               prefactor1, prefactor2, verbose = TRUE) {

  theta_x <- v[1] * 2
  theta_y <- v[1] * 2
  alpha_x <- v[2]
  alpha_y <- v[2]

  if (model == "D1") {
    alpha_y <- v[3]
  }

  if (alpha_x < 0 ||
     alpha_y < 0 ||
     theta_x < 1 ||
     theta_y < 1 ||
     alpha_x > (1 - (1e-8)) ||
     alpha_y > (1 - (1e-8))
     ) return(-Inf)

  ll <- logLikguilds(theta_x, theta_y, alpha_x, alpha_y, J,
                     Sx, Sy, Nx, Ny, kda_x, kda_y,
                     prefactor1, prefactor2, verbose)

  cond_ll <- calc_conditional(v, model, Nx, Ny)
  out <- ll - cond_ll
  return(out)
}

maxLikelihood.Guilds.Conditional <- function(init_vals, model,
                                             method,
                                             sadx, sady,
                                             verbose = TRUE) {
  incorrectlength <- 0
  if (model == "D0" && length(init_vals) != 2) incorrectlength <- 1
  if (model == "D1" && length(init_vals) != 3) incorrectlength <- 1
  
  if (incorrectlength == 1) { 
    stop("maxLikelihood.Guilds.Conditional: ",
         "Input vector is of incorrect length\n")
  }
  
  if(init_vals[1] < 1) {
    stop("maxLikelihood.Guilds.Conditional: ",
         "initial theta can not be below one")
  }
  if(init_vals[2] < 0) {
    stop("maxLikelihood.Guilds.Conditional: ",
         "initial alpha can not be below zero")
  }
  if(init_vals[2] > 1) {
    stop("maxLikelihood.Guilds.Conditional: ",
         "initial alpha can not be above 1")
  }
  if(model == "D1") {
    if(init_vals[3] < 0 ) {
      stop("maxLikelihood.Guilds.Conditional: ",
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
  for (i in seq_along(x)) freq_x[i] <- x[[i]]

  prefactor1 <- -1 * ( sum(log(sadx)) + sum(lgamma(1 + freq_x)) )

  x2 <- c(table(sady))
  freq_y <- c()
  for (i in seq_along(x2)) freq_y[i] <- x2[[i]]

  prefactor2 <- -1 * ( sum(log(sady)) + sum(lgamma(1 + freq_y)) )

  Sx <- length(sadx)
  Sy <- length(sady)
  Nx <- sum(sadx)
  Ny <- sum(sady)
  J  <- Nx + Ny

  g <- function(x) {
    out <- -1 * conditional.LogLik(x, model, J, Sx, Sy, Nx, Ny,
                                   kda_x, kda_y, prefactor1,
                                   prefactor2, verbose)
    return(out)
  }

  x
  if (method == "simplex") {
    x <- simplex(init_vals, g, verbose)
  }
#  if (method == "subplex") {
#    x <- subplex::subplex(init_vals, g)
#  }
  return(x)
}