generate.Guilds <- function(theta, alpha_x, alpha_y, J) {
  if(theta < 1) {
    stop("generate.Guilds: ",
         "theta can not be below one")
  }
  if(alpha_x < 0) {
    stop("generate.ESF: ",
         "alpha_x can not be below zero")
  }
  if(alpha_y < 0) {
    stop("generate.ESF: ",
         "alpha_y can not be below zero")
  }
  if(alpha_x > 1) {
    stop("generate.ESF: ",
         "alpha_x can not be above one")
  }
  if(alpha_y > 1) {
    stop("generate.ESF: ",
         "alpha_y can not be above one")
  }
  
  if(J < 0) {
    stop("generate.ESF: ",
         "J can not be below zero")
  }

  #first draw nx and ny from a beta distribution
  nx <- rbeta(1, theta, theta)
  ny <- 1 - nx

  #update I_X and I_Y accordingly
  I_X <- alpha_x * nx * (J-1) / (1 - alpha_x * nx - alpha_y * ny) 
  I_Y <- alpha_y * ny * (J-1) / (1 - alpha_x * nx - alpha_y * ny)

  probs <- c()
  allN <- 0:J
  if(is.infinite(I_X) && is.infinite(I_Y)) {
     probs <- exp( lgamma(J+1) - 
                 (lgamma(allN + 1) + 
                  lgamma(J - allN + 1)) + 
                  allN * log(nx) + 
                 (J-allN) * log(ny))
  } else {
    #set up a probability vector
    probs <- polyaeggenberger(I_X, I_Y, J, allN) 
  }

  NX <- sample(0:J, 1, replace = TRUE, prob = probs)
  NY <- J - NX

  sadx <- generate.ZSM(theta, I_X, NX)
  sady <- generate.ZSM(theta, I_Y, NY)

  output <- list( guildX = sadx, guildY = sady)
  
  return(output)
}

polyaeggenberger <- function(theta_x,theta_y,J,N) {
  a <- lgamma(J + 1) #J!
  b <- lgamma(theta_x + theta_y + J) - 
      lgamma(theta_x + theta_y) #(theta_x + theta_y)_J pochhammer
  
  c1 <- lgamma(theta_x + N) - lgamma(theta_x) #(theta_x)_Nx  pochhammer
  c2 <- lgamma(theta_y + J - N) - lgamma(theta_y) #(theta_y)_Ny  pochhammer
  d <- lgamma(N + 1) + lgamma(J-N+1)
  
  return( exp(a - b + c1 + c2 - d))
}

localComm <- function(alpha_x, alpha_y, Jx, Jy, px) {
  J <- Jx + Jy

  nx <- px
  ny <- 1 - nx

  I_X <- alpha_x * nx * (J-1) / (1 - alpha_x * nx - alpha_y * ny) 
  I_Y <- alpha_y * ny * (J-1) / (1 - alpha_x * nx - alpha_y * ny)

 if(is.infinite(I_X) && is.infinite(I_Y)) {
     output <- exp( lgamma(J+1) - 
                  (lgamma(Jx + 1) + 
                   lgamma(J - Jx + 1)) +
                    Jx * log(nx) + 
                    (J-Jx) * log(ny))
 } else {
  output <- polyaeggenberger(I_X, I_Y, J, Jx)
 }
  return(output)
}

rho <- function(theta,px) {
  a <- lgamma(2 * theta) - 2 * lgamma(theta)
  b <- (theta - 1) * log(px) + (theta - 1) * log((1-px))
  output <- exp(a + b)
  return(output)
}

getpx <- function(theta, alpha_x, alpha_y, JX, JY) {
  px <- (1:(1000-1)) / 1000
  calcLocal <- function(x) {
    a <- localComm(alpha_x, alpha_y, JX, JY, x) * 
               rho(theta,x)
    return(a)
  }
  divider <- integrate(f = calcLocal,
                       lower = 0,
                       upper = 1,
                       abs.tol = 1e-9)$value
  
  probs <- calcLocal(px) / divider

  return(sample(px, 1, prob = probs))
}

generate.Guilds.Cond <- function(theta, alpha_x, alpha_y, JX, JY) {

  J <- JX + JY

  nx <- getpx(theta, alpha_x, alpha_y, JX, JY)
  ny <- 1 - nx

  I_X <- alpha_x * nx * (J - 1) / (1 - alpha_x * nx - alpha_y * ny)
  I_Y <- alpha_y * ny * (J - 1) / (1 - alpha_x * nx - alpha_y * ny)

  sadx <- generate.ZSM(theta, I_X, JX)
  sady <- generate.ZSM(theta, I_Y, JY)

  output <- list( guildX = sadx, guildY = sady)
  
  return(output)
}


generate.ZSM <- function(theta,I,J) {

  #based on urn.gp code by Rampal Etienne available as a supplementary
  #online appendix for his 2005 paper in Ecology Letters
  if ( is.infinite(I) ) {
    I <- .Machine$double.xmax
  }

  anc_ct <- 0 
  spec_ct <- 0
  ancestor <- rep(NaN, J)
  species <- rep(NaN, J)
  ancestor_species <- rep(NaN, J)

  for (j in 1:J) { #loop over each invidual 
    #ancestor not previously encountered
     if ( runif(1,0,1) <= I / (I + j - 1) ) { 
        anc_ct <- anc_ct + 1
        ancestor[j] <- anc_ct
        #new ancestor is a new species
        if ( runif(1,0,1) <= theta / (theta + anc_ct - 1)) { 
            spec_ct <- spec_ct + 1
            species[j] <- spec_ct
            ancestor_species[anc_ct] <- spec_ct
         } else { 
           # new ancestor is an existing species but
           # new migration into local community
            prior_lineage <- pracma::ceil( runif(1, 0, 1) * (anc_ct - 1))
            s <- ancestor_species[prior_lineage]
            species[j] <- s
            ancestor_species[anc_ct] <- s
         }
     } else { #descendant of existing individual
       #select one individual at random
        descend_from <- pracma::ceil( runif(1, 0, 1) * (j - 1)) 
        ancestor[j] <- ancestor[descend_from]
        species[j] <- species[descend_from]
     }
  }
  
  x <- c( table( species))
  abund <- c()
  for (i in seq_along(x)) {
     abund[i] <- x[[i]]
  }
  abund <- sort( abund, decreasing = TRUE)
  return(abund)
}

generate.ESF <- function(theta,I,J) {
  if(theta < 1) {
    stop("generate.ESF: ",
         "theta can not be below one")
  }
  if(I < 0) {
    stop("generate.ESF: ",
         "I can not be below zero")
  }
  
  if(J < 0) {
    stop("generate.ESF: ",
         "J can not be below zero")
  }
  
  return( generate.ZSM( theta, I, J) )
}