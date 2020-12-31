aggregate.data.complete=function(dat,id){
  dat1=dat[order(dat[,id]),]
  nloc=max(dat1[,id])

  ind=which(colnames(dat1)==id)

  res=matrix(NA,nloc,ncol(dat1)-1)
  n=rep(NA,nloc)
  for (i in 1:nloc){
    cond=dat1[,id]==i
    dat.tmp=dat1[cond,-ind]
    n[i]=nrow(dat.tmp)
    if (n[i]>1)  dat.tmp=colSums(dat.tmp)
    if (n[i]==1) dat.tmp=as.numeric(dat.tmp)
    res[i,]=dat.tmp
  }
  colnames(res)=colnames(dat1[,-ind])
  loc.id=unique(dat1[,id])
  list(dat=res,loc.id=id,n=n)
}

#' @name rlda.bernoulli
#' @title Gibbs Sampling for LDA Presence and Absence
#' @description Compute the Gibbs Sampling for LDA Presence and Absence
#' @param data - DataFrame with Presence and Absecence (Zeros and Ones)
#' @param loc.id - Location ID variable
#' @param n_community - Number of communities
#' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
#' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
#' @param gamma - Hyperparameter  Beta(1,gamma)
#' @param n_gibbs - Total number of Gibbs Samples
#' @param ll_prior - Likelihood compute with Priors ?
#' @param bool display_progress=true - Should I Show the progressBar ?
#' @return Rlda object
#' @export
rlda.fastbernoulli <- function(data, loc.id, n_community, alpha0, alpha1, gamma, n_gibbs, ll_prior = TRUE, display_progress = TRUE) {
  #----------------------------------------------------------------
  get.theta=function(nlk,gamma,ncomm,nloc){
    vmat=matrix(NA,nloc,ncomm)
    for (i in 1:(ncomm-1)){
      if (i==(ncomm-1)) cumsoma=nlk[,ncomm]
      if (i< (ncomm-1)) cumsoma=rowSums(nlk[,(i+1):ncomm])
      vmat[,i]=rbeta(nloc,nlk[,i]+1,cumsoma+gamma)
    }
    vmat[,ncomm]=1
    convertVtoTheta(vmat,rep(1,nloc))
  }


  # Create a stop point
  stopifnot(inherits(data, "data.frame"))
  stopifnot(inherits(n_community, "numeric"))
  stopifnot(inherits(alpha0, "numeric"))
  stopifnot(inherits(alpha1, "numeric"))
  stopifnot(inherits(gamma, "numeric") | is.na(gamma))
  stopifnot(inherits(n_gibbs, "numeric"))
  stopifnot(inherits(ll_prior, "logical"))
  stopifnot(inherits(display_progress, "logical"))

  tmp=aggregate.data.complete(data,loc.id)

  y=tmp$dat
  n=tmp$n

  #useful stuff
  ind=which(colnames(data)==loc.id)
  dat1=data.matrix(data[,-ind])
  loc.id=data[,loc.id]
  nloc=max(loc.id)
  nspp=ncol(dat1)
  ncomm=n_community
  nlinhas=nrow(dat1)
  hi=0.999999
  lo=0.000001

  theta=matrix(1/ncomm,nloc,ncomm)
  phi=matrix(0.5,ncomm,nspp)
  z=matrix(sample(1:ncomm,size=nlinhas*nspp,replace=T),nlinhas,nspp)

  #priors
  a.phi=alpha0
  b.phi=alpha1

  #gibbs details
  ngibbs=n_gibbs
  theta.out=matrix(NA,ngibbs,ncomm*nloc)
  phi.out=matrix(NA,ngibbs,ncomm*nspp)
  llk=rep(NA,ngibbs)
  options(warn=2)
  if(display_progress) pb   <- txtProgressBar(1, ngibbs, style=3)
  for (i in 1:ngibbs){
    #sample z
    rand.u=matrix(runif(nlinhas*nspp),nlinhas,nspp)
    z=samplez(log(theta), log(1-phi), log(phi), dat1, loc.id,rand.u, ncomm, nloc)

    #calculate summaries
    tmp=getks(z=z, ncommun=ncomm, dat=dat1)
    nks1=tmp$nks1
    nks0=tmp$nks0
    nlk=getlk(z=z,locid=loc.id, ncommun=ncomm, nloc=nloc)

    #get parameters
    theta=get.theta(nlk,gamma,ncomm,nloc) #theta.true#
    theta[theta>hi]=hi; theta[theta<lo]=lo
    phi=matrix(rbeta(nspp*ncomm,nks1+a.phi,nks0+b.phi),ncomm,nspp) #phi.true#
    phi[phi>hi]=hi; phi[phi<lo]=lo
    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo
    prob1=prob[loc.id,]

    #store results
    llk[i]=sum(dat1*log(prob1)+(1-dat1)*log(1-prob1))
    theta.out[i,]=theta
    phi.out[i,]=phi

    if(display_progress) setTxtProgressBar(pb, i)
  }
  res<-list()
  #Log-Likelihood
  res$logLikelihood<-llk
  #Theta
  res$Theta<-theta.out
  #Phi
  res$Phi<-phi.out
  # Type distribution
  res$type <- "Bernoulli"
  # Number of communities
  res$n_community <- n_community
  # Sample size
  res$N <- nrow(data)
  # Covariates
  res$Species <- colnames(data)
  res$Species <- res$Species [! res$Species %in% loc.id]
  # Alpha0
  res$alpha0 <- alpha0
  # Alpha1
  res$alpha1 <- alpha1
  # Gamma
  res$gamma <- gamma
  # Number of gibbs
  res$n_gibbs <- n_gibbs
  # Species
  res$colnames <- colnames(data)
  # Locations
  res$rownames <- rownames(data)
  # Create the class
  class(res) <- c("rlda", "list")
  return(res)
}

#' @name rlda.bernoulli
#' @title Gibbs Sampling for LDA Presence and Absence
#' @description Compute the Gibbs Sampling for LDA Presence and Absence
#' @param data - DataFrame with Presence and Absecence (Zeros and Ones)
#' @param n_community - Number of communities
#' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
#' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
#' @param gamma - Hyperparameter  Beta(1,gamma)
#' @param n_gibbs - Total number of Gibbs Samples
#' @param ll_prior - Likelihood compute with Priors ?
#' @param bool display_progress=true - Should I Show the progressBar ?
#' @return Rlda object
#' @export
rlda.bernoulli <- function(data, n_community, alpha0, alpha1, gamma, n_gibbs, ll_prior = TRUE, display_progress = TRUE) {
  # Create a stop point
  stopifnot(inherits(data, "data.frame"))
  stopifnot(inherits(n_community, "numeric"))
  stopifnot(inherits(alpha0, "numeric"))
  stopifnot(inherits(alpha1, "numeric"))
  stopifnot(inherits(gamma, "numeric") | is.na(gamma))
  stopifnot(inherits(n_gibbs, "numeric"))
  stopifnot(inherits(ll_prior, "logical"))
  stopifnot(inherits(display_progress, "logical"))

  # Use a function not exported Execute the LDA for the Bernoulli entry
  res <- lda_bernoulli(data, n_community, alpha0, alpha1, gamma, n_gibbs, ll_prior, display_progress)
  # Type distribution
  res$type <- "Bernoulli"
  # Number of communities
  res$n_community <- n_community
  # Sample size
  res$N <- nrow(data)
  # Covariates
  res$Species <- colnames(data)
  # Alpha0
  res$alpha0 <- alpha0
  # Alpha1
  res$alpha1 <- alpha1
  # Gamma
  res$gamma <- gamma
  # Number of gibbs
  res$n_gibbs <- n_gibbs
  # Species
  res$colnames <- colnames(data)
  # Locations
  res$rownames <- rownames(data)
  # Create the class
  class(res) <- c("rlda", "list")
  return(res)
}

#' @name rlda.bernoulliMH
#' @title Gibbs Sampling for LDA Presence and Absence with Metropolis-Hasting
#' @description Compute the Gibbs Sampling for LDA Presence and Absence with Stick-Breaking Priors
#' @param data - DataFrame with Presence and Absecence (Zeros and Ones)
#' @param n_community - Number of communities
#' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
#' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
#' @param gamma - Hyperparameter  Beta(1,gamma)
#' @param n_gibbs - Total number of Gibbs Samples
#' @param ll_prior - Likelihood compute with Priors ?
#' @param bool display_progress=true - Should I Show the progressBar ?
#' @return Rlda object
#' @export
rlda.bernoulliMH <- function(data, loc.id, n_community, alpha0, alpha1, gamma, n_gibbs, nadapt, ll_prior = TRUE, display_progress = TRUE) {
  # Create a stop point
  stopifnot(inherits(data, "data.frame"))
  stopifnot(inherits(n_community, "numeric"))
  stopifnot(inherits(alpha0, "numeric"))
  stopifnot(inherits(alpha1, "numeric"))
  stopifnot(inherits(gamma, "numeric") | is.na(gamma))
  stopifnot(inherits(n_gibbs, "numeric"))
  stopifnot(inherits(ll_prior, "logical"))
  stopifnot(inherits(display_progress, "logical"))

  # Dictionary
  dat <- data
  a.phi <- alpha0
  b.phi <- alpha1
  nspp <- ncol(data)
  ncomm <- n_community
  nloc <- nrow(data)
  y <- as.matrix(data)
  ngibbs <- n_gibbs

  # initial values convert from a bunch of bernoulli to a single binomial per location
  tmp = aggregate.data(dat, loc.id)
  y = tmp$dat
  loc.id = tmp$loc.id
  nspp = ncol(y)
  nloc = length(unique(loc.id))
  n = tmp$n
  nmat = matrix(n, nloc, nspp)

  # initial values
  theta = matrix(1/ncomm, nloc, ncomm)
  vmat = theta
  vmat[, ncomm] = 1
  phi = matrix(0.2, ncomm, nspp, byrow = T)

  # gibbs stuff
  vec.theta = matrix(0, ngibbs, nloc * ncomm)

  # Remove the loc.id
  vec.phi = matrix(0, ngibbs, ncomm * nspp)
  vec.logl = matrix(NA, ngibbs, 1)
  param = list(theta = theta, phi = phi, vmat = vmat)

  # for MH algorithm
  jump1 = list(vmat = matrix(0.1, nloc, ncomm), phi = matrix(0.1, ncomm, nspp))
  accept1 = list(vmat = matrix(0, nloc, ncomm), phi = matrix(0, ncomm, nspp))
  accept.output = 50

  # create progress bar
  if (display_progress)
    pb <- txtProgressBar(min = 0, max = n_gibbs, style = 3)
  count = 0
  for (i in 1:ngibbs) {
    tmp = update.phiAbundanceSB(param = param, jump = jump1$phi, ncomm = ncomm, nspp = nspp, y = y, nmat = nmat, a.phi = a.phi, b.phi = b.phi)
    param$phi = tmp$phi
    accept1$phi = accept1$phi + tmp$accept

    tmp = update.thetaAbundanceSB(param = param, jump = jump1$vmat, nloc = nloc, ncomm = ncomm, y = y, nmat = nmat, gamma = gamma)
    param$theta = tmp$theta
    param$vmat = tmp$v
    accept1$vmat = accept1$vmat + tmp$accept

    if (i%%accept.output == 0 & i < nadapt) {
      k = print.adapt(parmAccept = accept1, parmJump = jump1, accept.output = accept.output)
      accept1 = k$accept1
      jump1 = k$jump1
    }

    # to assess convergence, examine logl
    prob = get.logl(theta = param$theta, phi = param$phi, y = y, nmat = nmat)
    loglikel = sum(prob)
    if (ll_prior) {
      loglikel = loglikel + sum(dbeta(param$phi, a.phi, b.phi, log = T)) + sum(dbeta(param$vmat[, -ncomm], 1, gamma, log = T))
    }

    vec.logl[i] = loglikel
    vec.theta[i, ] = param$theta

    # Remove the loc.id
    vec.phi[i, ] = param$phi

    # Progress Bar
    if (display_progress)
      setTxtProgressBar(pb, i)
  }
  if (display_progress)
    close(pb)

  # Use a function not exported Execute the LDA for the Bernoulli entry
  res <- list(Theta = vec.theta, Phi = vec.phi, logLikelihood = vec.logl)

  # Type distribution
  res$type <- "Bernoulli"
  # Number of communities
  res$n_community <- n_community
  # Sample size
  res$N <- (length(vec.theta[n_gibbs, ])/ncomm)
  # Covariates
  res$Species <- seq(1, nspp)
  # Alpha0
  res$alpha0 <- alpha0
  # Alpha1
  res$alpha1 <- alpha1
  # Gamma
  res$gamma <- gamma
  # Number of gibbs
  res$n_gibbs <- n_gibbs
  # Species
  res$colnames <- colnames(data)
  # Locations
  res$rownames <- unique(loc.id)
  # Create the class
  class(res) <- c("rlda", "list")
  return(res)
}






#' @name rlda.multinomial
#' @title Gibbs Sampling for LDA Abundance with Stick-Breaking
#' @description Compute the Gibbs Sampling for LDA Abundance with Stick-Breaking
#' @param data - dataFrame with Abundance
#' @param int n_community - Number of communities
#' @param beta - NumericVector for beta (Sx1)
#' @param gamma - Hyperparameter  Beta(1,gamma)
#' @param n_gibbs - Total number of Gibbs Samples
#' @param ll_prior - Likelihood compute with Priors ?
#' @param bool display_progress=true - Should I Show the progressBar ?
#' @return Rlda object
#' @export
rlda.multinomial <- function(data, n_community, beta, gamma, n_gibbs, ll_prior = TRUE, display_progress = TRUE) {
  # Create a stop point
  stopifnot(inherits(data, "data.frame"))
  stopifnot(inherits(n_community, "numeric"))
  stopifnot(inherits(beta, "numeric"))
  stopifnot(inherits(gamma, "numeric") | is.na(gamma))
  stopifnot(inherits(n_gibbs, "numeric"))
  stopifnot(inherits(ll_prior, "logical"))
  stopifnot(inherits(display_progress, "logical"))

  # Use a function not exported Execute the LDA for the Multinomial entry
  res <- lda_multinomial(data, n_community, beta, gamma, n_gibbs, ll_prior, display_progress)
  # Type distribution
  res$type <- "Multinomial"
  # Number of communities
  res$n_community <- n_community
  # Sample size
  res$N <- nrow(data)
  # Covariates
  res$Species <- colnames(data)
  # Beta
  res$beta <- beta
  # Gamma
  res$gamma <- gamma
  # Number of gibbs
  res$n_gibbs <- n_gibbs
  # Species
  res$colnames <- colnames(data)
  # Locations
  res$rownames <- rownames(data)
  # Create the class
  class(res) <- c("rlda", "list")
  return(res)
}

#' @name rlda.binomial
#' @title Compute the Gibbs Sampling for LDA Binomial
#' @description Compute the Gibbs Sampling for LDA Binomial
#' @param DATA - DataFrame with Presence and Absecence (Binomial)
#' @param POP - DataFrame with Population Size (Binomial)
#' @param int n_community - Number of communities
#' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
#' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
#' @param gamma - Hyperparameter  Beta(1,gamma)
#' @param n_gibbs - Total number of Gibbs Samples
#' @param ll_prior - Likelihood compute with Priors ?
#' @param bool display_progress=true - Should I Show the progressBar ?
#' @return Rlda object
#' @export
rlda.binomial <- function(data, pop, n_community, alpha0, alpha1, gamma, n_gibbs, ll_prior = TRUE, display_progress = TRUE) {
  # Create a stop point
  stopifnot(inherits(data, "data.frame"))
  stopifnot(inherits(pop, "data.frame"))
  stopifnot(inherits(n_community, "numeric"))
  stopifnot(inherits(alpha0, "numeric"))
  stopifnot(inherits(alpha1, "numeric"))
  stopifnot(inherits(gamma, "numeric") | is.na(gamma))
  stopifnot(inherits(n_gibbs, "numeric"))
  stopifnot(inherits(ll_prior, "logical"))
  stopifnot(inherits(display_progress, "logical"))
  if (nrow(data) != nrow(pop)) {
    stop("Both \"data\" and \"pop\" must have the same number of rows.")
  }
  # Execute the LDA for the Binomial entry
  res <- lda_binomial(data, pop, n_community, alpha0, alpha1, gamma, n_gibbs, ll_prior, display_progress)
  # Type distribution
  res$type <- "Binomial"
  # Maximum value
  res$max <- max(pop)
  # Number of communities
  res$n_community <- n_community
  # Sample size
  res$N <- nrow(data)
  # Covariates
  res$Species <- colnames(data)
  # Alpha0
  res$alpha0 <- alpha0
  # Alpha1
  res$alpha1 <- alpha1
  # Gamma
  res$gamma <- gamma
  # Number of gibbs
  res$n_gibbs <- n_gibbs
  # Locations
  res$rownames <- rownames(data)
  # Create the class
  class(res) <- c("rlda", "list")
  return(res)
}



rlda.binomialVB <- function(data, loc.id, n_community, alpha0, alpha1, gamma, maxit=1000, thresh=0.0001) {
  # Create a stop point
  stopifnot(inherits(data, "data.frame"))
  stopifnot(inherits(n_community, "numeric"))
  stopifnot(inherits(alpha0, "numeric"))
  stopifnot(inherits(alpha1, "numeric"))
  stopifnot(inherits(gamma, "numeric") | is.na(gamma))
  stopifnot(inherits(maxit, "numeric"))
  stopifnot(inherits(thresh, "numeric"))

  #Initialize
  delta_elbo<-Inf
  nobs<-max(table(data[,loc.id]))

  tmp<-as.matrix(data)
  dat=aggregate.data(tmp, loc.id)
  nloc=nrow(dat)
  nspp=ncol(dat)

  m0=m1=array(abs(rnorm(nloc*nspp*n_community)),dim=c(nloc,nspp,n_community),
              dimnames=list(paste('loc',1:nloc,sep=''),
                            paste('spp',1:nspp,sep=''),
                            paste('comm',1:n_community,sep='')))
  a=b=matrix(1,nloc,n_community)
  c=d=matrix(1,n_community,nspp)


  # Execute the LDA for the Binomial entry
  res <- lda_binomial_var(dat, n_community, maxit, nobs, gamma, alpha0, alpha1, thresh, delta_elbo, m1, m0)
  # Type distribution
  res$type <- "Binomial Variational"
  # Number of communities
  res$n_community <- n_community
  # Sample size
  res$N <- nrow(data)
  # Covariates
  res$Species <- colnames(data)[!colnames(data) %in% loc.id]
  # Alpha0
  res$alpha0 <- alpha0
  # Alpha1
  res$alpha1 <- alpha1
  # Gamma
  res$gamma <- gamma
  # Number of gibbs
  res$n_gibbs <- maxit
  # Locations
  res$rownames <- rownames(data)
  names(res)[3]<-"logLikelihood"
  # Create the class
  class(res) <- c("rlda", "list")
  return(res)
}




aggregate.data=function(dat,id){
  ind<-which(colnames(dat)==id)
  locid=dat[,ind]
  dat1=dat[,-ind]

  nloc=max(locid)
  nspp=ncol(dat1)
  res<- matrix(NA,nloc,nspp)
  for (i in 1:nloc){
    cond=locid==i
    res[i,]=colSums(dat1[cond,])
  }
  return(res)
}




#' @name rlda.binomialMH
#' @title Compute the Gibbs Sampling for LDA Binomial for Remote Sensing
#' @description Compute the Gibbs Sampling for LDA Binomial
#' @param DATA - DataFrame with Presence and Absecence (Binomial)
#' @param POP - DataFrame with Population Size (Binomial)
#' @param int n_community - Number of communities
#' @param alpha0 - Hyperparameter Beta(alpha0,alpha1)
#' @param alpha1 - Hyperparameter Beta(alpha0,alpha1)
#' @param gamma - Hyperparameter  Beta(1,gamma)
#' @param n_gibbs - Total number of Gibbs Samples
#' @param ll_prior - Likelihood compute with Priors ?
#' @param bool display_progress=true - Should I Show the progressBar ?
#' @return Rlda object
#' @export
rlda.binomialMH <- function(data, pop, n_community, alpha0, alpha1, gamma, n_gibbs, ll_prior = TRUE, display_progress = TRUE) {
  # Create a stop point
  stopifnot(inherits(data, "data.frame"))
  stopifnot(inherits(pop, "data.frame"))
  stopifnot(inherits(n_community, "numeric"))
  stopifnot(inherits(alpha0, "numeric"))
  stopifnot(inherits(alpha1, "numeric"))
  stopifnot(inherits(gamma, "numeric") | is.na(gamma))
  stopifnot(inherits(n_gibbs, "numeric"))
  stopifnot(inherits(ll_prior, "logical"))
  stopifnot(inherits(display_progress, "logical"))
  if (nrow(data) != nrow(pop)) {
    stop("Both \"data\" and \"pop\" must have the same number of rows.")
  }

  # Dictionary
  a.omega <- alpha0
  b.omega <- alpha1
  nbands <- ncol(data)
  ncommun <- n_community
  nloc <- nrow(data)
  ngibbs <- n_gibbs
  ndig.values <- as.matrix(pop)
  remote <- as.matrix(data)

  # initial values
  omega = matrix(runif(ncommun * nbands), ncommun, nbands)
  theta = matrix(1/ncommun, nloc, ncommun)
  v = theta
  v[, ncommun] = 1

  # stuff for gibbs sampling
  param = list(theta = theta, omega = omega, v = v, gamma = gamma)
  vec.theta = matrix(NA, ngibbs, nloc * ncommun)
  vec.omega = matrix(NA, ngibbs, ncommun * nbands)

  # stuff for MH algorithm
  jump1 = list(omega = matrix(1, ncommun, nbands), v = matrix(0.3, nloc, ncommun))
  accept1 = list(omega = matrix(0, ncommun, nbands), v = matrix(0, nloc, ncommun))
  accept.output = 50

  # create progress bar
  if (display_progress)
    pb <- txtProgressBar(min = 0, max = ngibbs, style = 3)

  for (i in 1:ngibbs) {
    tmp = update.thetaRemote(remote, param, jump1$v, ncommun, nloc, ndig.values)
    param$theta = tmp$theta
    param$v = tmp$v
    accept1$v = accept1$v + tmp$accept

    tmp = update.omegaRemote(remote, param, jump1$omega, ncommun, nbands, ndig.values, a.omega, b.omega)
    param$omega = tmp$omega
    accept1$omega = accept1$omega + tmp$accept

    if (i%%accept.output == 0 & i < 1000) {
      k = print.adapt(parmAccept = accept1, parmJump = jump1, accept.output = accept.output, FALSE)
      accept1 = k$accept1
      jump1 = k$jump1
    }

    vec.theta[i, ] = param$theta
    vec.omega[i, ] = param$omega
    # Progress Bar
    if (display_progress)
      setTxtProgressBar(pb, i)
  }
  if (display_progress)
    close(pb)

  vec.logl <- rep(0, nloc)

  # Use a function not exported Execute the LDA for the Bernoulli entry
  res <- list(Theta = vec.theta, Phi = vec.omega, logLikelihood = vec.logl)

  # Type distribution
  res$type <- "Binomial"
  # Maximum value
  res$max <- max(pop)
  # Number of communities
  res$n_community <- n_community
  # Sample size
  res$N <- nrow(data)
  # Covariates
  res$Species <- colnames(data)
  # Alpha0
  res$alpha0 <- alpha0
  # Alpha1
  res$alpha1 <- alpha1
  # Gamma
  res$gamma <- gamma
  # Number of gibbs
  res$n_gibbs <- n_gibbs
  # Locations
  res$rownames <- rownames(data)
  # Create the class
  class(res) <- c("rlda", "list")
  return(res)
}


#' Plot the Rlda object.
#'
#' @param x rlda object
#' @param ... ignored
#' @export
#'
#'

plot.rlda <- function(x, burnin = 0.1, maxCluster = NA, ...) {
  old.par <- par(no.readonly = T)
  stopifnot(inherits(burnin, "numeric"))
  stopifnot(!(burnin > 1 || burnin < 0))
  # Burn-in
  i <- ceiling(x$n_gibbs * burnin)
  # Plot the log-likelihood
  if(x$type == "Binomial Variational"){
    plot(x$logLikelihood[i:x$n_gibbs], type = "l", xlab = "Variational iteration", ylab = "ELBO", main = "ELBO")
  }
  else{
    plot(x$logLikelihood[i:x$n_gibbs], type = "l", xlab = "Gibbs iteration", ylab = "Log-Likelihood", main = "Log-Likelihood")
  }

  par(ask = T)
  # Plot the box-plot Theta
  if (is.na(maxCluster))
    maxCluster = x$n_community
  tmp <- colMeans(x$Theta[i:x$n_gibbs, ])
  theta <- matrix(tmp, x$N, maxCluster)
  colnames(theta) = paste("Cluster ", 1:maxCluster, sep = "")
  rownames(theta) = x$rownames
  boxplot(theta, main = "Theta matrix", ylab = "Probability")
  par(mar = c(5.1, 4.1, 4.1, 8.1), ask = T, xpd = TRUE)
  # Plot the box-plot Phi
  tmp <- colMeans(x$Phi[i:x$n_gibbs, ])
  phi <- matrix(tmp, maxCluster, length(x$Species))
  rownames(phi) = paste("Cluster ", 1:maxCluster, sep = "")
  colnames(phi) = x$Species
  #Add
  par(mar = c(5.1, 4.1, 4.1, 2.1), ask = T, xpd = FALSE)
  for (i in 1:maxCluster){
    plot(phi[i,],main=rownames(phi)[i],type='h',ylim=c(0,1),ylab='Probability',
         xaxt='n',xlab='')
    axis(1,at=1:ncol(phi),colnames(phi),las=2)
    abline(h=seq(0,1,by=0.2),lty=3,col='grey')
    abline(v=1:ncol(phi),lty=3,col='grey')
    par(new=TRUE)
    plot(phi[i,],main=rownames(phi)[i],type='h',ylim=c(0,1),
         xaxt='n',xlab='',ylab='')
  }
  invisible(x)
  par(old.par)
}

#' Summarize the Bayesian LDA.
#'
#' @param object rlda object
#' @param ... ignored
#' @export
summary.rlda <- function(object, burnin = 0.1, quantile = 0.95, silent = FALSE, ...) {
  stopifnot(inherits(object, "rlda"))
  stopifnot(inherits(burnin, "numeric"))
  stopifnot(!(burnin > 1 || burnin < 0))
  stopifnot(!(quantile > 1 || burnin < 0))
  stopifnot(inherits(silent, "logical"))

  # Burn-in
  i <- ceiling(object$n_gibbs * burnin)
  seq <- i:object$n_gibbs
  if (!silent) {
    cat(paste("Total number of gibbs sampling:", object$n_gibbs, "\nNumber of clusters:", object$n_community, "\nNumber of variables:", length(object$Species)))
  }
  # Summary Theta (Mean)
  tmp <- colMeans(object$Theta[i:object$n_gibbs, ])
  theta <- matrix(tmp, object$N, object$n_community)
  colnames(theta) = paste("Cluster ", 1:object$n_community, sep = "")
  rownames(theta) = object$rownames
  # Summary Theta (Var)
  tmp <- apply(as.data.frame(object$Theta[i:object$n_gibbs, ]),2,var)
  theta.var <- matrix(tmp, object$N, object$n_community)
  colnames(theta.var) = paste("Cluster ", 1:object$n_community, sep = "")
  rownames(theta.var) = object$rownames

  # Summary Theta (Lower Interval)
  alpha<- (1-quantile)/2
  tmp <- apply(as.data.frame(object$Theta[i:object$n_gibbs, ]),2,function(x) quantile(x,probs=c(alpha)))
  theta.li <- matrix(tmp, object$N, object$n_community)
  colnames(theta.li) = paste("Cluster ", 1:object$n_community, sep = "")
  rownames(theta.li) = object$rownames

  # Summary Theta (Upper Interval)
  tmp <- apply(as.data.frame(object$Theta[i:object$n_gibbs, ]),2,function(x) quantile(x,probs=c(1-alpha)))
  theta.ui <- matrix(tmp, object$N, object$n_community)
  colnames(theta.ui) = paste("Cluster ", 1:object$n_community, sep = "")
  rownames(theta.ui) = object$rownames

  # Summary Phi (Mean)
  tmp <- colMeans(object$Phi[i:object$n_gibbs, ])
  phi <- matrix(tmp, object$n_community, length(object$Species))
  rownames(phi) = paste("Cluster ", 1:object$n_community, sep = "")
  colnames(phi) = object$Species

  #Summary Phi (Var)
  tmp <- apply(as.data.frame(object$Phi[i:object$n_gibbs, ]),2,var)
  phi.var <- matrix(tmp, object$n_community, length(object$Species))
  rownames(phi.var) = paste("Cluster ", 1:object$n_community, sep = "")
  colnames(phi.var) = object$Species

  #Summary Phi (Lower Interval)
  tmp <- apply(as.data.frame(object$Phi[i:object$n_gibbs, ]),2,function(x) quantile(x,probs=c(alpha)))
  phi.li <- matrix(tmp, object$n_community, length(object$Species))
  rownames(phi.li) = paste("Cluster ", 1:object$n_community, sep = "")
  colnames(phi.li) = object$Species

  #Summary Phi (Upper Interval)
  tmp <- apply(as.data.frame(object$Phi[i:object$n_gibbs, ]),2,function(x) quantile(x,probs=c(1-alpha)))
  phi.up <- matrix(tmp, object$n_community, length(object$Species))
  rownames(phi.up) = paste("Cluster ", 1:object$n_community, sep = "")
  colnames(phi.up) = object$Species

  return(list("Theta.mean" = theta, "Phi.mean" = phi,
              "Theta.var" = theta.var, "Phi.var" = phi.var,
              "Theta.li" = theta.li, "Phi.li" = phi.li,
              "Theta.ui" = theta.ui, "Phi.ui" = phi.up))
}


#' Predict the Bayesian LDA.
#'
#' @param object rlda object
#' @param ... ignored
#' @export
predict.rlda <- function(object, data, nclus = NA, burnin = 0.1, places.round = 0, ...) {
  stopifnot(inherits(object, "rlda"))
  stopifnot(inherits(nclus, "numeric"))
  stopifnot(inherits(places.round, "numeric"))
  stopifnot(nclus > 0 || is.na(nclus))
  stopifnot(places.round >= 0)
  stopifnot(inherits(burnin, "numeric"))
  stopifnot(!(burnin > 1 || burnin < 0))

  # Matrix
  summ <- summary.rlda(object, burnin, T)
  phi <- summ$Phi

  if (is.na(nclus)) {
    nclus <- nrow(phi)
  }

  # Create a matrix with all possible combinations of proportions
  seq1 <- seq(from = 0, to = 1, by = 0.05)
  combo <- expand.grid(p1 = seq1)
  for (i in 2:(nclus - 1)) {
    temp <- expand.grid(var = seq1)
    colnames(temp) <- paste0("p", i)
    combo <- merge(combo, temp)
  }
  cond <- apply(combo, 1, sum) <= 1
  combo1 <- combo[cond, ]
  combo1[, paste0("p", nclus)] <- 1 - apply(combo1, 1, sum)


  # Calculate implied binomial probabilities
  probs <- data.matrix(combo1) %*% data.matrix(phi)

  # Import the data for the desired region
  dat1 <- data[, object$Species]
  nbands <- length(object$Species)

  # Let's change the range of our data to start at zero.
  tmp <- apply(dat1, 2, range)
  dat2 <- dat1 - matrix(tmp[1, ], nrow(dat1), nbands, byrow = T)
  tmp <- apply(dat2, 2, range)
  max1 <- tmp[2, ]
  max2 <- matrix(max1, nrow(probs), length(max1), byrow = T)

  # Divisor
  div <- 10^(places.round)
  max2 <- floor(max2/div)
  dat2 <- floor(dat2/div)
  df_args <- c(as.data.frame(dat2), sep = "")
  dat2full <- as.data.frame(dat2)
  dat2full$ID <- as.character(do.call(paste, df_args))
  dat2full$Sort <- seq(1, nrow(dat2full))
  dat2 <- unique(dat2)

  # Keep the same scale
  max2 <- max2 * div
  dat2 <- dat2 * div

  ncl <- detectCores()
  cl <- makeCluster(ncl)
  registerDoParallel(cl)
  # find which proportion of endmembers that yields the highest likelihood
  res2 <- foreach(i = 1:nrow(dat2), .combine = rbind) %dopar% {
    rasc = matrix(as.integer(dat2[i, ]), nrow = nrow(probs), ncol = ncol(dat2), byrow = T)
    llikel = dbinom(rasc, size = max2, prob = probs, log = T)
    fim = apply(llikel, 1, sum)
    ind = which(fim == max(fim))
    as.numeric(combo1[ind, ])
  }
  colnames(res2) = paste("prop", 1:nclus, sep = "")
  rownames(res2) = NULL
  # Stop clusters
  stopCluster(cl)

  # Convert to data.frame
  dat2 <- as.data.frame(dat2)
  df_args <- c(dat2, sep = "")
  res2 <- as.data.frame(res2)
  res2$ID <- as.character(do.call(paste, df_args))
  final <- merge(dat2full, res2, by = "ID", all = T)
  final <- final[, -which(names(final) %in% c("ID"))]
  final <- final[order(final$Sort), ]
  final <- final[, paste("prop", 1:nclus, sep = "")]
  return(final)
}

logLik.rlda <- function(object, ...) {
  return(object$logLikelihood)
}

print.rlda <- function(x, burnin = 0.1, ...) {
  stopifnot(inherits(x, "rlda"))
  stopifnot(inherits(burnin, "numeric"))
  stopifnot(!(burnin > 1 || burnin < 0))

  # Burn-in
  i <- ceiling(x$n_gibbs * burnin)
  seq <- i:x$n_gibbs

  # Summary Theta
  tmp <- colMeans(x$Theta[i:x$n_gibbs, ])
  theta <- matrix(tmp, x$N, x$n_community)
  colnames(theta) = paste("Cluster ", 1:x$n_community, sep = "")
  rownames(theta) = x$rownames
  # Summary Phi
  tmp <- colMeans(x$Phi[i:x$n_gibbs, ])
  phi <- matrix(tmp, x$n_community, length(x$Species))
  rownames(phi) = paste("Cluster ", 1:x$n_community, sep = "")
  colnames(phi) = x$Species
  #Get the logLik
  logLikeli= logLik.rlda(x)
  return(list(Theta = theta, Phi = phi, logLik=logLikeli))
}



getTheta.rlda <- function(object, burnin = 0.1, ...) {
  stopifnot(inherits(object, "rlda"))
  stopifnot(inherits(burnin, "numeric"))
  stopifnot(!(burnin > 1 || burnin < 0))
  # Burn-in
  i <- ceiling(object$n_gibbs * burnin)
  seq <- i:object$n_gibbs
  # Summary Theta
  tmp <- colMeans(object$Theta[i:object$n_gibbs, ])
  theta <- matrix(tmp, object$N, object$n_community)
  colnames(theta) = paste("Cluster ", 1:object$n_community, sep = "")
  rownames(theta) = object$rownames
  return(theta)
}

getPhi.rlda <- function(object, burnin = 0.1, ...) {
  stopifnot(inherits(object, "rlda"))
  stopifnot(inherits(burnin, "numeric"))
  stopifnot(!(burnin > 1 || burnin < 0))
  # Burn-in
  i <- ceiling(object$n_gibbs * burnin)
  seq <- i:object$n_gibbs
  # Summary Phi
  tmp <- colMeans(object$Phi[i:object$n_gibbs, ])
  phi <- matrix(tmp, object$n_community, length(object$Species))
  rownames(phi) = paste("Cluster ", 1:object$n_community, sep = "")
  colnames(phi) = object$Species
  return(phi)
}

rlda2mcmc.rlda<-function(object, ...){
  stopifnot(inherits(object, "rlda"))

  #Theta object
  thetaMCMC <- object$Theta
  ss1<- paste("Cluster ", 1:object$n_community, sep = "")
  ss2<- object$rownames
  colnames(thetaMCMC) <- apply(expand.grid(ss2,ss1),1,function(x) paste(x,collapse=" - "))
  rownames(thetaMCMC) <- paste("Gibbs ", 1:nrow(thetaMCMC), sep = "")
  #Casting
  Theta <- coda::mcmc(thetaMCMC)

  #Phi object
  phiMCMC <- object$Phi
  ss1<- paste("Cluster ", 1:object$n_community, sep = "")
  ss2<- object$Species
  colnames(phiMCMC) <- apply(expand.grid(ss1,ss2),1,function(x) paste(x,collapse=" - "))
  rownames(phiMCMC) <- paste("Gibbs ", 1:nrow(phiMCMC), sep = "")
  #Casting
  Phi <- coda::mcmc(phiMCMC)
  return(list("Theta"=Theta,"Phi"=Phi))
}

generateMultinomialLDA.rlda<-function(seed0, community, variables, observations, totalElements, beta, gamma, ...){
  #Generate the fake data
  set.seed(seed0)

  #Number of communities
  stopifnot(variables==length(beta))
  stopifnot(all(beta>0))
  stopifnot(length(gamma)==1)
  stopifnot(gamma>0)
  stopifnot(community>=2)
  stopifnot(variables>=community)
  stopifnot(totalElements>=1)
  stopifnot(observations>=1)
  #Number of Species
  species<-variables

  #Number of locations
  locations<-observations

  #Number of individuals in each location
  size<-rpois(locations,totalElements)

  #Generate Phi
  Phi<-gtools::rdirichlet(community,beta)

  #Generate V
  vMat<-matrix(rbeta(locations*community,1,gamma),nrow=locations,ncol=community)
  vMat[,community]<-1

  #Generate Theta
  Theta<-apply(vMat,1,function(x){
    prod<-1;
    Theta<-rep(NA,community)
    for(c in 1:community){
      vNumber <- x[c];
      if (c == 1) prod<-1;
      if (c >  1) prod<-prod*(1.0-x[c-1]);
      Theta[c]<-vNumber*prod;
    }
    Theta
  })
  Theta<-t(Theta)

  #Generate Z
  tmp<-cbind(size,Theta)
  Z<-t(apply(tmp,1,function(x)rmultinom(1, x[1], x[2:(community+1)])))

  #Generate S
  FIA<-matrix(NA,locations,species)
  for(l in 1:locations){
    tmp<-rep(0,species)
    for(c in 1:community){
      tmp<-tmp+t(as.data.frame(rmultinom(1, Z[l,c], Phi[c,])))
    }
    FIA[l,]<-tmp
  }

  #Create rownames colnames
  rownames(FIA)<-paste0("Location ",seq(1,nrow(FIA)))
  colnames(FIA)<-paste0("Specie ",seq(1,ncol(FIA)))
  return(list("Theta"=Theta,"Phi"=Phi, "Z"=Z, "Data"= FIA))
}

generateBernoulliLDA.rlda<-function(seed0, community, variables, observations, alpha0, alpha1, gamma, ...){
  #Generate the fake data
  set.seed(seed0)

  #Number of communities
  stopifnot(length(gamma)==1)
  stopifnot(gamma>0)
  stopifnot(length(alpha0)==1)
  stopifnot(alpha0>0)
  stopifnot(length(alpha1)==1)
  stopifnot(alpha1>0)
  stopifnot(community>=2)
  stopifnot(variables>=community)
  stopifnot(observations>=1)

  #Number of Species
  species<-variables

  #Number of locations
  locations<-observations

  Theta<-matrix(rbeta(species*community,alpha0,alpha1),nrow=species,ncol=community)

  #Generate V
  gamma<-0.01
  vMat<-matrix(rbeta(locations*community,1,gamma),nrow=locations,ncol=community)
  vMat[,community]<-1

  #Generate Phi
  Phi<-apply(vMat,1,function(x){
    prod<-1;
    Phi<-rep(NA,community)
    for(c in 1:community){
      vNumber <- x[c];
      if (c == 1) prod<-1;
      if (c >  1) prod<-prod*(1.0-x[c-1]);
      Phi[c]<-vNumber*prod;
    }
    Phi
  })
  Phi<-t(Phi)

  #Generate Z
  Z<-t(apply(Phi,1,function(x)rmultinom(1,1,x)))

  #Generate Data
  tmp<-as.data.frame(t(rep(0,species)))
  DATA<-data.frame()
  for(l in 1:locations){
    for(s in 1:species){
      tmp[s]<-rbinom(1,1,sum(Theta[s,]*Phi[l,]))
    }
    DATA<-rbind(DATA,tmp)
  }

  #Create rownames colnames
  rownames(DATA)<-paste0("Location ",seq(1,nrow(DATA)))
  colnames(DATA)<-paste0("Specie ",seq(1,ncol(DATA)))



  return(list("Theta"=Theta,"Phi"=Phi, "Z"=Z, "Data"= DATA))
}

generateBinomialLDA.rlda<-function(seed0, community, variables, observations, totalElements, alpha0, alpha1, ...){
  #Generate the fake data
  set.seed(seed0)

  #Number of communities
  stopifnot(length(alpha0)==1)
  stopifnot(alpha0>0)
  stopifnot(length(alpha1)==1)
  stopifnot(alpha1>0)
  stopifnot(community>=2)
  stopifnot(variables>=community)
  stopifnot(totalElements>=1)
  stopifnot(observations>=1)

  #Number of Species
  bands<-variables

  #Number of locations
  locations<-observations

  #Generate the Omega
  Omega<-matrix(rbeta(bands*community,alpha0,alpha1),nrow=bands,ncol=community)
  Theta=matrix(rbeta(locations*community,alpha0,alpha1),nrow=locations,ncol=community)
  Theta=t(apply(Theta,1,function(x) x/sum(x)))

  #Generate POP matrix
  pop<-matrix(totalElements,nrow=locations,ncol=bands)

  #Generate Data
  tmp<-rep(0,bands)
  DATA<-matrix(NA,locations,bands)
  for(l in 1:locations){
    for(b in 1:bands){
      prob=sum(Theta[l,]*Omega[b,])
      tmp[b]<-rbinom(1,pop[l,b],prob)
    }
    DATA[l,]<-tmp
  }

  #Create rownames colnames
  rownames(DATA)<-paste0("Location ",seq(1,nrow(DATA)))
  colnames(DATA)<-paste0("Bands ",seq(1,ncol(DATA)))
  POP<-as.data.frame(pop)
  rownames(POP)<-paste0("Location ",seq(1,nrow(POP)))
  colnames(POP)<-paste0("Bands ",seq(1,ncol(POP)))


  return(list("Theta"=Theta, "Phi"=Omega, "Pop"=POP, "Data"= DATA))
}





