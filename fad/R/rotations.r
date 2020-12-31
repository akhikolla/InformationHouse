
varimax <- function(x, normalize = TRUE, eps = 1e-5)
{
  nc <- ncol(x)
  if(nc < 2) return(x)
  if(normalize) {
    sc <- sqrt(drop(apply(x, 1L, function(x) sum(x^2))))
    x <- x/sc
  }
  p <- nrow(x)
  TT <- diag(nc)
  d <- 0
  for(i in 1L:1000L) {
    z <- x %*% TT
    B  <- crossprod(x, z^3 - z %*% diag(drop(rep(1, p) %*% z^2))/p)
    sB <- La.svd(B)
    TT <- sB$u %*% sB$vt
    dpast <- d
    d <- sum(sB$d)
    if(d < dpast * (1 + eps)) break
  }
  z <- x %*% TT
  if(normalize) z <- z * sc
  dimnames(z) <- dimnames(x)
  class(z) <- "loadings"
  list(loadings = z, rotmat = TT)
}

promax <- function(x, m = 4)
{
  if(ncol(x) < 2) return(x)
  dn <- dimnames(x)
  xx <- varimax(x)
  x <- xx$loadings
  Q <- x * abs(x)^(m-1)
  U <- lm.fit(x, Q)$coefficients
  d <- diag(solve(crossprod(U)))
  U <- .postmdiag(U , sqrt(d))
  dimnames(U) <- NULL
  z <- x %*% U
  U <- xx$rotmat %*% U
  dimnames(z) <- dn
  class(z) <- "loadings"
  list(loadings = z, rotmat = U)
}


quartimax <- function(x, normalize=FALSE, eps=1e-5, maxit=1000) {
  L <- x
  al <- 1
  Tmat <- diag(ncol(x))
  f <- -sum(colSums(L^4))/4
  Gq <- -L^3
  G <- crossprod(L,Gq)

  for (iter in 0:maxit){
    M <- crossprod(Tmat,G)
    S <- (M + t(M))/2
    Gp <- G - Tmat %*% S
    s <- sqrt(sum(colSums(Gp^2)))
    if (s < eps)  break
    al <- 2*al
    for (i in 0:10){
      X <- Tmat - al*Gp
      UDV <- svd(X)
      Tmatt <- UDV$u %*% t(UDV$v)
      L <- x %*% Tmatt
      ft <- -sum(colSums(L^4))/4
      if (ft < (f - 0.5*s^2*al)) break
      al <- al/2
    }
    Tmat <- Tmatt
    f <- ft
    G <- crossprod(x,-L^3)
  }
  convergence <- (s < eps)
  if ((iter == maxit) & !convergence)
    warning("convergence not obtained in quartimax. ", maxit, " iterations used.")
  dimnames(L) <- dimnames(x)

  r <- list(loadings=L, rotmat=solve(t(Tmat)), orthogonal=TRUE, convergence=convergence, Gq=-L^3)
  return(r);
}


equamax <- function(x, normalize=FALSE, eps=1e-5, maxit=1000) {
  L <- x
  al <- 1
  Tmat <- diag(ncol(x))
  #N <- matrix(1,q,q)-diag(q)
  #M <- matrix(1,p,p)-diag(p)
  kappa <- ncol(L)/(2*nrow(L))
  q <- ncol(L)
  p <- nrow(L)
  L2 <- L^2
  rL2 <- matrix(rep(rowSums(L2),q),nrow = p,byrow=F) - L2 # L2 %*% N
  f1 <- (1-kappa)*sum(colSums(L2*rL2))/4
  cL2 <- matrix(rep(colSums(L2),p),ncol = q,byrow=T) - L2 # M %*% L2
  f2 <- kappa*sum(colSums(L2*cL2))/4

  Gq <- (1-kappa)*L*rL2 + kappa*L*cL2
  f <- f1 + f2
  G <- crossprod(L,Gq)

  for (iter in 0:maxit){
    M <- crossprod(Tmat,G)
    S <- (M + t(M))/2
    Gp <- G - Tmat %*% S
    s <- sqrt(sum(colSums(Gp^2)))
    if (s < eps)  break
    al <- 2*al
    for (i in 0:10){
      X <- Tmat - al*Gp
      UDV <- svd(X)
      Tmatt <- UDV$u %*% t(UDV$v)
      L <- x %*% Tmatt
      L2 <- L^2
      rL2 <- matrix(rep(rowSums(L2),q),nrow = p,byrow=F) - L2 # L2 %*% N
      f1 <- (1-kappa)*sum(colSums(L2*rL2))/4
      cL2 <- matrix(rep(colSums(L2),p),ncol = q,byrow=T) - L2 # M %*% L2
      f2 <- kappa*sum(colSums(L2*cL2))/4
      ft <- f1+f2
      if (ft < (f - 0.5*s^2*al)) break
      al <- al/2
    }
    Tmat <- Tmatt
    f <- ft
    Gq <- (1-kappa)*L*rL2 + kappa*L*cL2
    G <- crossprod(x,Gq)
  }
  convergence <- (s < eps)
  if ((iter == maxit) & !convergence)
    warning("convergence not obtained in quartimax. ", maxit, " iterations used.")
  dimnames(L) <- dimnames(x)

  r <- list(loadings=L, rotmat=t(solve(Tmat)), orthogonal=TRUE, convergence=convergence, Gq=Gq)
  return(r);
}


