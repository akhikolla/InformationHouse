simLong <- function(n = 100,
                    ntest = 0,
                    N = 5,
                    rho = 0.8,
                    model = c(1, 2),
                    phi = 1,
                    q_x = 0,
                    q_y = 0,
                    type = c("corCompSym", "corAR1", "corSymm", "iid"))
{
  dta <- data.frame(do.call("rbind", lapply(1:(n+ntest), function(i) {
    Ni <- round(runif(1, 1, 3 * N))
    type <- match.arg(type, c("corCompSym", "corAR1", "corSymm", "iid"))
    if (type == "corCompSym") {
      corr.mat <- matrix(rho, nrow=Ni, ncol=Ni)
      diag(corr.mat) <- 1
    }
    if (type == "corAR1") {
      corr.mat <- diag(rep(1, Ni))
      if (Ni > 1) {
        for (ii in 1:(Ni - 1)) {
          corr.mat[ii, (ii + 1):Ni] <- rho^(1:(Ni - ii))
        }
        ind <- lower.tri(corr.mat)
        corr.mat[ind] <- t(corr.mat)[ind]
      }
    }
    if (type == "iid") {
      corr.mat <- diag(rep(1, Ni))
    }
    tm <- sort(sample((1:(3 * N))/N, size = Ni, replace = TRUE))
    if (model == 1) {
      x1 <- rnorm(Ni)
      x2 <- rnorm(Ni)
      x3 <- rnorm(Ni)
      x4 <- rnorm(Ni)
      x <- cbind(x1, x2, x3, x4)
      p <- ncol(x)
      if (q_x > 0) {
        xnoise <- matrix(rnorm(Ni*q_x),ncol = q_x)
        x <- cbind(x, xnoise)
      }
      eps1 <- sqrt(phi) * t(chol(corr.mat)) %*% rnorm(Ni)
      eps2 <- sqrt(phi) * t(chol(corr.mat)) %*% rnorm(Ni)
      eps3 <- sqrt(phi) * t(chol(corr.mat)) %*% rnorm(Ni)
      y1 <- 1.5 + 1.5 * x1 + 1 * x4 + eps1
      y2 <- 1.5 + 1.5 * x2 + 1 * x4 + eps2
      y3 <- 1.5 + 1.5 * x3 + 1 * x4 + eps3
      y <- cbind(y1,y2,y3)
      if (q_y > 0) {
        ynoise <- matrix(rnorm(Ni*q_y),ncol = q_y)
        y <- cbind(y,ynoise)
      }
    }
    if (model == 2) {
      x1 <- rnorm(Ni)
      x2 <- rnorm(Ni)
      x3 <- rnorm(Ni)
      x4 <- rnorm(Ni)
      x <- cbind(x1, x2, x3, x4)
      p <- ncol(x)
      if (q_x > 0) {
        xnoise <- matrix(rnorm(Ni*q_x),ncol = q_x)
        x <- cbind(x, xnoise)
      }
      eps <- sqrt(phi) * t(chol(corr.mat)) %*% rnorm(Ni)
      b.x1 <- ifelse(tm < 1.5, 0, 1.5 * x1 )
      y1 <- 1.5 + b.x1 - 0.65 * (x2^2) * (tm^2) - 2.5 * x3 - 1.5 * exp(x4) + eps
      y <- cbind(y1)
      if (q_y > 0) {
        ynoise <- matrix(rnorm(Ni*q_y),ncol = q_y)
        y <- cbind(y,ynoise)
      }
    }
    if (model == 3) {
      x1 <- sort(rnorm(Ni))
      x2_Temp <- sort(runif(n = Ni,min = 0,max = 1))
      x2 <- unlist(lapply(1:Ni,function(ii){ if(tm[ii] > 1 && tm[ii] < 2) -1*x2_Temp[ii] else x2_Temp[ii] }))
      x3 <- rnorm(1)
      x4 <- rnorm(1)
      x <- cbind(x1, x2, x3, x4)
      B1 <- 2
      B2 <- tm^2
      B3 <- exp(tm)
      B4 <- 1
      p <- ncol(x)
      if (q_x > 0) {
        xnoise <- matrix(rnorm(Ni*q_x),ncol = q_x)
        x <- cbind(x, xnoise)
      }
      eps <- sqrt(phi) * t(chol(corr.mat)) %*% rnorm(Ni)
      y1 <- 1.5 + B1 * x1 + B2 * x2 + B3 * x3 + B4 * x4 + eps
      y <- cbind(y1)
      if (q_y > 0) {
        ynoise <- matrix(rnorm(Ni*q_y),ncol = q_y)
        y <- cbind(y,ynoise)
      }
    }
    cbind(x,tm, rep(i, Ni), y)
  })))
  d_x <- q_x + 4
  if(model == 1){
    d_y <- q_y + 3
  }
  if(model == 2){
    d_y <- q_y + 1
  }
  if(model == 3){
    d_y <- q_y + 1
  }
  colnames(dta) <- c(paste("x", 1:d_x, sep = ""), "time", "id", paste("y", 1:d_y, sep = ""))
  dtaL <- list(features = dta[, 1:d_x], time = dta$time, id = dta$id, y = dta[, -c(1:(d_x+2)) ])
  trn <- c(1:sum(dta$id <= n))
  return(invisible(list(dtaL = dtaL, dta = dta, trn = trn)))
}
