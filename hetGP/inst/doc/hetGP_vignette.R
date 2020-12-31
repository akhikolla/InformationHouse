## ----include=FALSE-------------------------------------------------------
library("knitr")
## cache can be set to TRUE 
opts_chunk$set(
  engine='R', tidy=FALSE, cache=FALSE, autodep=TRUE
)
render_sweave() # For JSS when using knitr
knitr::opts_chunk$set(fig.pos = 'ht!')

## ----preliminaries, echo=FALSE, results='hide'----------------------
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE, scipen = 5)

## ----nl,message=FALSE-----------------------------------------------
library("hetGP")
nLL <- function(par, X, Y) { 
  theta <- par[1:ncol(X)] 
  tau2 <- par[ncol(X) + 1] 
  n <- length(Y) 
  K <- cov_gen(X1 = X, theta = theta) + diag(tau2, n) 
  Ki <- solve(K) 
  ldetK <- determinant(K, logarithm = TRUE)$modulus 
  (n / 2) * log(t(Y) %*% Ki %*% Y) + (1 / 2) * ldetK 
}

## ----gnl------------------------------------------------------------
gnLL <- function(par, X, Y) {
  n <- length(Y)
  theta <- par[1:ncol(X)]; tau2 <- par[ncol(X) + 1]
  K <- cov_gen(X1 = X, theta = theta) + diag(tau2, n)
  Ki <- solve(K); KiY <- Ki %*% Y
  dlltheta <- rep(NA, length(theta))
  for(k in 1:length(dlltheta)) {
    dotK <- K * as.matrix(dist(X[, k]))^2 / (theta[k]^2)
    dlltheta[k] <- n * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - 
      sum(diag(Ki %*% dotK))
  }
  dlltau2 <- n * t(KiY) %*% KiY / (t(Y) %*% KiY) - sum(diag(Ki))
  -c(dlltheta / 2, dlltau2 / 2)
}

## ----exp2d----------------------------------------------------------
library("lhs")
X <- 6 * randomLHS(40, 2) - 2
X <- rbind(X, X)
y <- X[, 1] * exp(-X[, 1]^2 - X[, 2]^2) + rnorm(nrow(X), sd = 0.01)

## ----exp2doptim-----------------------------------------------------
Lwr <- sqrt(.Machine$double.eps); Upr <- 10
out <- optim(c(rep(0.1, 2), 0.1 * var(y)), nLL, gnLL, method = "L-BFGS-B",
  lower = Lwr, upper = c(rep(Upr, 2), var(y)), X = X, Y = y)
out$par

## ----pred1----------------------------------------------------------
Ki <- solve(cov_gen(X, theta = out$par[1:2]) + diag(out$par[3], nrow(X)))
nuhat <- drop(t(y) %*% Ki %*% y / nrow(X))

## ----xx-------------------------------------------------------------
xx <- seq(-2, 4, length = 40)
XX <- as.matrix(expand.grid(xx, xx))

## ----pred-----------------------------------------------------------
KXX <- cov_gen(XX, theta = out$par[1:2]) + diag(out$par[3], nrow(XX))
KX <- cov_gen(XX, X, theta = out$par[1:2])
mup <- KX %*% Ki %*% y
Sigmap <- nuhat * (KXX - KX %*% Ki %*% t(KX))

## ----exp2dp, echo = FALSE, fig.height=6, fig.width=12, fig.align='center', fig.cap="\\label{fig:exp2d}Example predictive surface from a GP.  Open circles are the training locations."----
library("colorspace")
sdp <- sqrt(diag(Sigmap))
par(mfrow = c(1,2))
cols <- sequential_hcl(palette = "Viridis", n = 128, l = c(40, 90))
persp(xx, xx, matrix(mup, ncol = length(xx)), theta = -30, phi = 30,
  main = "mean surface", xlab = "x1", ylab = "x2", zlab = "y")
image(xx, xx, matrix(sdp, ncol = length(xx)), main = "variance", 
  xlab = "x1", ylab = "x2", col = cols)
points(X[, 1], X[, 2])

## ----library--------------------------------------------------------
fit <- mleHomGP(X, y, rep(Lwr, 2), rep(Upr, 2), known = list(beta0 = 0),
  init = c(list(theta = rep(0.1, 2), g = 0.1 * var(y))))
c(fit$theta, fit$g)

## ----motofit--------------------------------------------------------
library("MASS")
hom <- mleHomGP(mcycle$times, mcycle$accel, covtype = "Matern5_2")
het <- mleHetGP(mcycle$times, mcycle$accel, covtype = "Matern5_2")

## ----motopred-------------------------------------------------------
Xgrid <- matrix(seq(0, 60, length = 301), ncol = 1)
p <- predict(x = Xgrid, object = hom)
p2 <- predict(x = Xgrid, object = het)

## ----motofig, echo=FALSE,fig.height=6, fig.width=7, out.width='4in', fig.align='center', fig.cap="\\label{fig:moto1}Homoskedastic (solid red) versus heteroskedastic (dashed blue) GP fits to the motorcycle data via mean (thick) and 95\\% error bars (thin). Open circles mark the actual data, dots are averaged observations $\\yu$ with corresponding error bars from the empirical variance (when $a_i > 0$)."----
plot(mcycle, main = "Predictive Surface", ylim = c(-160, 90),
  ylab = "acceleration (g)", xlab = "time (ms)")
lines(Xgrid, p$mean, col = 2, lwd = 2)
lines(Xgrid, qnorm(0.05, p$mean, sqrt(p$sd2 + p$nugs)), col = 2)
lines(Xgrid, qnorm(0.95, p$mean, sqrt(p$sd2 + p$nugs)), col = 2)
lines(Xgrid, p2$mean, col = 4, lwd = 2, lty = 4)
lines(Xgrid, qnorm(0.05, p$mean, sqrt(p2$sd2 + p2$nugs)), col = 4, lty = 4)
lines(Xgrid, qnorm(0.95, p$mean, sqrt(p2$sd2 + p2$nugs)), col = 4, lty = 4)
empSd <- sapply(find_reps(mcycle[, 1], mcycle[, 2])$Zlist, sd)
points(het$X0, het$Z0, pch = 20)
arrows(x0 = het$X0, y0 = qnorm(0.05, het$Z0, empSd), 
  y1 = qnorm(0.95, het$Z0, empSd), code = 3, angle = 90, length = 0.01)

## ----sirdesign------------------------------------------------------
Xbar <- randomLHS(200, 2)
a <- sample(1:100, nrow(Xbar), replace = TRUE)
X <- matrix(NA, ncol = 2, nrow = sum(a))
nf <- 0
for(i in 1:nrow(Xbar)) {
  X[(nf + 1):(nf + a[i]),] <- matrix(rep(Xbar[i,], a[i]), ncol = 2,
    byrow = TRUE)
  nf <- nf + a[i]
}

## ----sireval--------------------------------------------------------
Y <- apply(X, 1, sirEval)

## ----sirfit---------------------------------------------------------
fit <- mleHetGP(X, Y, lower = rep(0.05, 2), upper = rep(10, 2), 
  settings = list(linkThetas = "none"), covtype = "Matern5_2", maxit = 1e4)

## ----sirpred, echo = FALSE------------------------------------------
xx <- seq(0, 1, length = 100)
XX <- as.matrix(expand.grid(xx, xx))
p <- predict(fit, XX)

## ----sirvis, echo = FALSE, fig.height=6, fig.width=12, fig.align='center', fig.cap="\\label{fig:sir}Heteroskedastic GP fit to SIR data.  Left panel shows the predictive mean surface; right panel shows the estimated standard deviation.  Text in both panels shows numbers of replicates."----
psd <- sqrt(p$sd2 + p$nugs)
par(mfrow = c(1, 2))
image(xx, xx, matrix(p$mean, 100), xlab = "S0", ylab = "I0", col = cols,
  main = "Mean Infected")
text(Xbar, labels = a, cex = 0.75)
image(xx, xx, matrix(psd, 100), xlab = "S0", ylab = "I0", col = cols,
  main = "SD Infected")
text(Xbar, labels = a, cex = 0.75)

## ----loadbf---------------------------------------------------------
data("bfs")
thetas <- matrix(bfs.exp$theta, ncol = 1)
bfs <- as.matrix(t(bfs.exp[, -1]))

## ----fitbf----------------------------------------------------------
bfs1 <- mleHetTP(X = list(X0 = log10(thetas), Z0 = colMeans(log(bfs)),
  mult = rep(nrow(bfs), ncol(bfs))), Z = log(as.numeric(bfs)),
  lower = 10^(-4), upper = 5, covtype = "Matern5_2")

## ----predbf, echo = FALSE-------------------------------------------
dx <- seq(0, 1, length = 100)
dx <- 10^(dx * 4 - 3)
p <- predict(bfs1, matrix(log10(dx), ncol = 1))

## ----visbf, echo=FALSE,fig.height=6, fig.width=12, fig.align='center', fig.cap="Left: heteroskedastic TP fit to the Bayes factor data under exponential hyperprior. Right: output given by the \\code{plot} method."----
par(mfrow = c(1, 2))
matplot(log10(thetas), t(log(bfs)), col = 1, pch = 21, ylab = "log(bf)", 
  main = "Bayes factor surface")
lines(log10(dx), p$mean, lwd = 2, col = 2)
lines(log10(dx), p$mean + 2 * sqrt(p$sd2 + p$nugs), col = 2, lty = 2,
  lwd = 2)
lines(log10(dx), p$mean + 2 * sqrt(p$sd2), col = 4, lty = 3, lwd = 2)
lines(log10(dx), p$mean - 2 * sqrt(p$sd2 + p$nugs), col = 2, lty = 2,
  lwd = 2)
lines(log10(dx), p$mean - 2 * sqrt(p$sd2), col = 4, lty = 3, lwd = 2)
legend("topleft", c("hetTP mean", expression(paste("hetTP interval on Y(x)|", D[N])), "hetTP interval on f(x)"), col = c(2,2,4), lty = 1:3,
  lwd = 2)
plot(bfs1)
par(mfrow = c(1,1))

## ----loadbf2--------------------------------------------------------
D <- as.matrix(bfs.gamma[, 1:2])
bfs <- as.matrix(t(bfs.gamma[, -(1:2)]))

## ----fitbf2---------------------------------------------------------
bfs2 <- mleHetTP(X = list(X0 = log10(D), Z0 = colMeans(log(bfs)), 
  mult = rep(nrow(bfs), ncol(bfs))), Z = log(as.numeric(bfs)), 
  lower = rep(10^(-4), 2), upper = rep(5, 2), covtype = "Matern5_2")

## ----predbf2,echo=FALSE---------------------------------------------
DD <- as.matrix(expand.grid(dx, dx))
p <- predict(bfs2, log10(DD))

## ----visbf2, echo = FALSE, fig.height=6, fig.width=12, fig.align='center', fig.cap="Heteroskedastic TP fit to the Bayes factor data under Gamma hyperprior."----
par(mfrow = c(1, 2))
mbfs <- colMeans(bfs)
image(log10(dx), log10(dx), t(matrix(p$mean, ncol=length(dx))),  
  col = cols, xlab = "log10 alpha", ylab = "log10 beta", 
  main = "mean log BF")
text(log10(D[, 2]), log10(D[, 1]), signif(log(mbfs), 2), cex = 0.75)
contour(log10(dx), log10(dx), t(matrix(p$mean, ncol = length(dx))),
  levels = c(-5, -3, -1, 0, 1, 3, 5), add = TRUE, col = 4)
image(log10(dx), log10(dx), t(matrix(sqrt(p$sd2 + p$nugs), 
  ncol = length(dx))), col = cols, xlab = "log10 alpha", 
  ylab = "log10 beta", main = "sd log BF")
text(log10(D[, 2]), log10(D[, 1]), signif(apply(log(bfs), 2, sd), 2),
  cex = 0.75)

## ----atoload--------------------------------------------------------
data("ato")  

## ----atotime--------------------------------------------------------
c(n = nrow(Xtrain), N = length(unlist(Ztrain)), time = out$time)

## ----atotestscore---------------------------------------------------
sc <- scores(out, Xtest, matrix(unlist(Ztest), byrow = TRUE, ncol = 10))

## ----atotrainscore--------------------------------------------------
sc.out <- scores(model = out, Xtest = Xtrain.out, Ztest = Ztrain.out)

## ----atobothscore---------------------------------------------------
c(test = mean(sc), train = mean(sc.out), combined = mean(c(sc, sc.out)))

## ----twors----------------------------------------------------------
rn <- c(4.5, 5.5, 6.5, 6, 3.5)
X0 <- matrix(seq(0.05, 0.95, length.out = length(rn)))
X1 <- matrix(c(X0, 0.2, 0.4))
Y1 <- c(rn, 5.2, 6.3)
r1 <- splinefun(x = X1, y = Y1, method = "natural")
X2 <- matrix(c(X0, 0.0, 0.3))
Y2 <- c(rn, 7, 4)
r2 <- splinefun(x = X2, y = Y2, method = "natural")

## ----twovarsXX------------------------------------------------------
XX <- matrix(seq(0, 1, by = 0.005))

## ----imspe.r--------------------------------------------------------
IMSPE.r <- function(x, X0, theta, r) {
  x <- matrix(x, nrow = 1)
  Wijs <- Wij(mu1 = rbind(X0, x), theta = theta, type = "Gaussian")
  K <- cov_gen(X1 = rbind(X0, x), theta = theta)
  K <- K + diag(apply(rbind(X0, x), 1, r))
  return(1 - sum(solve(K) * Wijs))
}

## ----twoimspe-------------------------------------------------------
imspe1 <- apply(XX, 1, IMSPE.r, X0 = X0, theta = 0.25, r = r1)
imspe2 <- apply(XX, 1, IMSPE.r, X0 = X0, theta = 0.25, r = r2)
xstar1 <- which.min(imspe1)
xstar2 <- which.min(imspe2)

## ----rx-------------------------------------------------------------
rx <- function(x, X0, rn, theta, Ki, kstar, Wijs) {
  x <- matrix(x, nrow = 1)
  kn1 <- cov_gen(x, X0, theta = theta)
  wn <- Wij(mu1 = x, mu2 = X0, theta = theta, type = "Gaussian")
  a <- kn1 %*% Ki %*% Wijs %*% Ki %*% t(kn1) - 2 * wn %*% Ki %*% t(kn1) 
  a <- a + Wij(mu1 = x, theta = theta, type = "Gaussian")
  Bk <- tcrossprod(Ki[, kstar], Ki[kstar,]) / 
    (2 / rn[kstar] - Ki[kstar, kstar])
  b <- sum(Bk * Wijs)
  sn <- 1 - kn1 %*% Ki %*% t(kn1) 
  return(a / b - sn)
}

## ----rxeval---------------------------------------------------------
bestk <- which.min(apply(X0, 1, IMSPE.r, X0 = X0, theta = 0.25, r = r1))
Wijs <- Wij(X0, theta = 0.25, type = "Gaussian")
Ki <- solve(cov_gen(X0, theta = 0.25, type = "Gaussian") + diag(rn))
rx.thresh <- apply(XX, 1, rx, X0 = X0, rn = rn, theta = 0.25, Ki = Ki,
  kstar = bestk, Wijs = Wijs)

## ----threersfig, echo=FALSE, fig.height=5, fig.width=5, fig.show='hide'----
plot(X0, rn, xlab = "x", ylab = "r(x)", xlim = c(0, 1), ylim = c(2, 8),
  col = 2, main = "Two variance hypotheses")
lines(XX, r1(XX), col = 3)
lines(XX, r2(XX), col = 4)
lines(XX, rx.thresh, lty = 2, col = "darkgrey")
points(XX[xstar1], r1(XX[xstar1]), pch = 23, bg = 3)
points(XX[xstar2], r2(XX[xstar2]), pch = 23, bg = 4)
points(X0, rn, col = 2)

## ----threeimspefig, echo = FALSE, fig.height=5, fig.width=5, fig.show='hide'----
plot(XX, imspe1, type = "l", col = 3, ylab = "IMSPE", xlab = "x", 
  ylim = c(0.6, 0.7), main = "IMSPE for two variances")
lines(XX, imspe2, col = 4)
abline(v = X0, lty = 3, col = 'red')
points(XX[xstar1], imspe1[xstar1], pch = 23, bg = 3)
points(XX[xstar2], imspe2[xstar2], pch = 23, bg = 4)

## ----forr-----------------------------------------------------------
fn <- function(x) { 1/3 * (exp(sin(2 * pi * x))) }
fr <- function(x) { f1d2(x) + rnorm(length(x), sd = fn(x)) }

## ----forrinit-------------------------------------------------------
X <- seq(0, 1, length = 10)
Y <- fr(X)
mod <- mleHetGP(X = X, Z = Y, lower = 0.0001, upper = 1)

## ----forrIMSPE------------------------------------------------------
opt <- IMSPE_optim(mod, h = 5)
c(X, opt$par)

## ----forrupdate-----------------------------------------------------
X <- c(X, opt$par)
Ynew <- fr(opt$par)
Y <- c(Y, Ynew) 
mod <- update(mod, Xnew = opt$par, Znew = Ynew, ginit = mod$g * 1.01)

## ----forr500--------------------------------------------------------
for(i in 1:489) {
  opt <- IMSPE_optim(mod, h = 5)
  X <- c(X, opt$par)
  Ynew <- fr(opt$par)
  Y <- c(Y, Ynew)
  mod <- update(mod, Xnew = opt$par, Znew = Ynew, ginit = mod$g * 1.01)
  if(i %% 25 == 0) { 
    mod2 <- mleHetGP(X = list(X0 = mod$X0, Z0 = mod$Z0, mult = mod$mult),
    Z = mod$Z, lower = 0.0001, upper = 1)
    if(mod2$ll > mod$ll) mod <- mod2
  } 
}

## ----forrn, echo=FALSE, results='hide'------------------------------
nrow(mod$X0)

## ----forrpred-------------------------------------------------------
xgrid <- seq(0, 1, length = 1000)
p <- predict(mod, matrix(xgrid, ncol = 1)) 
pvar <- p$sd2 + p$nugs

## ----forrfig, echo = FALSE, fig.height=5, fig.width=6, out.width="5in", out.height="4.2in", fig.align='center', fig.cap="\\label{fig:forr}Sequential design with horizon $h=5$.  The truth is in black and the predictive distribution in red."----
plot(xgrid, f1d2(xgrid), type = "l", xlab = "x", ylab = "y", 
  main="1d example, IMSPE h=5", ylim = c(-4, 5))
lines(xgrid, qnorm(0.05, f1d2(xgrid), fn(xgrid)), col = 1, lty = 2)
lines(xgrid, qnorm(0.95, f1d2(xgrid), fn(xgrid)), col = 1, lty = 2)
points(X, Y)
segments(mod$X0, rep(0, nrow(mod$X0)) - 4, mod$X0, mod$mult * 0.25 - 4, 
  col = "gray")
lines(xgrid, p$mean, col = 2)
lines(xgrid, qnorm(0.05, p$mean, sqrt(pvar)), col = 2, lty = 2)
lines(xgrid, qnorm(0.95, p$mean, sqrt(pvar)), col = 2, lty = 2)
legend("top", c("truth", "estimate"), col = 1:2, lty = 1:2)

## ----adapt, warning=FALSE,message=FALSE-----------------------------
X <- seq(0, 1, length = 10)
Y <- fr(X)
mod.a <- mleHetGP(X = X, Z = Y, lower = 0.0001, upper = 1)
h <- rep(NA, 500)

## ----adapt2---------------------------------------------------------
for(i in 1:490) {
  h[i] <- horizon(mod.a)
  opt <- IMSPE_optim(mod.a, h = h[i])
  X <- c(X, opt$par)
  Ynew <- fr(opt$par)
  Y <- c(Y, Ynew)
  mod.a <- update(mod.a, Xnew = opt$par, Znew = Ynew, ginit = mod.a$g * 1.01)
  if(i %% 25 == 0) { 
    mod2 <- mleHetGP(X = list(X0 = mod.a$X0, Z0 = mod.a$Z0,
      mult = mod.a$mult), Z = mod.a$Z, lower = 0.0001, upper = 1)
    if(mod2$ll > mod.a$ll) mod.a <- mod2
  } 
}

## ----adapt3, echo = FALSE-------------------------------------------
p.a <- predict(mod.a, matrix(xgrid, ncol = 1))
pvar.a <- p.a$sd2 + p.a$nugs

## ----adapfig, echo = FALSE, fig.height=4, fig.width=8, out.width="6in", out.height="3in", fig.align='center', fig.cap="\\label{fig:adapt}{\\em Left:} Horizons chosen per iteration; {\\em right:} final design and predictions versus the truth, similar to Figure \\ref{fig:forr}."----
par(mfrow = c(1, 2))
plot(h, main = "Horizon", xlab = "Iteration")
plot(xgrid, f1d2(xgrid), type = "l", xlab = "x", ylab = "y",
  main = "Adaptive Horizon Design", ylim = c(-4, 5))
lines(xgrid, qnorm(0.05, f1d2(xgrid), fn(xgrid)), col = 1, lty = 2)
lines(xgrid, qnorm(0.95, f1d2(xgrid), fn(xgrid)), col = 1, lty = 2)
points(X, Y)
segments(mod$X0, rep(0, nrow(mod$X0)) - 4, mod$X0, mod$mult * 0.25 - 4, 
  col = "gray")
lines(xgrid, p$mean, col = 2)
lines(xgrid, qnorm(0.05, p$mean, sqrt(pvar.a)), col = 2, lty = 2)
lines(xgrid, qnorm(0.95, p$mean, sqrt(pvar.a)), col = 2, lty = 2)

## ----adaptn, echo=FALSE, results='hide'-----------------------------
nrow(mod.a$X0)

## ----rmsescore------------------------------------------------------
ytrue <- f1d2(xgrid)
yy <- fr(xgrid)
rbind(rmse = c(h5 = mean((ytrue - p$mean)^2), 
  ha = mean((ytrue - p.a$mean)^2)), 
  score = c(h5 = - mean((yy - p$mean)^2 / pvar + log(pvar)), 
  ha = -mean((yy - p.a$mean)^2 / pvar.a + log(pvar.a))))

## ----atoatime-------------------------------------------------------
c(n = nrow(out.a$X0), N = length(out.a$Z), time = out.a$time)

## ----atoatestscore--------------------------------------------------
sc.a <- scores(out.a, Xtest = Xtest, Ztest = Ztest)
c(batch = sc, adaptive = sc.a)

## ----atorebuild-----------------------------------------------------
out.a <- rebuild(out.a)

## ----atoadapt-------------------------------------------------------
Wijs <- Wij(out.a$X0, theta = out.a$theta, type = out.a$covtype)
h <- horizon(out.a, Wijs = Wijs)
control <- list(tol_dist = 1e-4, tol_diff = 1e-4, multi.start = 30)
opt <- IMSPE_optim(out.a, h, Wijs = Wijs, control = control)

## ----atoopt---------------------------------------------------------
opt$par

## ----atooptunique---------------------------------------------------
opt$path[[1]]$new

## ----EIahead, warning=FALSE, message=FALSE--------------------------
X <- seq(0, 1, length = 10)
X <- c(X, X, X)
Y <- -fr(X)
mod <- mleHetGP(X = X, Z = Y)

## ----EIahead2-------------------------------------------------------
library("parallel")
ncores <- 1 # or: detectCores()
for(i in 1:470) {
  opt <- crit_optim(mod, crit = "crit_EI", h = 5, ncores = ncores)
  X <- c(X, opt$par)
  Ynew <- -fr(opt$par)
  Y <- c(Y, Ynew)
  mod <- update(mod, Xnew = opt$par, Znew = Ynew, ginit = mod$g * 1.01)
  if(i %% 25 == 0) {
    mod2 <- mleHetGP(X = list(X0 = mod$X0, Z0 = mod$Z0, mult = mod$mult),
      Z = mod$Z, lower = 0.0001, upper = 1)
    if(mod2$ll > mod$ll) mod <- mod2
  }
}

## ----EIahead3-------------------------------------------------------
p <- predict(mod, matrix(xgrid, ncol = 1))
pvar <- p$sd2 + p$nugs

## ----EIgraphs, echo = FALSE,fig.height=5, fig.width=6, out.width="5in", out.height="4.2in", fig.align='center', fig.cap="\\label{fig:ei} Sequential optimization with horizon $h = 5$. The truth is in black and the predictive distribution in red ."----
plot(xgrid, f1d2(xgrid), type = "l", xlab = "x", ylab = "y",
  ylim = c(-4, 5), main = "1d example with EI, h = 5")
lines(xgrid, qnorm(0.05, f1d2(xgrid), fn(xgrid)), col = 1, lty = 2)
lines(xgrid, qnorm(0.95, f1d2(xgrid), fn(xgrid)), col = 1, lty = 2)
points(X, -Y)
segments(mod$X0, rep(0, nrow(mod$X0)) - 4, mod$X0, mod$mult * 0.5 - 4,
  col = "gray")
lines(xgrid, -p$mean, col = 2)
lines(xgrid, qnorm(0.05, -p$mean, sqrt(pvar)), col = 2, lty = 2)
lines(xgrid, qnorm(0.95, -p$mean, sqrt(pvar)), col = 2, lty = 2)
legend("top", c("truth", "estimate"), col = 1:2, lty = 1:2)

## ----EIreps---------------------------------------------------------
nrow(mod$X0)

## ----Contour_ahead, warning=FALSE, message=FALSE--------------------
X <- seq(0, 1, length = 10)
X <- c(X, X, X)
Y <- fr(X)
mod <- mleHetGP(X = X, Z = Y)

for(i in 1:470) {
  opt <- crit_optim(mod, crit = "crit_cSUR", h = 5, ncores = ncores)
  X <- c(X, opt$par)
  Ynew <- fr(opt$par)
  Y <- c(Y, Ynew)
  mod <- update(mod, Xnew = opt$par, Znew = Ynew, ginit = mod$g * 1.01)
  if(i %% 25 == 0) {
    mod2 <- mleHetGP(X = list(X0 = mod$X0, Z0 = mod$Z0, mult = mod$mult),
      Z = mod$Z, lower = 0.0001, upper = 1)
    if(mod2$ll > mod$ll) mod <- mod2
  }
}

p <- predict(mod, matrix(xgrid, ncol = 1))
pvar <- p$sd2 + p$nugs

## ----contour_n------------------------------------------------------
nrow(mod$X0)

## ----cSURgraphs, echo = FALSE, fig.height=5, fig.width=6, out.width="5in", out.height="4.2in", fig.align='center', fig.cap="\\label{fig:contour} Sequential contour finding with horizon $h = 5$. The truth is in black and the predictive distribution in red."----
plot(xgrid, f1d2(xgrid), type = "l", xlab = "x", ylab = "y", 
  ylim = c(-4, 5), main="1d example with cSUR, h = 5")
lines(xgrid, qnorm(0.05, f1d2(xgrid), fn(xgrid)), col = 1, lty = 2)
lines(xgrid, qnorm(0.95, f1d2(xgrid), fn(xgrid)), col = 1, lty = 2)
points(X, Y)
segments(mod$X0, rep(0, nrow(mod$X0)) - 4, mod$X0, mod$mult * 0.5 - 4,
  col="gray")
lines(xgrid, p$mean, col = 2)
lines(xgrid, qnorm(0.05, p$mean, sqrt(pvar)), col = 2, lty = 2)
lines(xgrid, qnorm(0.95, p$mean, sqrt(pvar)), col = 2, lty = 2)
legend("top", c("truth", "estimate"), col = 1:2, lty = 1:2)
abline(h = 0, col = "blue")

