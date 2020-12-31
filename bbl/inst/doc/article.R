### R code from vignette source 'article.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: data
###################################################
titanic <- as.data.frame(Titanic)
titanic
freq <- titanic$Freq
titanic <- titanic[, 1:4]


###################################################
### code chunk number 3: raw
###################################################
library('bbl')
titanic_raw <- freq2raw(data = titanic, freq = freq)
head(titanic_raw)
summary(titanic_raw)


###################################################
### code chunk number 4: lr
###################################################
gfit0 <- glm(Survived ~ Class + Sex + Age, family = binomial(), 
  data = titanic, weights=freq)
gfit0
summary(gfit0)


###################################################
### code chunk number 5: glm2
###################################################
gfit1 <- glm(Survived ~ (Class + Sex + Age)^2, family = binomial(), 
  data = titanic, weights = freq)
summary(gfit1)


###################################################
### code chunk number 6: div
###################################################
set.seed(159)
nsample <- NROW(titanic_raw)
flag <- rep(TRUE, nsample)
flag[sample(nsample, nsample/2)] <- FALSE
dtrain <- titanic_raw[flag,]
dtest <- titanic_raw[!flag,]


###################################################
### code chunk number 7: lr
###################################################
gfit2 <- glm(Survived ~ Class * Sex + Sex * Age, family = binomial(), 
  data = dtrain)
prl <- predict(gfit2, newdata = dtest)
yhat <- ifelse(prl > 0, 'Yes', 'No')
mean(yhat == dtest$Survived)
gauc <- pROC::roc(response = dtest$Survived, predictor = prl, 
  direction = '<')$auc
gauc


###################################################
### code chunk number 8: glmnet
###################################################
if(!require('glmnet'))
  install.packages('glmnet')
library('glmnet')
xdat <- data.matrix(dtrain[, 1:3])
y <- dtrain[, 4]
gnet <- cv.glmnet(x = xdat, y = y, family = 'binomial', alpha = 1,
  nfolds = 5, type.measure = 'auc')
plot(gnet)


###################################################
### code chunk number 9: glmnet
###################################################
plot(gnet)


###################################################
### code chunk number 10: glmnet2
###################################################
head(xdat)


###################################################
### code chunk number 11: class
###################################################
bfit0 <- bbl(Survived ~ Class + Sex + Age, data = titanic, weights = freq,
  prior.count = 0)


###################################################
### code chunk number 12: print
###################################################
bfit0


###################################################
### code chunk number 13: summary
###################################################
summary(bfit0)


###################################################
### code chunk number 14: survival
###################################################
bfit <- bbl(Survived ~ Class * Sex + Sex * Age, data = titanic, 
  weights = freq)
bfit


###################################################
### code chunk number 15: plot
###################################################
oldpar <- par(mar = c(6, 4.5, 3, 4), tck = -0.05, cex.axis = 0.8)
plot(bfit)
par(oldpar)


###################################################
### code chunk number 16: predict.bbl
###################################################
bfit2 <- bbl(Survived ~ Class * Sex + Sex * Age, data = dtrain)
pr <- predict(bfit2, newdata = dtest, type = 'prob')
head(pr)
auc <- pROC::roc(response = dtest$Survived, predictor = pr[, 2], 
  direction = '<')$auc
auc


###################################################
### code chunk number 17: cvsim
###################################################
cv <- crossVal(Survived ~ .^2, data = dtrain, method = 'pseudo', 
  lambda = 10^seq(-5, -2, 0.2), verbose = 0)
cv
plot(cv, mar=c(4, 4, 3, 3), tck = -0.04, bty = 'n')


###################################################
### code chunk number 18: plotcv
###################################################
plot(cv, mar = c(4, 4, 3, 3), tck = -0.04, bty = 'n')


###################################################
### code chunk number 19: pr2
###################################################
model <- bbl(Survived ~ .^2, data = dtrain, lambda = cv$regstar)
pr2 <- predict(model, newdata = dtest)
bscore <- mean(dtest$Survived == pr2$yhat)
bscore
bauc <- pROC::roc(response = dtest$Survived, predictor = pr2[,2], 
  direction = '<')$auc
bauc


###################################################
### code chunk number 20: mcmc
###################################################
map <- mcSample(bfit, nstep = 1000, progress.bar = FALSE)
map


###################################################
### code chunk number 21: sim1
###################################################
predictors <- list()
m <- 5
L <- 3
for(i in 1:m) predictors[[i]] <- seq(0, L-1)
par <- randompar(predictors)
names(par)


###################################################
### code chunk number 22: sample
###################################################
xi <- sample_xi(nsample = 10000, predictors = predictors, h = par$h, 
  J = par$J, code_out = TRUE)
head(xi)


###################################################
### code chunk number 23: mle
###################################################
fit <- mlestimate(xi = xi, method = 'pseudo', lambda = 0)


###################################################
### code chunk number 24: par
###################################################
oldpar <- par(mar = c(4, 4, 1, 2), lwd = 0.5, cex.axis = 0.8, 
  cex.lab = 1.0, mgp = c(2.2 ,0.9, 0), tck = -0.03)
range <- range(par$h, par$J, fit$h, fit$J)
plot(x = unlist(par$h), y = unlist(fit$h), bg = 'cornflowerblue', 
  xlim = range, ylim = range, pch = 21, cex = 0.8, xlab = 'True', 
  ylab = 'Inferred', lwd = 0.7, xaxt = 'n', yaxt = 'n', bty = 'n')
axis(side = 1, at = seq(-1.5, 1.5, 0.5), lwd = 0.5, las = 1)
axis(side = 2, at = seq(-1.5, 1.5, 0.5), lwd = 0.5, las = 1)
segments(x0 = -1, x1 = 1, y0 = -1, y1 = 1, lty = 2, lwd = 0.7)
points(x = unlist(par$J), y = unlist(fit$J), pch = 24, bg = 'orange', 
  cex = 0.8, lwd = 0.7)
legend(x = 0.5, y = -0.5, legend = expression(italic(h), italic(J)), 
  cex = 0.8, pch = c(21, 24), pt.bg = c('cornflowerblue', 'orange'))
par(oldpar)


###################################################
### code chunk number 25: atgc
###################################################
nt <- c('a', 'c', 'g', 't')
set.seed(135)
for(i in 1:m) predictors[[i]] <- nt
names(predictors) <- paste0('v', 1:m)
par <- list()
par[[1]] <- randompar(predictors)
par[[2]] <- randompar(predictors, h0 = 0.1, J0 = 0.1)
dat <- randomsamp(predictors, response = c('ctrl', 'case'), par = par,
  nsample = 1000)


###################################################
### code chunk number 26: cr-mf
###################################################
cv <- crossVal(y ~ .^2, data = dat, method = 'mf', eps = seq(0, 1, 0.1),
  verbose=0)
cv


###################################################
### code chunk number 27: mf-par
###################################################
fit <- list()
eps <- c(0.2, 0.7, 1.0)
for(i in seq_along(eps))
  fit[[i]] <- bbl(y ~ .^2, data = dat, method = 'mf', eps = eps[i], 
    verbose = 0)


###################################################
### code chunk number 28: cv
###################################################
oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), lwd = 0.5, 
  cex.axis = 0.8, cex.lab = 0.9, mgp = c(2.2, 0.8, 0), tck = -0.03, 
  las = 1)
estar <- cv$regstar
plot(cv, xlab = expression(epsilon), ylab = 'AUC', lwd = 0.7, cex = 0.7,
  bty = 'n', log = '')
segments(x0 = estar, x1 = estar, y0 = 0, y1 = cv$maxscore, lty = 2, 
  lwd = 0.5, col = 'red')
title(adj = 0, cex.main = 1.2, font = 2, main = 'a')

for(i in 1:3){
  plot(x = c(unlist(par[[1]]$h), unlist(par[[2]]$h)), 
    y = unlist(coef(fit[[i]])$h), bg = 'cornflowerblue', 
    xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), pch = 21, cex = 0.7, 
    xlab = 'True', ylab = 'Inferred', lwd = 0.7, xaxt = 'n', yaxt = 'n',
    bty = 'n')
  axis(side = 1, at = seq(-1.5, 1.5, 0.5), lwd = 0.5, las = 1)
  axis(side = 2, at = seq(-1.5, 1.5, 0.5), lwd = 0.5, las = 1)
  segments(x0 = -2, x1 =2 , y0 = -2, y1 = 2, lty = 2, lwd = 0.7)
  points(x = c(unlist(par[[1]]$J), unlist(par[[2]]$J)), 
    y = unlist(coef(fit[[i]])$J), pch = 24, bg = 'orange', cex = 0.7, 
    lwd = 0.7)
  if(i==1) legend(x = 0.5, y = -0.5, legend = expression(italic(h), 
    italic(J)), cex = 0.8, pch = c(21, 24), 
    pt.bg = c('cornflowerblue', 'orange'))
  title(adj = 0,main = letters[i + 1], cex.main = 1.1, font = 2)
  mtext(side = 3, line = 1.0, cex = 0.8, bquote(epsilon == .(eps[i])),
    adj = 0.5)
}
par(oldpar)


###################################################
### code chunk number 29: ntaa
###################################################
set.seed(351)
n <- 2000
dat <- data.frame(b1 = sample(nt, size = n, replace = TRUE),
  b2 = sample(nt, size = n, replace = TRUE),
  b3 = sample(nt, size = n, replace = TRUE))
head(dat)


###################################################
### code chunk number 30: biostrings
###################################################
if(!require('Biostrings')){
  if(!require('BiocManager'))
    install.packages('BiocManager')
  BiocManager::install('Biostrings')
}
aa <- Biostrings::DNAString(paste(t(dat), collapse = ''))
aa
aa <- strsplit(as.character(Biostrings::translate(aa)), split = '')[[1]]
xdat <- cbind(data.frame(aa = aa), dat)
head(xdat)


###################################################
### code chunk number 31: aacv
###################################################
cv <- crossVal(aa ~ .^2, data = xdat, lambda = 10^seq(-3, 1, 0.5), 
  verbose = 0)
cv


###################################################
### code chunk number 32: codon
###################################################
panel <- expand.grid(b1 = nt, b2 = nt, b3 = nt)
head(panel)
dim(panel)
p <- predict(cv, panel)
ap <- Biostrings::DNAString(paste(t(panel), collapse = ''))
ap <- strsplit(as.character(Biostrings::translate(ap)), split = '')[[1]]
score <- mean(ap == p$yhat)
score


###################################################
### code chunk number 33: mnist
###################################################
dat0 <- read.csv(system.file('extdata/mnist_train.csv', package = 'bbl'))
dat <- removeConst(dat0)
dat[1:5, 1:10]


###################################################
### code chunk number 34: mnist_cval (eval = FALSE)
###################################################
## cv <- crossVal(y ~ .^2, data = dat, method = 'mf', eps = 0.05)


###################################################
### code chunk number 35: mnist3 (eval = FALSE)
###################################################
## mnist <- bbl(y ~ .^2, data = dat, method = 'mf', eps = 0.05)
## dtest <- read.csv(system.file('extdata/mnist_test.csv', package = 'bbl'))
## dtest <- dtest[, colnames(dtest) %in% colnames(dat)]
## pr <- predict(mnist, newdata = dtest[, -1], progress.bar = TRUE)
## score <- mean(pr$yhat == dtest$y)
## score


###################################################
### code chunk number 36: mnist_map1 (eval = FALSE)
###################################################
## mnist_map <- mcSample(mnist, nstep = 20, progress.bar = TRUE)
## oldpar <- par(mfrow = c(2, 5), mar = c(1, 1, 1, 1))
## xvar <- colnames(dat0[, -1])
## xmap <- apply(mnist_map$xmax, 1:2, as.numeric)
## xf <- matrix(0, nrow = length(xvar), ncol = 10)
## rownames(xf) <- xvar
## for(i in 1:10) xf[rownames(xmap), i] <- xmap[, i]
## for(i in 1:10){
##   mat <- matrix(t(xf[, i]), nrow = 28, ncol = 28)
##   image(x = 1:28, y = 1:28, z = mat[, 28:1], col = c('white', 'black'), 
##     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
## }
## par(oldpar)


###################################################
### code chunk number 37: jaspar
###################################################
seq <- readFasta(system.file('extdata/MA0014.3.fasta', package = 'bbl'))
head(seq)
dim(seq)


###################################################
### code chunk number 38: jaspar2
###################################################
set.seed(561)
nsample <- NROW(seq)
m <- NCOL(seq)
nt <- c('A', 'C', 'G', 'T')
ctrl <- as.matrix(seq)
for(k in seq_len(nsample))
  ctrl[k, sample(m, 3)] <- sample(nt, 3, replace = TRUE)
colnames(ctrl) <- 1:m
data <- rbind(data.frame(y = rep('Binding', nsample), seq), 
  data.frame(y = rep('Non-binding', nsample), ctrl))
data <- data[sample(NROW(data)), ]


###################################################
### code chunk number 39: jaspar3
###################################################
ps <- crossVal(y ~ .^2, data = data, method = 'pseudo', 
  lambda = 10^seq(-2, -1, 0.2), verbose = 0)
ps
mf <- crossVal(y ~ .^2, data = data, method = 'mf', 
  eps = seq(0.1, 0.4, 0.1), verbose = 0)
mf


