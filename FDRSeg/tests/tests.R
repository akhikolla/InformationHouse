library(FDRSeg)
library(stepR)

## independent Gaussian data

alpha <- 0.05 
r     <- 1e3   # number of Monte-Carlo simulations
n     <- 500   # number of observations
K     <- 50    # number of true change-points
sd    <- 0.3   # standard deviation of random noise

# generate teeth signal
set.seed(2)
u0 <- teethfun(n, K)
Y  <- rnorm(n, mean = u0, sd = sd)

# simulate quantiles
qs  <- simulQuantile(1-alpha, n, r, "smuce")
qfs <- simulQuantile(1-alpha, n, r, "fdrseg")

# change-point regression
us  <- smuce(Y, q = qs, sd = sd)
ufs <- fdrseg(Y, q = qfs, sd = sd)

# plot result
plot(Y, pch = 20, col = "grey", xlab = "", ylab = "", main = expression(alpha*" = 0.05"))
lines(u0, type = "s")
lines(evalStepFun(us),  type = "s", col = "blue")
lines(evalStepFun(ufs), type = "s", col = "red")
legend("topleft", c("Truth", "SMUCE", "FDRSeg"), lty = c(1, 1, 1), col = c("black", "blue", "red"))


## dependent Gaussian data
# simulate data (a continuous time Markov chain)
ts        <- 0.2  # sampling time
SNR       <- 3    # signal-to-noise ratio
sampling  <- 1e4  # sampling rate 10 kHz
over      <- 10   # tenfold oversampling
cutoff    <- 1e3  # 1 kHz 4-pole Bessel-filter, adjusted for oversampling
simdf     <- dfilter("bessel", list(pole=4, cutoff=cutoff/sampling/over))
transRate <- 50
rates     <- rbind(c(0, transRate), c(transRate, 0))
set.seed(123)
sim <- contMC(ts*sampling, c(0,SNR), rates, sampling = sampling, family = "gaussKern",  param = list(df=simdf, over=over, sd=1))
Y   <- sim$data$y
x   <- sim$data$x

# D-FDRseg 
convKern <- dfilter("bessel", list(pole=4, cutoff=cutoff/sampling))$kern
alpha    <- 0.1
qdfs     <- simulQuantile(1 - alpha, ts*sampling, type = "dfdrseg", convKern = convKern)
udfs     <- dfdrseg(Y, qdfs, convKern = convKern)

# plot results
plot(x, Y, pch = 20, col = "grey", xlab="", ylab = "", main = "Simulate Ion Channel Data")
lines(sim$discr, col = "blue")
lines(x, evalStepFun(udfs), col = "red")
legend("topleft", c("Truth", "D-FDRSeg"), lty = c(1, 1), col = c("blue", "red"))

