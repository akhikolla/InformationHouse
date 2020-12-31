set.seed(3)
data = data.frame(y = rnorm(10), x1 = runif(10),
                  x2 = runif(10))
d = as.matrix(dist(data[,c("x1", "x2")]))
x = cbind(1, data$x1)
n = nrow(d)
v = exp(-d) + diag(n)
coeff = solve(crossprod(x, solve(v, x)), crossprod(x, solve(v, data$y)))
m = x %*% coeff
ll_dmvnorm = mvtnorm::dmvnorm(data$y, mean = m,
                              sigma = v, log = TRUE)
m2 = rep(1, n)
ll_dmvnorm2 = mvtnorm::dmvnorm(data$y, mean = m2,
                               sigma = v, log = TRUE)

# save output
fpath = system.file("testdata",  package = "gear")
fname = paste(fpath, "/loglik_cmodMan_dmvtnorm.rda", sep = "")
save(data, v, ll_dmvnorm, ll_dmvnorm2,
     compress = "bzip2",
     file = fname,
     version = 2)
