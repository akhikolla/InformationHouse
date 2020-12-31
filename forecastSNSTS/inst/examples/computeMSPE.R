T <- 50
m <- 10

P <- 5
H <- 2

N_min <- 20

X <- rnorm(T)


mspe <- MSPE(X, m1 = T - m + 1, m2 = T, P = P, H = H, N = c(0, N_min:(T-m-H)))
N <- mspe$N
M <- mspe$mspe

h <- 1

plot(mspe, h, N_min = N_min, legend = (h == 1))

## Find minimums
idx1_s <- which(M[h, , N == 0] == min(M[h, , N == 0]), arr.ind = TRUE)[1]
abline(h = M[h, idx1_s, N == 0], col = idx1_s, lty="dashed", lwd = 2)

for (p in 1:p_max) {
  idx1_ls <- which(M[h, , N != 0] == min(M[h, , N != 0]), arr.ind = TRUE)[1,]
  idx1_ls_p <- which(M[h, p, N != 0] == min(M[h, p, N != 0]), arr.ind = TRUE)[1]
  abline(v = N[idx1_ls_p], col = p, lty = "dotted")
}
