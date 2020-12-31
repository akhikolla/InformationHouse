sigma12i <- function(sigma) {
    temp <- eigen(sigma, symmetric = T)
    sigma12i <- temp$vectors %*% diag(1/temp$values^0.5) %*% 
        t(temp$vectors)
    sigma12i
}
