fun_arglist <- function(lowerIC, upperIC, X, trunc, normalize, tol, niter) {

   l <- as.numeric(lowerIC)
   r <- as.numeric(upperIC)
   if(!is.null(trunc)) trunc <- as.numeric(trunc)
   n <- length(l)
   lr <- cbind(l, r)
   X <- as.matrix(X)

   if (normalize == TRUE) {
     z <- apply(X,2, function(x) (x-mean(x))/sqrt(sum((x-mean(x))^2)/n))
     true_sd <- sqrt(apply(X,2,var)*(n-1)/n)
   } else {
     z <- X
     true_sd <- 1
   }

   if (is.null(trunc)) {
     olr <- order(c(l + 1e-10, r))
     int0 <- rbind(cbind(0, l), cbind(1, r))[olr,]
   } else {
     olr <- order(c(l + 1e-10, r, trunc - 1e-10))
     int0 <- rbind(cbind(0, l), cbind(1, r), cbind(1, trunc))[olr,]
   }

   int <- cbind(cbind(int0[-nrow(int0), 1], c(int0[-1, 1])), cbind(int0[-nrow(int0), 2], c(int0[-1, 2])))
   set0 <- int[(int[, 1]==0)&(int[, 2]==1),c(3,4)]
   set <- set0[!is.infinite(rowSums(set0)), ]

   args <- list()
   args$l <- l
   args$r <- r
   args$trunc <- trunc
   args$n <- n
   args$z <- z
   args$true_sd <- true_sd
   args$set <- set
   args$tol <- tol
   args$niter <- niter

   return(args)

}


