# Naive mean field inference
naivemf <- function(xlevels, y, weights, data, prior.count=0){

   groups <- levels(factor(y))
   h <- vector('list',length(groups))
   names(h) <- groups
   h0 <- vector('list', length(xlevels))
   names(h0) <- names(xlevels)
   for(iy in groups){
     h[[iy]] <- vector('list',length(xlevels))
     names(h[[iy]]) <- names(xlevels)
   }
   if(is.null(weights)) weights <- rep(1, NROW(data))

   dev <- rep(0,length(xlevels))
   names(dev) <- names(xlevels)
   ntot <- sum(weights)
   for(i in seq_along(xlevels)){
     Li <- length(xlevels[[i]])
     f0 <- rep(0, Li)  # pooled inference
     names(f0) <- xlevels[[i]]
     for(w in xlevels[[i]])
       f0[as.character(w)] <- sum((data[,names(xlevels)[i]]==w)*weights)
     f0 <- f0 + prior.count/Li
     f0 <- f0 /(ntot + prior.count)
     h0[[i]] <- log(f0[-1]/f0[1])
  
     L <- 0
     for(iy in groups){
       ny <- sum(weights[y==iy])
       fv <- rep(0,Li)
       names(fv) <- xlevels[[i]]
       for(w in xlevels[[i]])
         fv[as.character(w)] <- 
          sum((data[y==iy,names(xlevels)[i]]==w)*weights[y==iy])
       fv <- fv + prior.count/Li
       fv <- fv /(ny+prior.count)
       h[[iy]][[i]] <- log(fv[-1]/fv[1])
       names(h[[iy]][[i]]) <- xlevels[[i]][-1]
       L <- L+ ny*sum(fv*log(fv))
     }
     dev[i] <- 2*(L - ntot*sum(f0*log(f0)))
   }

   ans <- list(h=h, h0=h0, dev=dev)
   return(ans)
}