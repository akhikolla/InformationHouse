wilcoxontau <- function(resd, p, delta = if ((length(resd)/p) > 5) 0.8 else 0.95, 
            param = 2, verbose=FALSE) {
    # Bootstrapping for large datasets
    if(length(resd) > 3000) {
      sample(resd, 1000, replace=TRUE)
      
      wilcoxontau = mean(sapply(1:100, function(x) {
        if(verbose == TRUE) { 
          cat(paste("wilcoxontau: bootstrap iteration", x, "\n"))
        }
        return(wilcoxontau(sample(resd, 1000, replace=TRUE), p, delta, param))
      }))
      
      return(wilcoxontau)
    }
    
    eps <- 1e-06
    n <- length(resd)
    temp <- pairup(resd)
    
    dresd = abs(temp[,1] - temp[,2])
    
    # this is a corrective measure
    dresd = remove_k_smallest(dresd, p)
    
    
    tdeltan <- quantile(dresd, delta)/sqrt(n)
    
    w <- rep(0, length(dresd))
    w[dresd <= tdeltan] <- 1
    cn <- 2/(n * (n - 1))
    scores = sqrt(12) * ((1:n)/(n + 1) - 0.5)
    mn = mean(scores)
    con = sqrt(sum((scores - mn)^2)/(n + 1))
    scores = (scores - mn)/con
    
    dn = max(scores) - min(scores)
    #dn = scores[n] - scores[1]
    
    wilcoxontau <- sqrt(n/(n - p - 1)) * ((2 * tdeltan)/(dn * sum(w) * cn))
    
    w <- rep(0, n)
    stan <- (resd - median(resd))/mad(resd)
    w[abs(stan) < param] <- 1
    hubcor <- sum(w)/n
    
    if (hubcor < eps) {
      hubcor <- eps
    }
    
    fincor <- 1 + (((p + 1)/n) * ((1 - hubcor)/hubcor))
    wilcoxontau <- fincor * wilcoxontau
    names(wilcoxontau) <- NULL
    wilcoxontau
  }

old.wilcoxontau <- function(resd, p, delta = if ((length(resd)/p) > 5) 0.8 else 0.95, 
            param = 2) {
    eps <- 1e-06
    n <- length(resd)
    temp <- pairup(resd)
    
    #dresd = abs(temp[,1] - temp[,2])
    # dresd = remove (p+1) smallest values
    
    #for(i in 1:p) {
    #  dresd = dresd[-which.min(dresd)]
    #}
    
    dresd <- sort(abs(temp[, 1] - temp[, 2]), method='quick')
    dresd = dresd[(p + 1):choose(n, 2)]
    
    tdeltan <- quantile(dresd, delta)/sqrt(n)
    
    w <- rep(0, length(dresd))
    w[dresd <= tdeltan] <- 1
    cn <- 2/(n * (n - 1))
    scores = sqrt(12) * ((1:n)/(n + 1) - 0.5)
    mn = mean(scores)
    con = sqrt(sum((scores - mn)^2)/(n + 1))
    scores = (scores - mn)/con
    
    #dn = max(scores) - min(scores)
    dn = scores[n] - scores[1]
    
    wilcoxontau <- sqrt(n/(n - p - 1)) * ((2 * tdeltan)/(dn * sum(w) * cn))
    
    w <- rep(0, n)
    stan <- (resd - median(resd))/mad(resd)
    w[abs(stan) < param] <- 1
    hubcor <- sum(w)/n
    
    if (hubcor < eps) {
      hubcor <- eps
    }
    
    fincor <- 1 + (((p + 1)/n) * ((1 - hubcor)/hubcor))
    wilcoxontau <- fincor * wilcoxontau
    names(wilcoxontau) <- NULL
    wilcoxontau
  }