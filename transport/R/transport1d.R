wasserstein1d <- function(a, b, p=1, wa=NULL, wb=NULL) {  	
  m <- length(a)
  n <- length(b)
  stopifnot(m > 0 && n > 0)
  if (m == n && is.null(wa) && is.null(wb)) {
  	return(mean(abs(sort(b)-sort(a))^p)^(1/p))
  }
  if (is.null(wa)) {wa <- rep(1,m)}
  if (is.null(wb)) {wb <- rep(1,n)}
  stopifnot(length(wa) == m && length(wb) == n)
  ua <- (wa/sum(wa))[-m]
  ub <- (wb/sum(wb))[-n]
  cua <- c(cumsum(ua))  
  cub <- c(cumsum(ub))  
  temp <- cut(cub,breaks=c(-Inf,cua,Inf))
  arep <- table(temp) + 1  
  temp <- cut(cua,breaks=c(-Inf,cub,Inf))
  brep <- table(temp) + 1
  # we sum over rectangles with cuts on the vertical axis each time one of the two ecdfs makes a jump
  # xrep and yrep tell us how many times each of the x and y data have to be repeated in order to get the points on the horizontal axis
  # note that sum(xrep)+sum(yrep) = m+n-1 (we do not count the height-zero final rectangle where both ecdfs jump to 1)

  aa <- rep(sort(a), times=arep)
  bb <- rep(sort(b), times=brep)

  uu <- sort(c(cua,cub))
  uu0 <- c(0,uu)
  uu1 <- c(uu,1)
  areap <- sum((uu1-uu0)*abs(bb-aa)^p)^(1/p)
#  print(rbind(uu1-uu0, pmax(aa,bb)-pmin(aa,bb)))
  return(areap)
}


# plots without weights
# 
# plot(ecdf(a),verticals=TRUE, xlim=c(-2,2))
# lines(ecdf(b),col=2,verticals=TRUE)


# plots with weights
#
# ua <- (wa/sum(wa))[-length(wa)]
# ub <- (wb/sum(wb))[-length(wb)]
# cua <- c(cumsum(ua))  
# cub <- c(cumsum(ub))  
# plot(sort(c(-2,a,2)),c(0,cua,1,1),type="s",xlim=c(-2,2))
# lines(sort(c(-2,b,2)),c(0,cub,1,1),type="s",col=2)
