factor.seg.alg <- function(x, r=NULL, bn.op=2, sig.lev=.05, max.q=NULL, q.seq=NULL,
                           qlen=5, qby=0, dw=NULL, p=NULL, B=200, scales=NULL,
                           rule=NULL, mby=NULL, tby=NULL, idio.diag=FALSE,
                           do.parallel=TRUE, no.proc=2){

	n <- nrow(x); T <- ncol(x)
	if(is.null(max.q)) max.q <- max(round(20, sqrt(min(n, T))))
	if(is.null(dw)) dw <- round(min(log(T)^2, T^(6/7)/4))
	if(do.parallel){
		cl <- parallel::makeCluster(no.proc); doParallel::registerDoParallel(cl)
	}

	gfm <- get.factor.model(x, bn.op = bn.op, max.q = max.q)
	if(is.null(r)) q.hat <- gfm$q.hat else q.hat <- r
	if(q.hat==0) q.hat <- 1
	if(is.null(q.seq)){
	  if(qlen > 0) q.seq <- sort(unique(c(q.hat, max.q,  round(seq(q.hat, max(min(n-1, 2*q.hat), max.q), length.out=qlen)))), decreasing=FALSE)
	  if(qlen==0 & qby > 0) q.seq <- sort(unique(c(q.hat, max.q,  round(seq(q.hat, max(min(n-1, 2*q.hat), max.q), by=qby)))), decreasing=FALSE)
    if(qlen==0 & qby==0) q.seq <- sort(unique(c(q.hat, max.q,  q.hat:max(min(n-1, 2*q.hat), max.q))), decreasing=FALSE)
	}

	est.cps <- cs.list <- list()
	for(qq in q.seq){
		cs <- common.seg(gfm, q=qq, scales=scales, sig.lev=sig.lev, rule=rule, B=B, p=p, dw=dw, mby=mby, tby=tby, do.parallel=do.parallel)
		est.cps <- c(est.cps, list(cs$est.cps))
		cs.list <- c(cs.list, list(cs))
	}
	qq <- max(which(unlist(lapply(est.cps, length))==max(unlist(lapply(est.cps, length)))))
	q <- q.seq[qq]
	cs <- cs.list[[qq]]
  common.est.cps <- est.cps[[qq]]

	is <- idio.seg(gfm, q=q, scales=scales, diag=idio.diag, sig.lev=sig.lev, rule=rule, B=B, p=p, dw=dw, mby=mby, tby=tby, do.parallel=do.parallel)
	idio.est.cps <- is$est.cps
	if(do.parallel) parallel::stopCluster(cl)

	return(list(cs.list=cs.list, r=q, common.est.cps=common.est.cps, idio.seg.res=is, idio.est.cps=idio.est.cps, gfm=gfm, q.seq=q.seq))

}

##########################################

get.factor.model <- function(x, max.q=NULL, q=NULL, bn=TRUE, bn.op=2, normalisation=TRUE){
	T <- ncol(x); n <- nrow(x)
	cnt <- min(n, T)
	if(is.null(max.q)) max.q <- round(sqrt(cnt))

	if(normalisation){
		mx <- matrix(rep(apply(x, 1, mean), each=T), byrow=TRUE, nrow=n)
		x <- x-mx
		sdx <- apply(x, 1, sd)
		x <- x/sdx
	} else{
		mx <- rep(0, n); sdx <- rep(1, n)
	}
	xx <- x%*%t(x)/T
	eig <- eigen(xx, symmetric=TRUE)
	lam <- eig$vectors[, 1:(cnt-1), drop=FALSE]*sqrt(n)
	f <- t(eig$vectors[, 1:(cnt-1), drop=FALSE])%*%x/sqrt(n)

	if(bn){
		ic <- rep(0, 1+max.q)
		ic[1] <- (bn.op <= 4)*log(mean(x^2)) + (bn.op==5)*mean(x^2)
		l <- 1
		while(l<=max.q){
			hchi <- lam[, 1:l, drop=FALSE]%*%f[1:l, , drop=FALSE]
			ic[l+1] <- (bn.op <= 4)*log(mean((x-hchi)^2)) +
				(bn.op==1)*l*(n+T)/(n*T)*log(n*T/(n+T)) +
				(bn.op==2)*l*(n+T)/(n*T)*log(cnt) +
				(bn.op==3)*l*log(cnt)/cnt +
				(bn.op==4)*l*((n+T-l)*log(n*T)/(n*T) + (n+T)/(n*T)*log(cnt))/2 +
				(bn.op==5)*(mean((x-hchi)^2)+l*mean((x-hchi)^2)*(n+T-l)*log(n*T)/(n*T))
			l <- l+1
		}
		q.hat <- which(ic==min(ic))-1
	} else{
		ic <- rep(0, max.q)
		q.hat <- q
	}

	return(list(lam = lam, f = f, norm.x=x, q.hat=q.hat, max.q=max.q, ic=ic))
}

make.tree <- function(y, dw, rule=NULL){
	len <- ncol(y)
	if(is.null(rule)) rule <- round(log(len, 2)/2)
	tree <- list(matrix(0, 6, 1))
	mat <- c()

	fd <- func_dc(y)
	stat <- fd$res
	test.stat <- max(stat[-c((1:dw), (len-dw):len)])
	hat.chp <- min(which(stat==test.stat))

	tree[[1]][1, 1] <- 1
	tree[[1]][2, 1] <- 1
	tree[[1]][3, 1] <- hat.chp
	tree[[1]][4, 1] <- len
	tree[[1]][5, 1] <- 0
	tree[[1]][6, 1] <- test.stat
	mat <- cbind(mat, c(tree[[1]][-5, ], 1, 1))

	j <- 1
	while(length(tree)==j & j < rule){
		npc <- dim(tree[[j]])[2]
		if(sum(tree[[j]][4,]-tree[[j]][2,]-rep(4*dw, npc)>0)){
			ncc <- 0; i <- 1
			while(i <= npc){
				if(tree[[j]][3, i]-tree[[j]][2, i]+1>4*dw){
					s <- tree[[j]][2, i]; e <- tree[[j]][3, i]
					fd <- func_dc(y[, s:e])
					stat <- fd$res
					test.stat <- max(stat[-c((1:dw), (e-s+1-dw):(e-s+1))])
					hat.chp <- s+min(which(stat==test.stat))-1

					if(length(tree)==j) tree <- c(tree, list(matrix(0, 6, 0)))
					ncc <- ncc+1
					tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 6, 1)), 6, ncc)
					tree[[j+1]][1, ncc] <- 2*tree[[j]][1, i]-1
					tree[[j+1]][2, ncc] <- s
					tree[[j+1]][3, ncc] <- hat.chp
					tree[[j+1]][4, ncc] <- e
					tree[[j+1]][5, ncc] <- 0
					tree[[j+1]][6, ncc] <- test.stat
					mat <- cbind(mat, c(tree[[j+1]][-5, ncc], j+1, ncc))
				}
				if(tree[[j]][4, i]-tree[[j]][3, i]>4*dw){
					s <- tree[[j]][3, i]+1; e <- tree[[j]][4, i]
					fd <- func_dc(y[, s:e])
					stat <- fd$res
					test.stat <- max(stat[-c((1:dw), (e-s+1-dw):(e-s+1))])
					hat.chp <- s+min(which(stat==test.stat))-1

					if(length(tree)==j) tree <- c(tree, list(matrix(0, 6, 0)))
					ncc <- ncc+1
					tree[[j+1]] <- matrix(c(tree[[j+1]], matrix(0, 6, 1)), 6, ncc)
					tree[[j+1]][1, ncc] <- 2*tree[[j]][1, i]
					tree[[j+1]][2, ncc] <- s
					tree[[j+1]][3, ncc] <- hat.chp
					tree[[j+1]][4, ncc] <- e
					tree[[j+1]][5, ncc] <- 0
					tree[[j+1]][6, ncc] <- test.stat
					mat <- cbind(mat, c(tree[[j+1]][-5, ncc], j+1, ncc))
				}
				i <- i+1
			}
			j <- j+1
		} else{
			break
		}
	}
	list(tree=tree, mat=mat)
}

##########################################

common.seg <- function(gfm, q, scales=NULL, sig.lev=0.05, rule=NULL, B=200, p=NULL,
                       dw=NULL, mby=NULL, tby=NULL, do.parallel=TRUE){
	lam <- gfm$lam; f <- gfm$f; nx <- gfm$norm.x
	T <- dim(nx)[2]
	if(is.null(p) || (p<=0 | p>=1)) p.seq <- apply(f[1:q,,drop=FALSE], 1, function(z){g <- get.gg(z); min(.5, ((g[2]/g[1])^2)^(-1/3)*T^(-1/5))}) else p.seq <- rep(p, q)
	if(is.null(scales)) scales <- -(1:floor(log(log(T, 2), 2)))

	y <- c()
	hat.chi <- lam[, 1:q, drop=FALSE]%*%f[1:q, , drop=FALSE]
	for(sc in scales){
	  cc <- func_coef(hat.chi, sc)
	  if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
	  chi.input <- t(func_input_on(cc))
	  y <- rbind(y, chi.input)
	}

	d <- nrow(y); len <- ncol(y)
	if(is.null(mby)) mby <- ceiling(log(d))
	if(is.null(tby)) tby <- ceiling(log(len))
	if(is.null(rule)) rule <- round(log(len, 2)/2)

	mat <- make.tree(y, dw, rule)$mat
	ns <- common.null.stat(mat=mat, p.seq=p.seq, B=B, q=q, lam=lam, f=f, scales=scales, mby=mby, tby=tby, do.parallel=do.parallel)

	k <- 1; pos.seq <- c()
	while(k <= ncol(mat)){
		pos <- mean(mat[5, k] > ns[k, ])
		pos.seq <- c(pos.seq, pos)
		if(pos < 1-sig.lev){
			l <- mat[6, k]+1; rm.ind <- mat[1, k]
			while(l <= max(mat[6, ])){
				ind <- which(mat[6, ]==l & is.element(mat[1, ], c(rm.ind*2-1, rm.ind*2)))
				if(length(ind) > 0){
					rm.ind <- mat[1, ind]
					mat <- mat[, -ind, drop=FALSE]
					l <- l+1
				} else{
					break
				}
			}
		}
		k <- k+1
	}
	tree <- list(matrix(0, 6, 1))
	mat <- rbind(mat, pos.seq); mat <- mat[, pos.seq >= 1-sig.lev, drop=FALSE]
	if(dim(mat)[2] > 0){
		for(l in 1:length(unique(mat[6, ]))){
			j <- unique(mat[6, ])[l]
			for(ncc in 1:sum(mat[6, ]==j)){
				k <- sum(mat[6, ] < j) + ncc
				if(length(tree)<j) tree <- c(tree, list(matrix(0, 6, 0)))
				tree[[l]] <- matrix(c(tree[[l]], matrix(0, 6, 1)), 6, ncc)
				tree[[l]][1, ncc] <- mat[1, k]
				tree[[l]][2, ncc] <- mat[2, k]
				tree[[l]][3, ncc] <- mat[3, k]
				tree[[l]][4, ncc] <- mat[4, k]
				tree[[l]][5, ncc] <- mat[8, k]
				tree[[l]][6, ncc] <- mat[5, k]
			}
		}
	}
	est.cps <- sort(mat[3, ])

	ls <- list(tree=tree, est.cps=est.cps+2^(-min(scales))-1, p=p.seq)
	return(ls)
}

common.null.stat <- function(mat, p.seq, B, q, lam, f, scales, mby, tby, do.parallel){
	len <- ncol(f)
	if(do.parallel){
		null.stat <- foreach::foreach(l=iterators::iter(1:B), .combine=cbind, .packages=c('Rcpp', 'RcppArmadillo', 'factorcpt')) %dopar% {
			boot.chi <- 0
			for(qq in 1:q){
				ind <- c()
				while(length(ind) < len){
					L <- rgeom(1, p.seq[qq]); I <- sample(len, 1); ind <- c(ind, rep(1:len, 1+ceiling(L/len))[I:(I+L-1)])
				}
				ind <- ind[1:len]
				boot.chi <- boot.chi + lam[, qq, drop=FALSE]%*%f[qq, ind, drop=FALSE]
			}
			by <- c()
			for(sc in scales){
			  cc <- func_coef(boot.chi, sc)
			  if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
			  by <- rbind(by, t(func_input_on(cc)))
			}
			tmp <- c()
			for(i in 1:ncol(mat)){
				s <- mat[2, i]; e <- mat[4, i]
				bfd <- func_dc_by(by[, s:e], mby, tby)
				tmp <- c(tmp, max(bfd$res))
				rm(bfd)
			}
			rm(by)
			tmp
	  	}
	} else{
		null.stat <- foreach::foreach(l=iterators::iter(1:B), .combine=cbind, .packages=c('Rcpp', 'RcppArmadillo', 'factorcpt')) %do% {
 			boot.chi <- 0
			for(qq in 1:q){
				ind <- c()
				while(length(ind)<len){
					L <- rgeom(1, p.seq[qq]); I <- sample(len, 1); ind <- c(ind, rep(1:len, 1+ceiling(L/len))[I:(I+L-1)])
				}
				ind <- ind[1:len]
				boot.chi <- boot.chi + lam[, qq, drop=FALSE]%*%f[qq, ind, drop=FALSE]
			}
			by <- c()
			for(sc in scales){
			  cc <- func_coef(boot.chi, sc)
			  if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
			  by <- rbind(by, t(func_input_on(cc)))
			}
			tmp <- c()
			for(i in 1:ncol(mat)){
				s <- mat[2, i]; e <- mat[4, i]
				bfd <- func_dc_by(by[, s:e], mby, tby)
				tmp <- c(tmp, max(bfd$res))
				rm(bfd)
			}
			rm(by)
			tmp
	  	}
	}
	null.stat
}

##########################################

idio.seg <- function(gfm, q, scales=NULL, diag=TRUE, sig.lev=0.05, rule=NULL, B=200,
                     p=NULL, dw=NULL, mby=NULL, tby=NULL, do.parallel=TRUE){
	lam <- gfm$lam[, 1:q, drop=FALSE]; f <- gfm$f[1:q, , drop=FALSE]; nx <- gfm$norm.x
	hat.vep <- nx - lam%*%f
	T <- dim(nx)[2]
	if(is.null(dw)) dw <- round(2*log(T))
	if(is.null(p) || (p<=0 | p>=1)) p <- min(.5, 1/mean(apply(hat.vep, 1, function(z){
		g <- get.gg(z); ((g[2]/g[1])^2)^(1/3)*T^(1/5)})))
	if(is.null(scales)) scales <- -(1:floor(log(log(T, 2), 2)))

	y <- NULL
	for(sc in scales){
		cc <- func_coef(hat.vep, sc)
		if(sc > min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
		if(diag){
			y <- rbind(y, t(func_input_on(cc)))
		} else{
			sgn <- sign(cc%*%t(cc))
			y <- rbind(y, t(func_input(cc, sgn)))
		}
	}
	d <- nrow(y); len <- ncol(y)
	if(is.null(mby)) mby <- ceiling(log(d))
	if(is.null(tby)) tby <- ceiling(log(len))
	if(is.null(rule)) rule <- round(log(len, 2)/2)

	mat <- make.tree(y, dw, rule)$mat
	ns <- idio.null.stat(mat=mat, p=p, B=B, hat.vep=hat.vep, diag=diag, scales=scales, mby=mby, tby=tby, do.parallel=do.parallel)

	k <- 1; pos.seq <- c()
	while(k <= ncol(mat)){
		pos <- mean(mat[5, k] > ns[k, ])
		pos.seq <- c(pos.seq, pos)
		if(pos < 1-sig.lev){
			l <- mat[6, k]+1; rm.ind <- mat[1, k]
			while(l <= max(mat[6, ])){
				ind <- which(mat[6, ]==l & is.element(mat[1, ], c(rm.ind*2-1, rm.ind*2)))
				if(length(ind) > 0){
					rm.ind <- mat[1, ind]
					mat <- mat[, -ind, drop=FALSE]
					l <- l+1
				} else{
					break
				}
			}
		}
		k <- k+1
	}
	tree <- list(matrix(0, 6, 1))
	mat <- rbind(mat, pos.seq); mat <- mat[, pos.seq >= 1-sig.lev, drop=FALSE]
	if(dim(mat)[2]>0){
		for(l in 1:length(unique(mat[6, ]))){
			j <- unique(mat[6, ])[l]
			for(ncc in 1:sum(mat[6, ]==j)){
				k <- sum(mat[6, ] < j) + ncc
				if(length(tree)<j) tree <- c(tree, list(matrix(0, 6, 0)))
				tree[[l]] <- matrix(c(tree[[l]], matrix(0, 6, 1)), 6, ncc)
				tree[[l]][1, ncc] <- mat[1, k]
				tree[[l]][2, ncc] <- mat[2, k]
				tree[[l]][3, ncc] <- mat[3, k]
				tree[[l]][4, ncc] <- mat[4, k]
				tree[[l]][5, ncc] <- mat[8, k]
				tree[[l]][6, ncc] <- mat[5, k]
			}
		}
	}
	est.cps <- sort(mat[3, ])

	ls <- list(tree=tree, est.cps=est.cps+2^(-min(scales))-1, p=p)
	return(ls)
}

idio.null.stat <- function(mat, p, B, hat.vep, diag, scales, mby, tby, do.parallel){
	len <- ncol(hat.vep)
  if(do.parallel){
    null.stat <- foreach::foreach(l=iterators::iter(1:B), .combine=cbind, .packages=c('Rcpp', 'RcppArmadillo', 'factorcpt')) %dopar% {
  		ind <- c()
	  	while(length(ind)<len){
		  		L <- rgeom(1, p); I <- sample(len, 1)
			  	ind <- c(ind, rep(1:len, 1+ceiling(L/len))[I:(I+L-1)])
      }
		  ind <- ind[1:len]
		  boot.vep <- hat.vep[, ind]
		  by <- NULL
		  for(sc in scales){
		    cc <- func_coef(boot.vep, sc)
		    if(sc>min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
		    if(diag){
		      by <- rbind(by, t(func_input_on(cc)))
		    } else{
		      sgn <- sign(cc%*%t(cc))
		      by <- rbind(by, t(func_input(cc, sgn)))
		    }
		  }
		  tmp <- c()
		  for(i in 1:ncol(mat)){
		    s <- mat[2, i]; e <- mat[4, i]
		    bfd <- func_dc_by(by[, s:e], mby, tby)
		    tmp <- c(tmp, max(bfd$res))
		    rm(bfd)
		  }
		  rm(by)
		  tmp
	  }
  } else{
    null.stat <- foreach::foreach(l=iterators::iter(1:B), .combine=cbind, .packages=c('Rcpp', 'RcppArmadillo', 'factorcpt')) %do% {
      ind <- c()
      while(length(ind)<len){
        L <- rgeom(1, p); I <- sample(len, 1)
        ind <- c(ind, rep(1:len, 1+ceiling(L/len))[I:(I+L-1)])
      }
      ind <- ind[1:len]
      boot.vep <- hat.vep[, ind]
      by <- NULL
      for(sc in scales){
        cc <- func_coef(boot.vep, sc)
        if(sc>min(scales)) cc <- cc[, -(1:(2^(-min(scales))-2^(-sc))), drop=FALSE]
        if(diag){
          by <- rbind(by, t(func_input_on(cc)))
        } else{
          sgn <- sign(cc%*%t(cc))
          by <- rbind(by, t(func_input(cc, sgn)))
        }
      }
      tmp <- c()
      for(i in 1:ncol(mat)){
        s <- mat[2, i]; e <- mat[4, i]
        bfd <- func_dc_by(by[, s:e], mby, tby)
        tmp <- c(tmp, max(bfd$res))
        rm(bfd)
      }
      rm(by)
      tmp
    }
	}
	null.stat
}

##########################################

tri.kern <- function(h){
    filter <- rep(0, h+1)
    i <- 0
    while (i <= h) {
        u <- i/h
        if (u < 1/2)
            filter[i+1] <- 1
        if (u >= 1/2 & u < 1)
            filter[i+1] <- 2 * (1 - u)
        if (u > 1)
            break
        i <- i + 1
    }
    filter
}

get.gg <- function(z, M=NULL, C=2, max.K=5){
	len <- length(z)
	max.K <- max(max.K, sqrt(log(len)))
	acv <- acf(z, type="covariance", lag.max=len-1, plot=FALSE)$acf[,,1]
	if(is.null(M)){
		l <- 1; ind <- 0
		while(l < sqrt(len)){
			if(abs(acv[l+1])/acv[1] < C*sqrt(log(len)/len)){
				ind <- ind+1
			} else{
				if(ind>0) ind <- 0
			}
			if(ind==max.K) break
			l <- l+1
		}
		lam <- max(1/2, l-max.K); M <- 2*lam
	}
	k <- tri.kern(M)
	c(acv[1]+2*sum(k[-1]*acv[2:(M+1)]), 2*sum(k[-1]*(1:M)*acv[2:(M+1)]))
}

##########################################
post.cpts.analysis <- function(x, est.cps, cutoff.seq = seq(.5, .95, by=.05), do.plot=TRUE){

  T <- dim(x)[2]
  B <- length(est.cps)
  brks <- c(0, est.cps, T)
  C <- length(cutoff.seq)

  heat.mat <- matrix(0, ncol=C, nrow=T)
  for(b in 1:(B+1)){
    s <- brks[b]+1; e <- brks[b+1]
    len <- e-s+1
    z <- x[, s:e, drop=FALSE]
    z <- z - apply(z, 1, mean)
    z <- z[!apply(z, 1, sd)==0,]
    z <- t(scale(t(z)))
    eig <- eigen(z%*%t(z)/len)
    z.eval <- eig$values[1:min(dim(z)-1)]
    for(c in 1:C) heat.mat[s:e, c] <- min(which(cumsum(z.eval)/sum(z.eval) > cutoff.seq[c]))
  }
  if(do.plot){
    par(mar=c(4, 4.5, 2, 7))
    image(heat.mat, col=tim.colors(diff(range(heat.mat))+1), breaks=min(heat.mat):(max(heat.mat)+1),
          xlab='time', ylab='c', axes=FALSE)
    abline(v=est.cps/T, col=8, lwd=2)
    axis(side=1, at=seq(0, 1, length.out=10), labels=round(seq(1, T, length.out=10)))
    axis(side=2, at=seq(0, 1, length.out=10), labels=seq(min(cutoff.seq), max(cutoff.seq), length.out=10))
    image.plot(heat.mat, col=tim.colors(diff(range(heat.mat))+1), breaks=min(heat.mat):(max(heat.mat)+1),
           legend.only=TRUE)
  }
  return(heat.mat)
}

