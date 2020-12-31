HWTriExact <- function(x,y=NULL,eps=1e-10,nperm=17000,verbose=TRUE) {
    test.type <- NULL
    x <- unlist(x)
    if(!is.null(y)) y <- unlist(y)
  if(length(x)==6 & is.null(y))   test.type <- 1 # ordinary autosomal
  if(length(x)==6 & length(y)==6) test.type <- 2 # autosomal accounting for sex
  if(length(x)==6 & length(y)==3) test.type <- 3 # x-chromosomal accounting for sex 
    if(is.null(test.type)) stop("incorrect number of genotype counts in x or y")
    pseudodist <- NULL
    
  if(test.type==1) {
    if(is.vector(x)) G <- toTriangular(x) else G <- x
    pofthesample <- dlevene(G) 
    n.a <- sort(colSums(G) + rowSums(G))
    na <- n.a[1]; nb <- n.a[2]; nc <- n.a[3]
    O <- outcomes.3(n.a)
    ntab <- nrow(O)
    pr <- numeric(nrow(O))
    for(i in 1:nrow(O)) {
      pr[i] <- dlevene(toTriangular(O[i,]))  
    }

    ii <- nearlyEqual(pr, rep(pofthesample, nrow(O)), eps)
    pval <- sum(pr[ii]) # sum of all tied samples
    iii <- ((!ii) & (pr < pofthesample))
    pval <- pval + sum(pr[iii])

    if(verbose) {
      cat("Tri-allelic Exact test for HWE (autosomal).\n")
      cat("Allele counts: A =",na,"B =",nb,"C =",nc,"\n"); 
      cat("sum probabilities all outcomes",sum(pr),"\n")
      cat("probability of the sample",pofthesample,"\n")
      cat("p-value = ",pval,"\n")  
    }
  } # end test.type==1
  if(test.type==2) { # autosomal accounting for gender
    out <- HWPerm.mult(toTriangular(x),toTriangular(y),nperm=nperm)
    pval <- out$pval
    pofthesample <- out$pofthesample
    pseudodist <- out$pseudodist
  }
  if(test.type==3) { # x-chromosomal test
    m <- y
    f <- x
    na <- m[1] + 2*f[1] + f[2] + f[3]
    nb <- m[2] + 2*f[4] + f[2] + f[5]
    nc <- m[3] + 2*f[6] + f[3] + f[5]
    f <- toTriangular(f)
    pofthesample <- density.ma.gender(m,f)
    z <- sort(c(na,nb,nc))
    na <- z[1]
    nb <- z[2]
    nc <- z[3]
    nm <- sum(m)
    nf <- sum(f)
    X <- gen.outcomesXtri(na,nb,nc,nm,nf)
    pr <- numeric(nrow(X))
    for (i in 1:nrow(X)) {
      ma <- X[i,1:3]
      fe <- toTriangular(X[i,4:9])
      pr[i] <- density.ma.gender(ma,fe)  
    }

    ii <- nearlyEqual(pr, rep(pofthesample, nrow(X)), eps)
    pval <- sum(pr[ii]) # sum of all tied samples
    iii <- ((!ii) & (pr < pofthesample))
    pval <- pval + sum(pr[iii])

    if(verbose) {
      cat("Tri-allelic Exact test for HWE and EAF (X-chromosomal)\n")
      cat("Allele counts: na = ",na,"nb = ",nb,"nc =",nc,"\n")
      cat("Sample contains: ",nm,"males and",nf,"females\n")
      cat("sum probabilities all outcomes",sum(pr),"\n")
      cat("probability of the sample",pofthesample,"\n")
      cat("p-value = ",pval,"\n")  
    }
  }
  list(pval=pval,pseudodist=pseudodist,pofthesample=pofthesample)
}


