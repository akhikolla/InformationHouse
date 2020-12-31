HWPerm.mult <- function(x, y = NULL, nperm = 17000, eps = 1e-10, verbose = TRUE, ...) {
  testtype <- NULL
  
  if(is.null(y)) { 
    if(is.vector(x) & length(x)==3) {
      testtype <- 1   # ordinary autosomal bi-allelic  
    }
    if(is.matrix(x)) {
      testtype <- 5 # Autosomal; k alleles; no gender.
    }
    
  }
  
  if((is.vector(x) & length(x)==2) & (is.vector(y) & length(y)==3)) {
    testtype <- 2  # x-chromosomal bi-allelic
  }
  
  if( (is.vector(x) & is.vector(y) & length(x)==length(y)) | 
      (is.matrix(x) & is.matrix(y) & all(dim(x)==dim(y)))) {
    testtype <- 3  # autosomal k-allelic stratified for gender
  }
 
  if(!is.null(y)) {
     if(is.vector(x) & is.matrix(y)) {
          if(length(x)==nrow(y)) {
              testtype <- 4  # X-chromosomal k-allelic
          }
      }
     if(is.vector(x) & length(x) > 2 & is.vector(y) & length(y) > 3) {
              testtype <- 4  # X-chromosomal k-allelic
    }   
  } 
 
  if(is.null(testtype)) {
    stop("data vectors x and/or y not correctly specified.")
  }
  
  if(testtype==1) {   # ordinary autosomal bi-allelic  
    nAA <- x["AA"]
    nAB <- x["AB"]
    nBB <- x["BB"]
    
    n <- sum(x)
    nA <- 2 * nAA + nAB
    nB <- 2 * n - nA
    
    pofthesample <- dlevene.bi(x)
    
    pseudodist <- numeric(nperm)
    i1 <- seq(1, 2 * n, 2)
    i2 <- seq(2, 2 * n, 2)
    for (i in 1:nperm) {
      xx <- sample(c(rep("A", nA), rep("B", nB)))
      A1 <- xx[i1]
      A2 <- xx[i2]
      Geno <- paste(A1, A2, sep = "")
      Geno[Geno == "BA"] <- "AB"
      nAA <- sum(Geno == "AA")
      nAB <- sum(Geno == "AB")
      nBB <- sum(Geno == "BB")
      y <- c(AA = nAA, AB = nAB, BB = nBB)
      stat.psuedo <- dlevene.bi(y)
      pseudodist[i] <- stat.pseudo
    }

    ii <- nearlyEqual(pseudodist, rep(pofthesample, nperm), eps)
    nsmaller <- sum(ii)
    iii <- !ii & pseudodist < pofthesample
    nsmaller <- nsmaller + sum(iii)
    pval <- nsmaller/nperm
    
    if (verbose) {
        cat("Permutation test for Hardy-Weinberg equilibrium (autosomal).\n")
        cat("2 alleles detected.\n")
        cat("Observed statistic:", pofthesample, " ", nperm, 
          "permutations. p-value:", pval, "\n")
    }
  }
  
  if(testtype == 2) { # x-chromosomal bi-allelic
    if(length(x)==2 & length(y)==3) {
      xlab <- names(x)
      ylab <- names(y)
    }
    nm <- sum(x)
    nf <- sum(y)
    nfAA <- y[ylab == "AA"]
    nfAB <- y[ylab == "AB"]
    nfBB <- y[ylab == "BB"]
    nmA <- x[xlab == "A"]
    nmB <- x[xlab == "B"]
    
    nA <- nmA + 2 * nfAA + nfAB
    nB <- nmB + 2 * nfBB + nfAB
    nt <- nA + nB
    pofthesample <- dgraffelmanweir.bi(x,y)
    pseudodist <- numeric(nperm)
    for (i in 1:nperm) {
      xx <- sample(c(rep("A", nA), rep("B", nB)))
      males <- xx[1:nm]
      nmAsim <- sum(males == "A")
      nmBsim <- sum(males == "B")
      females <- xx[(nm + 1):nt]
      i1 <- seq(1, 2 * nf, 2)
      i2 <- seq(2, 2 * nf, 2)
      A1 <- females[i1]
      A2 <- females[i2]
      Geno <- paste(A1, A2, sep = "")
      Geno[Geno == "BA"] <- "AB"
      nfAAsim <- sum(Geno == "AA")
      nfABsim <- sum(Geno == "AB")
      nfBBsim <- sum(Geno == "BB")
      x.s <- c(A = nmAsim, B = nmBsim)
      y.s <- c(AA = nfAAsim, AB = nfABsim, BB = nfBBsim)
      stat.pseudo <- dgraffelmanweir.bi(x.s,y.s)
      pseudodist[i] <- stat.pseudo
    }

    ii <- nearlyEqual(pseudodist, rep(pofthesample, nperm), eps)
    nsmaller <- sum(ii)
    iii <- !ii & pseudodist < pofthesample
    nsmaller <- nsmaller + sum(iii)
    pval <- nsmaller/nperm

    if (verbose) {
        cat("Permutation test for Hardy-Weinberg equilibrium and equality of allele frequencies (X-chromosomal).\n")
        cat("2 alleles detected.\n")
      cat("Observed statistic:", pofthesample, " ", nperm, 
          "permutations. p-value:", pval, "\n")
    }
  }
  
  if(testtype==3) { # autosomal k-allelic stratified for gender
  
    nm <- sum(x) # number of males
    nf <- sum(y) # number of females
    
    if(is.vector(x) & is.vector(y)) { # if vectors convert to triangular
      ngm <- length(x)
      ngf <- length(y)
      if(ngm!=ngf) stop("number of genotypes not equal for males and females")
      k <- round(-0.5+0.5*sqrt(1+8*ngf),digits=0) # number of alleles
      m <- toTriangular(x)
      f <- toTriangular(y)
    }
    
    if(is.matrix(x) & is.matrix(y)) {
        k <- nrow(x)
        m <- x
        f <- y
    } 
    
    pofthesample <- dens.auto(m,f)

    alsnames.m <- colnames(m)
    alsnames.f <- colnames(f)
      
    if(any(alsnames.m!=alsnames.f)) stop("different male and female alleles")
      
    mm.counts <- rowSums(m)+colSums(m) # male allele sample counts
      
    als.m <- rep(alsnames.m,times=mm.counts) # vector with all male alleles
      
    ff.counts <- rowSums(f)+colSums(f) # female allele sample counts
      
    als.f <- rep(alsnames.f,times=ff.counts) # vector with all female alleles
      
    alvec <- c(als.m,als.f) # vector with all alleles
      
    pseudodist <- numeric(nperm)
      
    for(l in 1:nperm) {
      alvec.scrambled <- sample(alvec)
      males   <- alvec.scrambled[1:(2*nm)]
      females <- alvec.scrambled[(2*nm+1):length(alvec)]
        #    print(males)
      ind1.m <- seq(1,2*nm,2)
      ind2.m <- seq(2,2*nm,2)
      genotypes.m <- paste(males[ind1.m],males[ind2.m],sep="") # paste alleles off into genotypes  
      genotypes.m.s <- sapply(genotypes.m,strsort) # sort the alleles in the male genotypes
      ind1.f <- seq(1,2*nf,2)
      ind2.f <- seq(2,2*nf,2)
      genotypes.f <- paste(females[ind1.f],females[ind2.f],sep="")
      genotypes.f.s <- sapply(genotypes.f,strsort) # sort the alleles in the female genotype
      tab.m <- table(genotypes.m.s) # table of observed male genotype counts 
      tab.f <- table(genotypes.f.s) # table of observed female genotype counts 
        
      gmatrix <- matrix(NA,nrow=k,ncol=k) 
        
      for(i in 1:k) {  # matrix of all possible genotypes
        for(j in 1:k) {
          gmatrix[j,i] <- paste(alsnames.m[i],alsnames.m[j],sep="")
        }
      }
        
      gmatrix.f <- matrix(0,nrow=k,ncol=k) # matrix of observed female genotype counts
      rownames(gmatrix.f) <- alsnames.f
      colnames(gmatrix.f) <- alsnames.f
        
      for(i in 1:k) {
        for(j in 1:k) {
          ind <- gmatrix[i,j]==names(tab.f)
          if(any(ind)) gmatrix.f[i,j] <- tab.f[ind]
        }
      }
        
      #    print(gmatrix.n)
        
      gmatrix.m <- matrix(0,nrow=k,ncol=k) # matrix of observed female genotype counts
      rownames(gmatrix.m) <- alsnames.m
      colnames(gmatrix.m) <- alsnames.m
        
      for(i in 1:k) {
        for(j in 1:k) {
          ind <- gmatrix[i,j]==names(tab.m)
          if(any(ind)) gmatrix.m[i,j] <- tab.m[ind]
        }
      }
        
      #    dens.auto(gmatrix.m,gmatrix.f)
      pseudodist[l] <- dens.auto(gmatrix.m,gmatrix.f)
    } # end for
      
    #  print(pseudodist)
      
    ii <- nearlyEqual(pseudodist, rep(pofthesample, nperm), eps)
    nsmaller <- sum(ii)
    iii <- !ii & pseudodist < pofthesample
    nsmaller <- nsmaller + sum(iii)
    pval <- nsmaller/nperm
    
    if (verbose) {
        cat("Permutation test for Hardy-Weinberg equilibrium and equality of allele frequencies (autosomal).\n")
        cat(k,"alleles detected.\n")
      cat("Observed statistic:", pofthesample, " ", nperm, 
          "permutations. p-value:", pval, "\n")
    }
  } # end testtype 3
  if(testtype==4) {
    nm <- sum(x) # number of males
    nf <- sum(y) # number of females

    if(is.vector(x) & is.vector(y)) { # if vectors convert females to triangular
      nal.m <- length(x)
      nal.f <- nrow(toTriangular(y))
      if(nal.m!=nal.f) stop("number of alleles not equal for males and females")
      m <- x
      f <- toTriangular(y)
  } else {
      m <- x
      f <- y
  }
        
    k <- length(m) # number of alleles
          
    pofthesample <- density.ma.gender(m,f)

    als <- names(m)
    als.m <- rep(als,times=m) # vector with all male alleles
    ff.counts <- rowSums(f)+colSums(f) # female allele sample counts
    als.f <- rep(als,times=ff.counts) # vector with all female alleles
    alvec <- c(als.m,als.f) # vector with all alleles
        
    pseudodist <- numeric(nperm)
        
    for(l in 1:nperm) {
      alvec.scrambled <- sample(alvec)   
      males   <- alvec.scrambled[1:nm]
      females <- alvec.scrambled[(nm+1):length(alvec)]
        #    print(males)
          #    print(females)
          
      ind1 <- seq(1,2*nf,2)
      ind2 <- seq(2,2*nf,2)
      genotypes.f <- paste(females[ind1],females[ind2],sep="")
      genotypes.s <- sapply(genotypes.f,strsort) # sort the alleles in the female genotypes
      m.counts <- numeric(k)
      for(i in 1:length(als)) {
        m.counts[i] <- sum(males==als[i]) # determine male genotype counts   
      }
          
      tab <- table(genotypes.s) # table of observed female genotype counts 
          #    print(tab)  
      gmatrix <- matrix(NA,nrow=k,ncol=k) 
          
      for(i in 1:k) {  # matrix of all possible genotypes
        for(j in 1:k) {
          gmatrix[j,i] <- paste(als[i],als[j],sep="")
        }
      }
          
      gmatrix.n <- matrix(0,nrow=k,ncol=k) # matrix of observed female genotype counts
      rownames(gmatrix.n) <- als
      colnames(gmatrix.n) <- als
          
        
      for(i in 1:k) {
        for(j in 1:k) {
          ind <- gmatrix[i,j]==names(tab)
          if(any(ind)) gmatrix.n[i,j] <- tab[ind]
        }
      }
          
        #    print(gmatrix.n)
      pseudodist[l] <- density.ma.gender(m.counts,gmatrix.n)
          
    } # end for
        
      #  print(pseudodist)
        
    ii <- nearlyEqual(pseudodist, rep(pofthesample, nperm), eps)
    nsmaller <- sum(ii)
    iii <- !ii & pseudodist < pofthesample
    nsmaller <- nsmaller + sum(iii)
    pval <- nsmaller/nperm

    if (verbose) {
        cat("Permutation test for Hardy-Weinberg equilibrium and equality of allele frequencies (X-chromosomal).\n")
        cat(k,"alleles detected.\n")
      cat("Observed statistic:", pofthesample, " ", nperm, 
            "permutations. p-value:", pval, "\n")
    }
  }
  if(testtype ==5) { # autosomal; k alleles; no gender.
    #
    # MC test for the autosomal case k alleles (sexes not distinguished)
    #
    
#    mcarlotest.auto <- function(x,nperm=100,verbose=TRUE,eps=1e-10,logscale=FALSE) {
    n <- sum(x) # total sample size

    k <- nrow(x) # number of alleles
    if(k != ncol(x)) stop("error x is not square.")
   
    pofthesample <- dlevene(x)

    als <- colnames(x)
    ac.counts <- rowSums(x)+colSums(x) # sample allele counts
    alvec <- rep(als,times=ac.counts) # vector with all alleles
      
    pseudodist <- numeric(nperm)
      
    for(l in 1:nperm) {
      scrambled <- sample(alvec)   
      ind1 <- seq(1,2*n,2)
      ind2 <- seq(2,2*n,2)
      genotypes <- paste(scrambled[ind1],scrambled[ind2],sep="")
        # pair female alleles off into genotypes
        #    print(genotypes.f)
      genotypes.s <- sapply(genotypes,strsort) # sort the alleles in the female genotypes
        #      print(genotypes.s)
      tab <- table(genotypes.s) # table of observed female genotype counts 
        #       print(tab)  
      gmatrix <- matrix(NA,nrow=k,ncol=k) 
      for(i in 1:k) {  # matrix of all possible genotypes
        for(j in 1:k) {
          gmatrix[j,i] <- paste(als[i],als[j],sep="")
        }
      }
        
      gmatrix.n <- matrix(0,nrow=k,ncol=k) # matrix of observed genotype counts
      rownames(gmatrix.n) <- als
      colnames(gmatrix.n) <- als
        
      for(i in 1:k) {
        for(j in 1:k) {
          ind <- gmatrix[i,j]==names(tab)
          if(any(ind)) gmatrix.n[i,j] <- tab[ind]
        }
      }
      
      pseudodist[l] <- dlevene(gmatrix.n)
        
    } # end for
      
    ii <- nearlyEqual(pseudodist, rep(pofthesample, nperm), eps)
    nsmaller <- sum(ii)
    iii <- !ii & pseudodist < pofthesample
    nsmaller <- nsmaller + sum(iii)
    pval <- nsmaller/nperm
    
    if (verbose) {
        cat("Permutation test for Hardy-Weinberg equilibrium (autosomal).\n")
        cat(k,"alleles detected.\n")
      cat("Observed statistic:", pofthesample, " ", nperm, 
          "permutations. p-value:", pval, "\n")
    }
    
  }
  return(list(pofthesample = pofthesample, pseudodist=pseudodist, pval = pval))
}
