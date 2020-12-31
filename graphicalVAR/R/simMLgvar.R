rewire <- function(x,p,directed){
  if (missing(directed)){
    directed <- !all(x == t(x))
  }
  
  if (directed){
    ind <- diag(1,ncol(x)) != 1
  } else {
    ind <- upper.tri(x)
  }
  
  # Current edges:
  curEdges <- which(x!=0 & ind, arr.ind=TRUE)
  
  # Select to rewire:
  toRewire <- which(runif(nrow(curEdges)) < p)
  
  for (i in seq_along(toRewire)){
    curZeros <- which(x==0 & ind, arr.ind=TRUE)
    dest <- sample(seq_len(nrow(curZeros)),1)
    
    x[curZeros[dest,1],curZeros[dest,2]] <- x[curEdges[toRewire[i],1],curEdges[toRewire[i],2]]
    x[curEdges[toRewire[i],1],curEdges[toRewire[i],2]] <- 0
    
    if (!directed){
      x[curZeros[dest,2],curZeros[dest,1]] <- x[curEdges[toRewire[i],2],curEdges[toRewire[i],1]]
      x[curEdges[toRewire[i],2],curEdges[toRewire[i],1]] <- 0
    }
  }
  
  return(x)
}

simMLgvar <- function(
  nTime,
  nVar,
  nPerson,
  propPositive = 0.5,
  kappaRange = c(0.25,0.5), #c(0.5,1),
  betaRange = c(0.25,0.5),
  betweenRange = c(0.25,0.5),
  rewireWithin = 0,
  betweenVar = 1,
  withinVar = .25,
  temporalOffset = 2
){
  repeat{
    
    # Obtain true signed fixed structures:
    trueKappa <- as.matrix(igraph::get.adjacency(igraph::watts.strogatz.game(1,nVar,1,0)))
    trueKappa[upper.tri(trueKappa)] <- trueKappa[upper.tri(trueKappa)] *  sample(c(-1,1),sum(upper.tri(trueKappa)),TRUE,prob=c(propPositive,1-propPositive))
    # Symmetrize:
    trueKappa[lower.tri(trueKappa)] <- t(trueKappa)[lower.tri(trueKappa)]  
    
    # Temporal effects:
    trueBeta <- diag(1,nVar)
    for (i in 1:nVar){
      trueBeta[(i+(temporalOffset-1))%%nVar+1,i ] <- sample(c(-1,1),1,propPositive)
    }
    
    # Between subjects:
    trueBetween <- as.matrix(igraph::get.adjacency(igraph::watts.strogatz.game(1,nVar,1,1)))
    trueBetween[upper.tri(trueBetween)] <- trueBetween[upper.tri(trueBetween)] *  sample(c(-1,1),sum(upper.tri(trueBetween)),TRUE,prob=c(propPositive,1-propPositive))
    # Parameterize:
    trueBetween[upper.tri(trueBetween)] <- runif(sum(upper.tri(trueBetween)),betweenRange[1],betweenRange[2]) * trueBetween[upper.tri(trueBetween)]
    # Symmetrize:
    trueBetween[lower.tri(trueBetween)] <- t(trueBetween)[lower.tri(trueBetween)]  
    
    # Add diagonal:
    diag(trueBetween) <- 1
    
    # Check them all:
    evK <- round(eigen(trueKappa)$values,10)
    evB <- round(eigen(trueBeta)$values,10)
    evBet <- round(eigen(trueBetween)$values,10)
    
    if (all(evBet > 0)){
      break
    }
    
    
  }
  # Simulate means:
  Sigma <- cov2cor(solve(trueBetween))
  D <- diag(sqrt(betweenVar),nVar)
  Means <- mvtnorm::rmvnorm(nPerson,sigma = D%*%Sigma%*%D)
  
  # Simulate for every subject:
  SubjectData <- lapply(1:nPerson,function(i){
    try <- 1
    maxtry <- 10
    repeat{
      kappa <- trueKappa
      kappa[upper.tri(kappa)] <- runif(sum(upper.tri(kappa)),kappaRange[1],kappaRange[2]) * kappa[upper.tri(kappa)]
      # Symmetrize:
      kappa[lower.tri(kappa)] <- t(kappa)[lower.tri(kappa)]  
      diag(kappa) <- 1
      
      # Temporal:
      beta <- trueBeta * runif(nVar^2, betaRange[1], betaRange[2])
      # diag(beta) <- Vmin
      
      # Rewire:
      kappa <- rewire(kappa,rewireWithin)
      beta <- rewire(beta,rewireWithin)
      
      evK <- eigen(kappa)$values
      evB <- eigen(beta)$values
      
      while(any(Re(evB)^2 + Im(evB)^2 > 1)){
        warning("Shrinking parameters")
        beta <- 0.95*beta
        evB <- eigen(beta)$values
      }
      
      
      if (all(evK > 0) & all(Re(evB)^2 + Im(evB)^2 < 1)){
        break
      }
      
      try <- try + 1
      if (try > maxtry){
        stop("Maximum number of tries reached.")
      }
    }
    
    D <- diag(sqrt(withinVar), nVar)
    Delta <- diag(1/sqrt(diag(solve(kappa))))
    kappa <- solve(D)%*%solve(Delta)%*%kappa%*%solve(Delta)%*%solve(D)
    # Simulate data:
    Data <- as.data.frame(graphicalVARsim(nTime,beta,kappa,mean = Means[i,]))
    Data$ID <- i
    
    # Return:
    return(list(
      kappa=kappa,
      beta=beta,
      PCC = computePCC(kappa),
      PDC = computePDC(beta,kappa),
      data=Data
    ))
  })
  
  # Compute fixed effects:
  fixedKappa <- Reduce("+",lapply(SubjectData,"[[","kappa")) / nPerson
  fixedBeta <- Reduce("+",lapply(SubjectData,"[[","beta")) / nPerson
  
  # Aggregate data:
  allData <- do.call(rbind,lapply(SubjectData,"[[","data"))
  
  Results <- list(
    data = allData,
    fixedKappa = fixedKappa,
    fixedPCC = computePCC(fixedKappa),
    fixedBeta = fixedBeta,
    fixedPDC = computePDC(fixedBeta,fixedKappa),
    between = trueBetween,
    means=Means,
    personData = SubjectData,
    idvar = "ID",
    vars = names(allData)[names(allData)!="ID"]
  )
  
  class(Results) <- "simMLgvar"
  return(Results)
}
