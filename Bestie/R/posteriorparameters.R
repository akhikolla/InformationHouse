
DAGparameters <- function(incidence, dataParams) {
  # add parameters to an incidence matrix

  if (dataParams$type != "bde") {
    stop("The implementation is currently only for the BDe score.")
  }
  
  if (dataParams$DBN) {

    n <- dataParams$n # number of nodes including background
    nsmall <- dataParams$nsmall # number of nodes without background
    bgn <- dataParams$bgn # number of background nodes
    slices <- dataParams$slices # number of time slices

### First slice parameters    
# remember background nodes were placed at the end of first slice, which we need to undo!
    
    dataParams_first <- dataParams$firstslice
    if(bgn > 0){ # place backgound nodes at the start again
      reorder <- c(1:bgn + nsmall, 1:nsmall)
      dataParams_first$data <- dataParams_first$data[, reorder]
      dataParams_first$d1 <- dataParams_first$d1[, reorder]
      dataParams_first$d0 <- dataParams_first$d0[, reorder]
    }
    params_first <- DAGparameters(incidence[1:n, 1:n], dataParams_first)

## Other slices parameters    
# remember the later times were placed first, which we need to undo!     

    dataParams_other <- dataParams$otherslices

    reorder <- c(1:n + nsmall, 1:nsmall) # put them back in order
    dataParams_other$data <- dataParams_other$data[, reorder]
    dataParams_other$d1 <- dataParams_other$d1[, reorder]
    dataParams_other$d0 <- dataParams_other$d0[, reorder]
        
    params <- DAGparameters(incidence, dataParams_other)

### Combine parameters from the slices
    
    allalphas <- params$alphas
    allalphas[1:n] <- params_first$alphas # update the first slice
    allbetas <- params$betas
    allbetas[1:n] <- params_first$betas # update the first slice
    allpmeans <- params$pmeans
    allpmeans[1:n] <- params_first$pmeans # update the first slice

### Need to unroll the DBN
    
    if (slices > 2) {
      nbig <- n + nsmall*(slices - 1)
      
      incidence_unroll <- matrix(0, nbig, nbig)
      incidence_unroll[1:nrow(incidence), 1:ncol(incidence)] <- incidence

      allalphas_unroll <- vector("list", nbig)
      allbetas_unroll <- vector("list", nbig)
      allpmeans_unroll <- vector("list", nbig)
      
      allalphas_unroll[1:ncol(incidence)] <- allalphas
      allbetas_unroll[1:ncol(incidence)] <- allbetas
      allpmeans_unroll[1:ncol(incidence)] <- allpmeans
      
      for (ii in 1:(slices - 2)) {
        block_rows <- n - nsmall + 1:(2*nsmall)
        block_cols <- n + 1:nsmall

        incidence_unroll[block_rows + nsmall*ii, block_cols + nsmall*ii] <- incidence[block_rows, block_cols]

        allalphas_unroll[block_cols + nsmall*ii] <- allalphas[block_cols]
        allbetas_unroll[block_cols + nsmall*ii] <- allbetas[block_cols]
        allpmeans_unroll[block_cols + nsmall*ii] <- allpmeans[block_cols]
        
        if (bgn > 0) { # if there are background nodes, repeat across slices
          block_rows <- 1:bgn
          incidence_unroll[block_rows, block_cols + nsmall*ii] <- incidence[block_rows, block_cols]
        }
      }
      
      incidence <- incidence_unroll
      allalphas <- allalphas_unroll
      allbetas <- allbetas_unroll
      allpmeans <- allpmeans_unroll
    }
 
  } else {
  
  n <- nrow(incidence) # number of nodes in DAG

  allalphas <- vector("list", n)
  allbetas <- vector("list", n)
  allpmeans <- vector("list", n)

  for (j in 1:n) {
    parentNodes <- which(incidence[, j]==1)
    tempResult <- DAGparametersCore(j, parentNodes, dataParams)
    allalphas[[j]] <- tempResult$alphas
    allbetas[[j]] <- tempResult$betas
    allpmeans[[j]] <- tempResult$pmeans
  }

  }
  
  posteriorParams <- list()
  posteriorParams$DAG <- incidence
  posteriorParams$alphas <- allalphas
  posteriorParams$betas <- allbetas
  posteriorParams$pmeans <- allpmeans

  return(posteriorParams)
}



DAGparametersCore <- function(j, parentNodes, param) {
  # this function computes the parameters and their posterior distribution
  # for a given node with parent set and for the data

  switch(param$type,
    "bge" = {
      stop("This option is not implemented")
    },
    "bde" = {
      lp <- length(parentNodes) # number of parents
      noParams <- 2^lp # number of binary states of the parents
      chi <- param$chi

      alphas <- rep(NA, noParams)
      betas <- rep(NA, noParams)

      switch(as.character(lp),
        "0"={ # no parents
          N1 <- sum(param$d1[, j])
          N0 <- sum(param$d0[, j])
          NT <- N0 + N1
          alphas <- (N1 + chi/(2*noParams))
          betas <- (N0 + chi/(2*noParams))
        },
        "1"={ # one parent
          summys <- param$data[, parentNodes]
          for (i in 1:noParams-1) {
            totest <- which(summys==i)
            N1 <- sum(param$d1[totest, j])
            N0 <- sum(param$d0[totest, j])
            NT <- N0 + N1
            alphas[i+1] <- (N1 + chi/(2*noParams))
            betas[i+1] <- (N0 + chi/(2*noParams))
           }
         },
         { # more parents
           summys <- colSums(2^(c(0:(lp-1)))*t(param$data[, parentNodes]))
           N1s <- collectC(summys, param$d1[, j], noParams)
           N0s <- collectC(summys, param$d0[, j], noParams)
           NTs <- N1s + N0s
           alphas <- (N1s + chi/(2*noParams))
           betas <- (N0s + chi/(2*noParams))
         }
      )

      coreParams <- list()
      coreParams$alphas <- alphas
      coreParams$betas <- betas
      coreParams$pmeans <- alphas/(alphas + betas)

      return(coreParams)
    }
  )
}



SampleParameters <- function(DAGparams) {
  # this function resamples the probability parameters from the posterior beta distributions

  sampledps <- DAGparams$pmeans

  n <- length(sampledps)

  for(jj in 1:n){
    as <- DAGparams$alphas[[jj]]
    bs <- DAGparams$betas[[jj]]
    ps <- rep(0,length(as))
    for(ii in 1:length(as)){
      ps[ii] <- stats::rbeta(1, as[ii], bs[ii])
    }
    sampledps[[jj]] <- ps
  }

  return(sampledps)
}



BinaryScoreAgainstDAG <- function(DAGparams, dataToScore) {
  # score a set of binary vectors against a DAG with parameters

  n <- nrow(DAGparams$DAG) # number of nodes

  logscoresagainstDAG<-matrix(NA,nrow(dataToScore),n)
  for (j in 1:n)  {
    parentNodes <- which(DAGparams$DAG[, j]==1)
    logscoresagainstDAG[, j] <- BinaryScoreAgainstDAGcore(j, parentNodes, DAGparams, dataToScore)
  }

  return(logscoresagainstDAG)
}



BinaryScoreAgainstDAGcore <- function(j, parentNodes, DAGparams, dataToScore) {
  # score of a single node of binary vectors against a DAG

  sampleNodeScores <- rep(NA, nrow(dataToScore)) # store the log scores

  lp <- length(parentNodes) # number of parents
  noParams <- 2^lp # number of binary states of the parents

  switch(as.character(lp),
    "0"={ # no parents
      theta <- DAGparams$pmeans[[j]] # the probability of each state
      sampleNodeScores[which(dataToScore[, j]==1)] <- log(theta) # log scores of 1s
      sampleNodeScores[which(dataToScore[, j]==0)] <- log(1-theta) # log scores of 0s
    },
    "1"={ # one parent
      summysfull<-dataToScore[,parentNodes]

      for (i in 1:noParams-1) {
        theta <- DAGparams$pmeans[[j]][i+1] # the probability of each state
        toScore <- which(summysfull==i)
        sampleNodeScores[toScore[which(dataToScore[toScore,j]==1)]] <- log(theta) # log scores of 1s
        sampleNodeScores[toScore[which(dataToScore[toScore,j]==0)]] <- log(1-theta) # log scores of 0s
      }
    },
    { # more parents
      summysfull <- colSums(2^(c(0:(lp-1)))*t(dataToScore[, parentNodes]))

      # find the entries where the child is 1
      toScore<-which(dataToScore[, j]==1)

      #Ns <- collectC(summysfull[toScore], rep(1, length(toScore)), noParams) # works like the table command

      Ns <- tabulate(summysfull[toScore]+1,noParams) # can use tabulate instead of collectC, but we need to add one
      tempScoreVec <- rep(log(DAGparams$pmeans[[j]]), Ns) # make relevant number of copies of each log score
      sampleNodeScores[toScore] <- tempScoreVec[rank(summysfull[toScore], ties.method="first")] # use the rank function to map scores to entries

      # find the entries where the child is 0
      toScore<-which(dataToScore[,j]==0)

      Ns<-tabulate(summysfull[toScore]+1,noParams) # again we need to add one
      tempScoreVec<-rep(log(1-DAGparams$pmeans[[j]]),Ns) # make relevant number of copies of each log score
      sampleNodeScores[toScore]<-tempScoreVec[rank(summysfull[toScore],ties.method="first")] # use the rank function to map scores to entries
    }
  )

  return(sampleNodeScores)
}


