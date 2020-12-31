#' Monte Carlo estimation of intervention effects for a DAG or chain of sampled DAGs
#'
#' \code{DAGinterventionMC} takes a DAG or a sampled chain of DAGs (for example from
#' the \code{\link[BiDAG]{partitionMCMC}} function of the BiDAG package) and computes a Monte Carlo
#' estimate of the intervention effect of each node on all others by simulating data
#' from the DAG. By default each node is intervened upon and the downstream effects
#' estimated by further sampling. A faster but less robust and accurate version is
#' also offered which reweights a single simulated dataset.
#'
#' @param incidences a single adjacency matrix of a list of adjacency matrices
#' of sampled DAGs, with entry [i,j] equal to 1 when a directed edge exists from
#' node i to node j
#' @param dataParams the data and parameters used to learn the DAGs derived from the
#' \code{\link[BiDAG]{scoreparameters}} function of the BiDAG package
#' @param sampleSize the number of Monte Carlo samples to draw
#' @param sample logical indicating whether to sample the parameters of each node
#' from the posterior (TRUE, default) or to take the expectation (FALSE)
#' @param fixNode logical indicating whether to intervene on each node (TRUE, default)
#' and resample downstream nodes or to sample once and reweight the sample (FALSE)
#' @param reducedVarianceSampling logical indicating whether to perform Bernoulli
#' samping for each node (FALSE) or to sample from a distribution with the same mean
#' and lower variance (TRUE, default)
#'
#' @return a single matrix or a list of matrices containing the full set of
#' intervention effects for each input DAG. Entry [i,j] is the downstream
#' effect on node j of intervening on node i
#' (the difference observed at node j when setting node i to 1 and 0)
#'
#' @examples
#'
#' scoreParam <- BiDAG::scoreparameters(8, "bde", BiDAG::Asia)
#' causalmatMC <- DAGinterventionMC(BiDAG::Asiamat, scoreParam, 1e4)
#'
#' @export
#'
#' @seealso \code{\link[BiDAG]{scoreparameters}}

DAGinterventionMC <- function(incidences, dataParams, sampleSize, sample = TRUE, fixNode = TRUE, reducedVarianceSampling = TRUE){
  # this wrapper takes in a chain of DAG, computes their parameters and returns all MC intervention effects

  if(sampleSize < 1e3) {
    sampleSize <- 1e3
    warning("Setting sample size to a minimum of 1000.")
  }

  if (!is.list(incidences)) { # turn to list internally
    incidences <- list(incidences)
  }

  numDAGs <- length(incidences)
  interventionMats <- vector("list", numDAGs) # to store the intervention effects

  for(kk in 1:numDAGs){
    DAGparams <- DAGparameters(incidences[[kk]], dataParams)
    interventionMats[[kk]] <- DAGinterventionMCparams(DAGparams, sampleSize, sample, fixNode, reducedVarianceSampling)
  }

  if (numDAGs == 1) { # turn back to single matrix
    interventionMats <- interventionMats[[1]]
  }

  return(interventionMats)
}



DAGinterventionMCparams <- function(DAGparams, sampleSize, sample = FALSE, fixNode = TRUE, reducedVarianceSampling = TRUE){
  # this wrapper takes in a DAG with parameters and returns all intervention effects estimated via Monte Carlo
  # when fixNode is FALSE this will simply sample binary vectors from the DAG which can give problems with low probability events!

  DAGparamsInternal <- DAGparams
  if(sample==TRUE){ # then we take a sample of parameters from the posterior instead of taking the expectation
    DAGparamsInternal$pmeans <- SampleParameters(DAGparams)
  }

  if(fixNode==FALSE){ # we sample the set of binary vectors and reweight
    sampledBinaryVecs <- BinarySampleFromDAG(DAGparamsInternal, sampleSize, reducedVarianceSampling)
    sampledLogScores <- BinaryScoreAgainstDAG(DAGparamsInternal, sampledBinaryVecs)
    interventionMat <- InterventionEstimationFromSample(sampledLogScores, sampledBinaryVecs)
  } else { # we fix each node in turn and calculate the downstream effects
    interventionMat <- InterventionEstimationMCfix(DAGparamsInternal, sampleSize, reducedVarianceSampling)
  }

  return(interventionMat)
}



BinarySampleFromDAG <- function(DAGparams, sampleSize, reducedVarianceSampling=TRUE){
  # sample a set of binary vectors from a DAG with parameters

  n <- ncol(DAGparams$DAG)

  binarySample <- matrix(NA, sampleSize, n)
  nodeOrder <- DAGtoorder(DAGparams$DAG) # get topological order

  for (j in rev(nodeOrder)) { # start with outpoints etc
    parentNodes <- which(DAGparams$DAG[, j]==1)
    binarySample[, j] <- BinarySampleFromDAGcore(j, parentNodes, DAGparams, binarySample, reducedVarianceSampling)
  }

  return(binarySample)
}



BinarySampleFromDAGcore <- function(j, parentNodes, DAGparams, binarySample, reducedVarianceSampling=TRUE){
  # sample one variable for a set of binary vectors from a DAG with parameters

  sampleNode <- rep(NA, nrow(binarySample)) # store the sampled values

  lp <- length(parentNodes) # number of parents
  noParams <- 2^lp # number of binary states of the parents

  switch(as.character(lp),
    "0"={# no parents
      theta <- DAGparams$pmeans[[j]] # the probability of each state
      sampleNode <- SampleBinaryVec(theta, nrow(binarySample), reducedVarianceSampling)
    },
    "1"={# one parent
      summysfull <- binarySample[, parentNodes]

      for (i in 1:noParams - 1) {
        theta <- DAGparams$pmeans[[j]][i+1] # the probability of each state
        toScore <- which(summysfull==i)
        if (length(toScore) > 0) {
          sampleNode[toScore] <- SampleBinaryVec(theta, length(toScore), reducedVarianceSampling)
        }
      }
    },
    { # more parents
      summysfull <- colSums(2^(c(0:(lp-1)))*t(binarySample[, parentNodes]))

      Ns <- tabulate(summysfull+1,noParams) # can use tabulate instead of collectC, but we need to add one
      tempScoreVec <- rep(NA, length(summysfull))
      localCounter <- 0

      for (i in 1:noParams) { # we run over the data size once
        if (Ns[i]>0) { # only if there are states to consider
          theta <- DAGparams$pmeans[[j]][i] # the probability of each state
          tempScoreVec[localCounter + 1:Ns[i]] <- SampleBinaryVec(theta, Ns[i], reducedVarianceSampling)
          localCounter <- localCounter + Ns[i]
        }
      }
      sampleNode <- tempScoreVec[rank(summysfull, ties.method="first")] # use the rank function to map scores to entries
    }
  )

  return(sampleNode)
}



SampleBinaryVec <- function(theta, vecSize, reducedVarianceSampling=TRUE){
  # sample with same expectation as Bernoulli sampling, but with much less variance.
  # or straightforward Bernoulli sampling if reducedVarianceSampling=FALSE

  if (reducedVarianceSampling==FALSE || vecSize==1) {
    binarySample <- stats::rbinom(vecSize, 1, theta)
  } else {
    binaryVecTemp <- rep(0, vecSize)
    expectedOnes <- vecSize*theta
    certainOnes <- floor(expectedOnes)
    if (certainOnes > 0) {
      binaryVecTemp[1:certainOnes] <- 1
    }
    if (certainOnes < vecSize) {
      binaryVecTemp[certainOnes + 1] <- stats::rbinom(1, 1, expectedOnes-certainOnes)
    }
    binarySample <- sample(binaryVecTemp) # we need to have a random order!!!
  }

  return(binarySample)
}



InterventionEstimationFromSample <- function(logScores, dataToScore){
  # estimate intervention effects from a set of binary vectors
  # which have been sampled proportionally to their score

  n <- ncol(dataToScore)

  interventionMatrix <- matrix(0, nrow=n, ncol=n)
  for (j in 1:n) {
    interventionMatrix[j,] <- InterventionEstimationFromSampleCore(j, logScores, dataToScore)
  }

  return(interventionMatrix)
}



InterventionEstimationFromSampleCore <- function(fixedNode, logScores, dataToScore){
  # core function involves setting one node to 0 and to 1 and removing its probability component

  upRows <- which(dataToScore[, fixedNode]==1)

  if (length(upRows)>1) {
    weightVecTemp <- -(logScores[upRows, fixedNode]) # probability of choosing fixednode as one for each vector
    weightVec <- exp(weightVecTemp - max(weightVecTemp))
    weightVec <- weightVec/sum(weightVec)
    probs1 <- colSums(weightVec*dataToScore[upRows, ])
  } else {
    warning("Sample size too small!")
    probs1 <- rep(NA, ncol(logScores))
  }

  downRows <- which(dataToScore[, fixedNode]==0)

  if (length(downRows)>1) {
    weightVecTemp <- -(logScores[downRows, fixedNode]) # probability of choosing fixednode as zero for each vector
    weightVec <- exp(weightVecTemp - max(weightVecTemp))
    weightVec <- weightVec/sum(weightVec)
    probs0 <- colSums(weightVec*dataToScore[downRows, ])
  } else {
    warning("Sample size too small!")
    probs0 <- rep(NA, ncol(logScores))
  }

  return(probs1-probs0)
}



NodeDescendantsFromDAG <- function(fixedNode, DAG){
  # compute the descendants of a node from the adjacency matrix
  # this could be more efficient
  descendants <- c()
  children <- which(DAG[fixedNode, ]==1) # the children of the node
  newDescendants <- children
  while(length(newDescendants) > 0 && length(children) > 0) {
    newDescendants <- setdiff(children, descendants)
    if(length(newDescendants) > 1) {
      children <- which(colSums(DAG[newDescendants, ]) > 0)
    } else if (length(newDescendants)==1) {
      children <- which(DAG[newDescendants, ]==1)
    } else {
      children <- c()
    }
    descendants <- union(newDescendants, descendants)
  }

  return(descendants)
}



InterventionEstimationMCfix <- function(DAGparams, sampleSize, reducedVarianceSampling=TRUE){
  # estimate intervention effects by sampling binary vectors in topological order
  # while fixing one node at a time to compute intervention effects

  n <- ncol(DAGparams$DAG) # number of nodes

  interventionMatrix <- matrix(0, nrow=n, ncol=n)

  binarySample <- matrix(NA, sampleSize, n)
  nodeOrder <- DAGtoorder(DAGparams$DAG) # get topological order

  for (j in rev(nodeOrder)) { # start with outpoints etc
    parentNodes <- which(DAGparams$DAG[, j]==1)
    binarySample[, j] <- BinarySampleFromDAGcore(j, parentNodes, DAGparams, binarySample, reducedVarianceSampling)
  }

  binarySampleOriginal <- binarySample # we reuse the original sample to match across different MC estimates
  for (j in nodeOrder) { # start with the inpoints
    # use graph descent to get descendents, but would be more efficient to build the descendant matrix
    descendants <- NodeDescendantsFromDAG(j, DAGparams$DAG)
    if(length(descendants) > 0) {
      binarySample[, j] <- 1 # upregulate
      for (kk in intersect(rev(nodeOrder), descendants)) { # resample descendents in correct order
        parentNodes <- which(DAGparams$DAG[, kk]==1)
        binarySample[, kk] <- BinarySampleFromDAGcore(kk, parentNodes, DAGparams, binarySample, reducedVarianceSampling)
      }
      probs1 <- colMeans(binarySample)
      binarySample[, j] <- 0 # downregulate
      for (kk in intersect(rev(nodeOrder), descendants)) { # resample descendents in correct order
        parentNodes <- which(DAGparams$DAG[, kk]==1)
        binarySample[, kk] <- BinarySampleFromDAGcore(kk, parentNodes, DAGparams, binarySample, reducedVarianceSampling)
      }
      probs0 <- colMeans(binarySample)
      interventionMatrix[j, ] <- probs1 - probs0
      binarySample <- binarySampleOriginal # reset binary vectors
    } else {
      interventionMatrix[j, j] <- 1
    }

  }

  return(interventionMatrix)
}



# this function takes in an adjacency matrix and returns the permutation

DAGtoorder <- function(incidence) {
  n <- ncol(incidence) # number of nodes
  permy <- numeric(n) # to store the permutation
  m <- n # counter
  while (m>0) {
    topnodes <- which(colSums(incidence)==0) # find the outpoints
    incidence[topnodes, ] <- 0 # remove their edges
    incidence[cbind(topnodes, topnodes)] <- 1 # add a one to their columns so they are no longer counted
    l <- length(topnodes) # the number of outpoints
    m <- m - l
    permy[m + 1:l]<-topnodes
  }

  return(permy)
}

