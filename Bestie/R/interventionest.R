#' Exact estimation of intervention effects for a single DAG or a chain of sampled DAGs
#'
#' \code{DAGintervention} takes a DAG or a sampled chain of DAGs (for example from
#' the \code{\link[BiDAG]{partitionMCMC}} function of the BiDAG package) and computes the
#' intervention effect of each node on all others by exhaustively
#' examining all possible binary states. This is exponentially complex in the
#' number of variables which should therefore be limited to around 20 or fewer.
#'
#' @param incidences a single adjacency matrix of a list of adjacency matrices of
#' sampled DAGs, with entry [i,j] equal to 1 when a directed edge exists from
#' node i to node j
#' @param dataParams the data and parameters used to learn the DAGs derived from the
#' \code{\link[BiDAG]{scoreparameters}} function of the BiDAG package
#' @param sample logical indicating whether to sample the parameters of each node
#' from the posterior (TRUE, default) or to take the expectation (FALSE)
#'
#' @return a single matrix or a list of matrices containing the full set of
#' intervention effects for each input DAG. Entry [i,j] is the downstream
#' effect on node j of intervening on node i
#' (the difference observed at node j when setting node i to 1 and 0)
#'
#' @examples
#'
#' scoreParam <- BiDAG::scoreparameters(8, "bde", BiDAG::Asia)
#' causalmat <- DAGintervention(BiDAG::Asiamat, scoreParam)
#'
#' @export
#'
#' @seealso \code{\link[BiDAG]{scoreparameters}}

DAGintervention <- function(incidences, dataParams, sample = TRUE){
  # this wrapper takes in a chain of DAG, computes their parameters and returns all intervention effects

  if (!is.list(incidences)) { # turn to list internally
    incidences <- list(incidences)
  }

  n <- ncol(incidences[[1]]) # number of nodes in DAG
  
  if (n > 20) {
    warning("Exhaustive enumeration may not be feasible")
  }

  allBinaryVecs <- matrix(0, 2^n, n)
  for (ii in 1:2^n) {
    allBinaryVecs[ii, ] <- as.integer(intToBits(ii-1))[n:1]
  }

  numDAGs <- length(incidences)
  interventionMats <- vector("list", numDAGs) # to store the intervention effects

  for(kk in 1:numDAGs){
    DAGparams <- DAGparameters(incidences[[kk]], dataParams)

    DAGparamsInternal <- DAGparams
    if(sample==TRUE){ # then we take a sample of parameters from the posterior instead of taking the expectation
      DAGparamsInternal$pmeans <- SampleParameters(DAGparams)
    }

    allLogScores <- BinaryScoreAgainstDAG(DAGparamsInternal, allBinaryVecs)
    interventionMats[[kk]] <- InterventionEstimation(allLogScores, allBinaryVecs)
  }

  if (numDAGs == 1) { # turn back to single matrix
    interventionMats <- interventionMats[[1]]
  }

  return(interventionMats)
}



DAGinterventionparams <- function(DAGparams, sample = TRUE){
  # this wrapper takes in a DAG with parameters and returns all intervention effects

  n <- ncol(DAGparams$DAG) # number of nodes in DAG

  allBinaryVecs <- matrix(0, 2^n, n)
  for (ii in 1:2^n) {
    allBinaryVecs[ii, ] <- as.integer(intToBits(ii-1))[n:1]
  }

  DAGparamsInternal <- DAGparams
  if(sample==TRUE){ # then we take a sample of parameters from the posterior instead of taking the expectation
    DAGparamsInternal$pmeans <- SampleParameters(DAGparams)
  }

  allLogScores <- BinaryScoreAgainstDAG(DAGparamsInternal, allBinaryVecs)
  interventionMat <- InterventionEstimation(allLogScores, allBinaryVecs)

  return(interventionMat)
}



InterventionEstimation <- function(logScores, dataToScore){
  # estimate intervention effects from the scores of all binary vectors

  n <- ncol(dataToScore) # number of nodes

  interventionMatrix <- matrix(0, nrow=n, ncol=n)
  for (j in 1:n) {
    interventionMatrix[j,] <- InterventionEstimationCore(j, logScores, dataToScore)
  }

  return(interventionMatrix)
}



InterventionEstimationCore <- function(fixedNode, logScores, dataToScore){
  # core function involves setting one node to 0 and to 1 and removing its probability component

  upRows <- which(dataToScore[, fixedNode]==1)

  weightVecTemp <- rowSums(logScores[upRows, -fixedNode]) # probability of each vector with fixedNode set to one
  weightVec <- exp(weightVecTemp - max(weightVecTemp))
  weightVec <- weightVec/sum(weightVec)
  probs1 <- colSums(weightVec*dataToScore[upRows, ])

  downRows <- which(dataToScore[, fixedNode]==0)

  weightVecTemp <- rowSums(logScores[downRows, -fixedNode]) # probability of each vector with fixedNode set to zero
  weightVec <- exp(weightVecTemp - max(weightVecTemp))
  weightVec <- weightVec/sum(weightVec)
  probs0 <- colSums(weightVec*dataToScore[downRows, ])

  return(probs1-probs0)
}


