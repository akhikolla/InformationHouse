#' Make state transitions using Rcpp.
#'
#' Take in the matrix of the states of synthetic population (created by \code{syn_pop} function)
#' and calculate the transitions from one state to other state(s) using the transition probabilities [not rate(s)].
#' The major difference from the R alone version was that instead of having the transition rate(s),
#' transition probabilities are used. These probabilities will thus be calculated with another function.
#'
#' @param origin A number which represents the column index \code{s.matrix} you want to do the transition from
#' @param new.states A numeric vector or a number which represents the column index \code{s.matrix} you want
#'     as the destination(s) for the transition
#' @param params A numeric vector of similar length to \code{new.states} which serves as the transition rate(s)
#' @param s.matrix A state matrix created from \code{syn_pop} function
#' @return A transition matrix of the same dimension as \code{s.matrix}. -1 indicates that the individual has left
#'     the corresponding state. +1 indicates that the individual has become the corresponding state.
#' @examples
#' pop <- syn_pop(c(19,1,0,0))
#' stRCPP(1,2,.1,pop)
#'
#' @useDynLib ibmcraftr
#'
#' @export
#' @importFrom stats runif
#' @importFrom Rcpp sourceCpp


#sourceCpp("D:/OneDrive/Rcpp/stRCPP.cpp") #same version as in the package folder
stRCPP <- function(origin, new.states, params, s.matrix){
  #origin   #single number
  #new.states  #a vector of length n (to index the matrix)
  #params #a vector of length m (to calculate the probabilities)
  #s.matrix  #state.matrix #a matrix cut from the data frame

  #dimension check
  if(ncol(s.matrix) <  max(c(origin, new.states))) stop("no such states in the input matrix") #stop if the dim requested is higher than input matrix

  org.s.matrix <- s.matrix    #keeping the original matrix for calculating the movement

  #probs <- 1-exp(-params*1) # calc probs from rates #this will be done in the "run_state_trans" function

  #cummulative probability # a seperate function, cumprob, is now used
  cum_probs <- cumprob(params)

  #load and run the Rcpp codes here
  s.matrix <- stateT(origin, new.states, cum_probs, s.matrix)

  s.matrix-org.s.matrix
}
