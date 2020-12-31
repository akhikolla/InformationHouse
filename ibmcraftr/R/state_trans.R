#' Make state transitions.
#'
#' Take in the matrix of the states of synthetic population (created by \code{syn_pop} function)
#' and calculate the transitions from one state to other state(s) using the transition rate(s).
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
#' state_trans(1,2,.1,pop)
#' state_trans(1,4,100,pop)
#'
#' @export
#' @importFrom stats runif

state_trans <- function(origin, new.states, params, s.matrix){
  #origin   #single number
  #new.states  #a vector of length n (to index the matrix)
  #params #a vector of length m (to calculate the probabilities)
  #s.matrix  #state.matrix #a matrix cut from the data frame

  #dimension check
  if(ncol(s.matrix) <  max(c(origin, new.states))) stop("no such states in the input matrix") #stop if the dim requested is higher than input matrix

  origin_v <- s.matrix[,origin] #initializing a new vector for calculation
  lo <- length(origin_v) #length of origin
  org.s.matrix <- s.matrix    #keeping the original matrix ??? for what?

  #cummulative probability
  #probs <- 1-exp(-params*1) # calc probs from rates

  cum_probs <- cumprob(params) #cumprob is a seperate function to calculate the cumulative probabilities

  last_prob <- cum_probs[1]

  for(i in new.states){
    probs_for <- cum_probs[which(new.states==i)+1] #calculating probs_for for transition
    rand <- runif(lo)

    s.matrix[,i] <- s.matrix[,i]+(s.matrix[,origin]*(rand<probs_for)*(rand>last_prob)) #origin is used here since ??
    s.matrix[,origin] <- s.matrix[,origin]-(s.matrix[,origin]*(rand<probs_for)*(rand>last_prob))
    #s.matrix[,-c(new.states,origin)] <- 0 #might not do this afterall!
    last_prob <- probs_for
  }
  s.matrix-org.s.matrix
}
