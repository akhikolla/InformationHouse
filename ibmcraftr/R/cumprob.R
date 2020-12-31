#' Calculate cumulative probabilities for state transitions.
#'
#' This function takes in a vector of probabilities of states transitions and calculate the probability of staying in the original state and output the cumulative probabilities for all possibilities.
#'
#' @param probs A numeric vector of the probabilities of transition to states.
#' @param actual A logical value, if TRUE, will calculate actual cumulative probabilities which may surpass 1!.
#' @return A numeric vector of cumulative probabilites inclusive of the probability of having the same state in the next timestep.
#' @examples
#' cumprob(c(.2,.2,.9))
#' cumprob(c(.2,.2,.9), actual=TRUE)
#' cumprob(c(.2,.2,.2))
#'
#' @export



#####function####
cumprob <- function(probs, actual=FALSE){
  #
  #probability of stay will be calculated automatically
  #actual, default FALSE, will calculate a cumulative probability within 0-1
  #
  #probs <- x

  compliments <- 1-probs

  #stay <- compliments[1]-probs[2]
  stay <- abs(compliments[1] - sum(probs[-1]))
  result <- cumsum(append(stay, probs))

  if(actual){
    result
  } else {
    result/result[length(result)]
  }
}

