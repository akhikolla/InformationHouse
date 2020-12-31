#' Create a synthetic population having several states.
#'
#' Populate a matrix in which columns represent the states of the individuals and rows represent the individuals.
#'
#' @param states A numeric vector with each element representing the number of individuals in a particular state
#'     its index corresponds to.
#' @param shuffle A logical value to enable shuffling of the individuals (rows) in the resulting matrix.
#' @return A matrix of 0s, and 1s. The rows representing the individuals and the columns representing the states
#'     the individuals are in
#' @examples
#' syn_pop(c(3,2,1))
#' syn_pop(c(0,0,1,5), shuffle=TRUE)
#'
#' @export

syn_pop <- function(states, shuffle=FALSE){ #states is the vector variable, each element represent the number of individuals belonging to that state (indexed)
  #this is assuming that an individual can have only one state at a single timestep

  total.pop <- sum(states) #total population size

  result <- c(rep(1,states[1]), rep(0,total.pop - states[1])) #to store the resulting matrix, columns will be binded here in the following for loop


  for(i in 2:length(states)){
    tmp <- c(rep(0,sum(states[1:(i-1)])), rep(1, states[i]) , rep(0,total.pop - sum(states[1:i]) )) #tmp to store the result before cbinding

    result <- cbind(result, tmp, deparse.level = 0)
  }

  if(shuffle==T){ result[sample(nrow(result)),]}
  else result
}
