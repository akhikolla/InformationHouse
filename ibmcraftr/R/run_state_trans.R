#' Run state_trans function over a given number of timesteps.
#'
#' Organize population data and transition parameters to run state_trans function over the given number of timesteps.
#'
#' @param timesteps A numeric scalar based on which the state_trans function will run for that specific no. of \code{timesteps} and accumulate the results.
#' @param param A list of lists. Each low-level list must contain transition parameters required by the \code{state_trans} function.
#' @param pop A state matrix created from \code{syn_pop} function. This matrix represents the states of the population.
#' @param transient A character vector. Each element must include formula(e)/expression(s) to evaluate dynamic parameters after each timestep.
#' @param useC A logical value, which is TRUE by default, will run \code{state_transition} function written in RCPP, \code{stRCPP}.
#' @return A summary matrix of the states all individuals in the population are in.
#' @examples
#' pop <- syn_pop(c(19,1,0,0,0)) #synthesizing population
#' b <- 2 #effective contact rate
#' param <- list(
#' list(1,c(2,5),c(NA,.1)), #transition from state 1 to 2 using FOI lambda
#' list(2,3,100), #transition from state 2 to 3,
#' list(3,4,100)  #the 3rd term ensures the transition to the next stage
#' )
#'
#' timesteps <- 10
#' transient <- c("param[[1]][[3]][1] <- rate2prob(b*sum(pop[,2],pop[,3])/sum(pop))")
#' eval(parse(text=transient))
#'
#' run_state_trans(timesteps, param, pop, transient)
#' run_state_trans(timesteps, param, pop, transient, useC = FALSE)
#'
#' @export
#'


run_state_trans <- function(timesteps, param, pop, transient = "", useC = TRUE){
  Matrix.List <- list() #master matrix list initiazilation, to store the transition values each timestep
  sim.table <- matrix(NA, timesteps, ncol(pop))  #table to record the summaries each time step
  if(useC){
    for(i in 1:timesteps){
      for(j in 1:length(param)){
        Matrix.List[[j]] <- stRCPP(param[[j]][[1]], param[[j]][[2]], param[[j]][[3]], pop)
      }
      pop <- Reduce('+', Matrix.List) + pop #Population after transition

      eval(parse(text=transient))
      sim.table[i,] <- colSums(pop) #getting summaries of the population
    }
  }
  else {
    for(i in 1:timesteps){
      for(j in 1:length(param)){
        Matrix.List[[j]] <- state_trans(param[[j]][[1]], param[[j]][[2]], param[[j]][[3]], pop)
      }
      pop <- Reduce('+', Matrix.List) + pop #Population after transition

      eval(parse(text=transient))
      sim.table[i,] <- colSums(pop) #getting summaries of the population
    }
  }

  sim.table
}



