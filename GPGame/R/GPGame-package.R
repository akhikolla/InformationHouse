##' Sequential strategies for finding game equilibria in a black-box setting (expensive pay-off evaluations, no derivatives).
##' Handles noiseless or noisy evaluations. Two acquisition functions are available. Graphical outputs can be generated automatically.
##' @title Package GPGame
##' @author Victor Picheny, Mickael Binois
##' @name GPGame
##' @references
##' V. Picheny, M. Binois, A. Habbal (2016+), A Bayesian Optimization approach to find Nash equilibria,
##' \emph{https://arxiv.org/abs/1611.02440}.
##' 
##' M. Binois, V. Picheny, A. Habbal, "The Kalai-Smorodinski solution for many-objective Bayesian optimization", 
##' NIPS BayesOpt workshop, December 2017, Long Beach, USA,
##' \emph{https://bayesopt.github.io/papers/2017/28.pdf}.
##' 
##' @details
##' Important functions: \cr
##' \code{\link[GPGame]{solve_game}} \cr
##' \code{\link[GPGame]{plotGame}} \cr
##'
##' @seealso \code{\link[DiceKriging]{DiceKriging}}, \code{\link[DiceOptim]{DiceOptim}}, \code{\link[KrigInv]{KrigInv}}, \code{\link[GPareto]{GPareto}}
##' @examples
##' \dontrun{
##' # To use parallel computation (turn off on Windows)
##' library(parallel)
##' parallel <- FALSE # TRUE # 
##' if(parallel) ncores <- detectCores() else ncores <- 1
##'
##' ##############################################
##' # 2 variables, 2 players, Nash equilibrium
##' # Player 1 (P1) wants to minimize fun1 and player 2 (P2) fun2
##' # P1 chooses x2 and P2 x2
##'
##' ##############################################
##' # First, define objective function fun: (x1,x2) -> (fun1,fun2)
##' fun <- function (x)
##' {
##'   if (is.null(dim(x)))    x <- matrix(x, nrow = 1)
##'   b1 <- 15 * x[, 1] - 5
##'   b2 <- 15 * x[, 2]
##'   return(cbind((b2 - 5.1*(b1/(2*pi))^2 + 5/pi*b1 - 6)^2 + 10*((1 - 1/(8*pi)) * cos(b1) + 1),
##'                -sqrt((10.5 - b1)*(b1 + 5.5)*(b2 + 0.5)) - 1/30*(b2 - 5.1*(b1/(2*pi))^2 - 6)^2-
##'                 1/3 * ((1 - 1/(8 * pi)) * cos(b1) + 1)))
##' }
##'
##' ##############################################
##' # x.to.obj indicates that P1 chooses x1 and P2 chooses x2
##' x.to.obj   <- c(1,2)
##'
##' ##############################################
##' # Define a discretization of the problem: each player can choose between 21 strategies
##' # The ensemble of combined strategies is a 21x21 cartesian grid
##'
##' # n.s is the number of strategies (vector)
##' n.s <- rep(21, 2)
##' # gridtype is the type of discretization
##' gridtype <- 'cartesian'
##'
##' integcontrol <- list(n.s=n.s, gridtype=gridtype)
##'
##' ##############################################
##' # Run solver with 6 initial points, 14 iterations
##' n.init <- 6 # number of initial points (space-filling)
##' n.ite <- 14 # number of iterations (sequential infill points)
##'
##' res <- solve_game(fun, equilibrium = "NE", crit = "sur", n.init=n.init, n.ite=n.ite,
##'                   d = 2, nobj=2, x.to.obj = x.to.obj, integcontrol=integcontrol,
##'                   ncores = ncores, trace=1, seed=1)
##'
##' ##############################################
##' # Get estimated equilibrium and corresponding pay-off
##' NE <- res$Eq.design
##' Poff <- res$Eq.poff
##'
##' ##############################################
##' # Draw results
##' plotGame(res)
##'
##' ##############################################
##' # See solve_game for other examples
##' ##############################################
##' }
NULL
