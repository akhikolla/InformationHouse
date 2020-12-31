## ' Probability for a strategy of being a Nash Equilibrium by simulation
## ' @title Simulated Pnash probability
## ' @param s strategy considered
## ' @param integcontrol is a list containing: \code{integ.pts}, a [npts x dim] matrix defining the grid,
## ' \code{expanded.indices} a matrix containing the indices of the integ.pts on the grid and \code{n.s},
## ' a nobj vector containting the number of strategies per player
## ' @param model list of km
## ' @param nsim number of conditional simulations for computation
## ' @param eps numerical jitter
## ' @export
## ' @examples
## ' \dontrun{
## ' Pnashs <- apply(expanded.indices, 1, Pnash, expanded.indices = expanded.indices,
## '                 integ.pts = integ.pts, model = model, nsim = 100)
## ' }
Pnash <- function(s, integcontrol, model, nsim, eps = 1e-6){
  if(is.null(eps))
    eps <- 1e-6

  expanded.indices <- integcontrol$expanded.indices
  integ.pts <- integcontrol$integ.pts

  ## Loop over objectives
  p <- 1
  for(i in 1:length(s)){
    idx <- get_alternatives(s, i, expanded.indices)
    preds <- predict(model[[i]], newdata = integ.pts[idx,], type = 'UK', checkNames = F, cov.compute = TRUE, se.compute = FALSE)
    sims <- mvrnorm(n = nsim, mu = preds$mean, Sigma = preds$cov + diag(eps, length(preds$mean)))
    sisBest <- apply(sims, 1, function(x) min(x) == x[s[i]])
    p <- p * sum(sisBest)/nsim
  }
  return(p)
}

## ' Probability for a stragy of being a Nash Equilibrium with mnormt
## ' @title Exact Pnash probability
## ' @param s strategy considered
## ' @param integcontrol is a list containing: \code{integ.pts}, a [npts x dim] matrix defining the grid,
## ' \code{expanded.indices} a matrix containing the indices of the integ.pts on the grid and \code{n.s},
## ' a nobj vector containting the number of strategies per player
## ' @param model list of km
## ' @param eps numerical jitter
## ' @examples
## ' \dontrun{
## ' Pnashs2 <- apply(expanded.indices, 1, Pnash_exact, expanded.indices = expanded.indices,
## '                  integ.pts = integ.pts, model = model)
## ' }
## ' @importFrom mnormt pmnorm
## ' @importFrom mvtnorm pmvnorm
## ' @details inspired from DiceOptim qEI function
## ' @export
Pnash_exact <- function(s, integcontrol, model, eps = 1e-8){
  p <- 1

  expanded.indices <- integcontrol$expanded.indices
  integ.pts <- integcontrol$integ.pts

  ## Loop over objectives
  for(i in 1:length(s)){
    idx <- get_alternatives(s, i, expanded.indices)

    preds <- predict(model[[i]], newdata = integ.pts[idx,], type = 'UK', checkNames = F, cov.compute = TRUE, se.compute = FALSE)

    mu_k <- preds$mean[s[i]] - preds$mean[-s[i]]
    Sigma_k <- covZk(preds$cov, s[i])
    Sigma_k <- Sigma_k[-s[i], -s[i]]

    if(length(idx) < 20){
      pNashsi <- pmnorm(x = -mu_k, varcov = Sigma_k + diag(eps, length(preds$mean) - 1)) # Error if more than 20 dimensions
    }else{
      pNashsi <- pmvnorm(upper = -mu_k, sigma = Sigma_k + diag(eps, length(preds$mean) - 1))
    }
    p <- p * pNashsi[1]
  }
  return(p)
}

## Just add argument i for compatibility with mclapply
Pnash_wrap <- function(i, s, integcontrol, type = 'simu', model, control = list(nsim = 100, eps = 1e-6)){
  if(is.null(nrow(s)))
    s <- matrix(s, nrow = 1)
  if(type == 'simu')
    return(Pnash(s = s[i,], integcontrol=integcontrol, model = model, nsim = control$nsim, eps = control$eps))
  if(type == 'exact')
    return(Pnash_exact(s = s[i,], integcontrol=integcontrol, model = model, eps = control$eps))
}

##' Acquisition function for solving game problems based on the probability for a strategy of being a Nash Equilibrium.
##' The probability can be computed exactly using the mutivariate Gaussian CDF (\code{mnormt}, \code{pmvnorm}) or by Monte Carlo.
##' @title Probability for a strategy of being a Nash Equilibrium
##' @param idx is the index on the grid of the strategy evaluated
##' @param integcontrol is a list containing: \code{integ.pts}, a [\code{npts x dim}] matrix defining the grid,
##' \code{expanded.indices} a matrix containing the indices of the \code{integ.pts} on the grid and \code{n.s},
##' a \code{nobj} vector containting the number of strategies per player
##' @param type '\code{exact}' or '\code{simu}'
##' @param model is a list of nobj \code{\link[DiceKriging]{km}} models
##' @param control list with slots \code{nsim} (number of conditional simulations for computation) and \code{eps}
##' @param ncores \code{\link[parallel]{mclapply}} is used if \code{> 1} for parallel evaluation
##' @param eps numerical jitter for stability
##' @export
##' @seealso \code{\link[GPGame]{crit_SUR_Eq}} for an alternative infill criterion
##' @references
##' V. Picheny, M. Binois, A. Habbal (2016+), A Bayesian optimization approach to find Nash equilibria,
##' \emph{https://arxiv.org/abs/1611.02440}.
##' @examples
##' \dontrun{
##' ##############################################
##' # Example 1: 2 variables, 2 players, no filter
##' ##############################################
##' library(DiceKriging)
##' set.seed(42)
##'
##' # Define objective function (R^2 -> R^2)
##' fun <- function (x)
##' {
##'   if (is.null(dim(x)))    x <- matrix(x, nrow = 1)
##'  b1 <- 15 * x[, 1] - 5
##'  b2 <- 15 * x[, 2]
##'  return(cbind((b2 - 5.1*(b1/(2*pi))^2 + 5/pi*b1 - 6)^2 + 10*((1 - 1/(8*pi)) * cos(b1) + 1),
##'                -sqrt((10.5 - b1)*(b1 + 5.5)*(b2 + 0.5)) - 1/30*(b2 - 5.1*(b1/(2*pi))^2 - 6)^2-
##'                 1/3 * ((1 - 1/(8 * pi)) * cos(b1) + 1)))
##' }
##'
##' # Grid definition
##' n.s <- rep(11, 2)
##' x.to.obj   <- c(1,2)
##' gridtype <- 'cartesian'
##' integcontrol <- generate_integ_pts(n.s=n.s, d=2, nobj=2, x.to.obj = x.to.obj, gridtype=gridtype)
##'
##' test.grid <- integcontrol$integ.pts
##' expanded.indices <- integcontrol$expanded.indices
##' n.init <- 11
##' design <- test.grid[sample.int(n=nrow(test.grid), size=n.init, replace=FALSE),]
##' response <- t(apply(design, 1, fun))
##' mf1 <- km(~., design = design, response = response[,1], lower=c(.1,.1))
##' mf2 <- km(~., design = design, response = response[,2], lower=c(.1,.1))
##' model <- list(mf1, mf2)
##'
##' crit_sim <- crit_PNash(idx=1:nrow(test.grid), integcontrol=integcontrol,
##'                        type = "simu", model=model, control = list(nsim = 100))
##' crit_ex <- crit_PNash(idx=1:nrow(test.grid), integcontrol=integcontrol, type = "exact", model=model)
##'
##' filled.contour(seq(0, 1, length.out = n.s[1]), seq(0, 1, length.out = n.s[2]), zlim = c(0, 0.7),
##'                matrix(pmax(0, crit_sim), n.s[1], n.s[2]), main = "Pnash criterion (MC)",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors,
##'                plot.axes = {axis(1); axis(2);
##'                             points(design[,1], design[,2], pch = 21, bg = "white")
##'                            }
##' )
##'
##' filled.contour(seq(0, 1, length.out = n.s[1]), seq(0, 1, length.out = n.s[2]), zlim = c(0, 0.7),
##'                matrix(pmax(0, crit_ex), n.s[1], n.s[2]), main = "Pnash criterion (exact)",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors,
##'                plot.axes = {axis(1); axis(2);
##'                             points(design[,1], design[,2], pch = 21, bg = "white")
##'                            }
##' )
##' }
##' @importFrom mnormt pmnorm
##' @importFrom mvtnorm pmvnorm
##'
crit_PNash <- function(idx, integcontrol, type = 'simu', model, ncores = 1, control = list(nsim = 100, eps = 1e-6)){

  if (is.null(control$nsim)) control$nsim <- 100
  if (is.null(control$eps)) control$eps <- 1e-6

  # if(is.null(nrow(x)))
  #   x <- matrix(x, nrow = 1)
  #
  # if (is.null(idx) && is.null(x)) {
  #   cat("At least one of the inputs idx and x should be provided \n")
  #   return(NA)
  # }

  s <- integcontrol$expanded.indices[idx,,drop=FALSE]

  if(ncores > 1)
    return(unlist(mclapply(X = 1:nrow(s), s = s, FUN = Pnash_wrap, model = model, integcontrol=integcontrol,
                           mc.cores = ncores, type = type, control = control)))

  return(apply(matrix(1:nrow(s), ncol = 1), 1, Pnash_wrap, s = s, integcontrol=integcontrol,
               type = type, model = model, control = control))
}

## ' Return the designs s[-i]
## ' @param s .
## ' @param i .
## ' @param expanded.indices .
## ' @export
get_alternatives <- function(s, i, expanded.indices){
  return(which(apply(expanded.indices[,-i, drop = FALSE], 1, function(x) all(x == s[-i]))))
}

covZk <- function(sigma,index){
  result <- sigma
  q <- nrow(sigma)
  I <- (1:q)[-index]
  result[I,index] <- result[index,I] <- sigma[index,index] - sigma[index,I]
  result[index,I] <- sigma[index,index] - sigma[I,index]
  M <- crossprod(sigma[index,I,drop=FALSE], rep(1,q-1))
  result[I,I] <- sigma[I,I] + sigma[index,index] - M - t(M)
  return(result)
}
