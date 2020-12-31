#----------------------------------------------------------------
##' Computes the SUR criterion associated to an equilibrium for a given \code{xnew} and a set of trajectories of objective functions
##' on a predefined grid.
##' @title SUR criterion for equilibria
##' @param idx is the index on the grid of the strategy evaluated
## ' #@param idx,x is the strategy evaluated (at least one should be provided). idx is the index on the grid (faster)
## ' #while x is the strategy value (vector of size dim)
##' @param model is a list of \code{nobj} \code{\link[DiceKriging]{km}} models
##' @param integcontrol is a list containing: \code{integ.pts}, a [\code{npts x dim}] matrix defining the grid,
##' \code{expanded.indices} a matrix containing the indices of the \code{integ.pts} on the grid and \code{n.s},
##' a \code{nobj} vector containing the number of strategies per player
##' @param Simu is a matrix of size [\code{npts x nsim*nobj}] containing the trajectories of the
##' objective functions (one column per trajectory,
##' first all the trajectories for obj1, then obj2, etc.)
##' @param precalc.data is a list of length \code{nobj} of precalculated data (based on kriging models at integration points)
##' for faster computation - computed if not provided
##' @param equilibrium equilibrium type: either "\code{NE}", "\code{KSE}", "\code{CKSE}" or "\code{NKSE}"
##' @param n.ynew is the number of \code{ynew} simulations (if not provided, equal to the number of trajectories)
##' @param cross if \code{TRUE}, all the combinations of trajectories are used (increases accuracy but also cost)
##' @param IS if \code{TRUE}, importance sampling is used for ynew
##' @param plot if \code{TRUE}, draws equilibria samples (should always be turned off)
## ' #@param cand.pts,J necessary if a filter has been applied (to avoid mismatch between idx and expanded.indices).
## ' #cand.pts is a [ncand x dim] matrix (without filter, equal to integ.pts) and J the vector of indices of cand.pts on the grid.
##' @param kweights kriging weights for \code{CKS} (TESTING)
##' @param Nadir optional vector of size \code{nobj}. Replaces the nadir point for \code{KSE}. If only a subset of values needs to be defined,
##' the other coordinates can be set to \code{Inf}.
##' @export
##' @seealso \code{\link[GPGame]{crit_PNash}} for an alternative infill criterion
##' @references
##' V. Picheny, M. Binois, A. Habbal (2016+), A Bayesian optimization approach to find Nash equilibria,
##' \emph{https://arxiv.org/abs/1611.02440}.
##' @importFrom stats cov
##' @importFrom GPareto checkPredict plotParetoEmp
##' @importFrom KrigInv predict_nobias_km computeQuickKrigcov precomputeUpdateData
##' @importFrom grDevices rainbow
##' @importFrom graphics pairs points
## ' @importFrom emoa is_dominated nondominated_points
##' @examples
##' \dontrun{
##' ##############################################
##' # 2 variables, 2 players
##' ##############################################
##' library(DiceKriging)
##' set.seed(42)
##'
##' # Objective function (R^2 -> R^2)
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
##' n.s <- rep(14, 2)
##' x.to.obj   <- c(1,2)
##' gridtype <- 'cartesian'
##' integcontrol <- generate_integ_pts(n.s=n.s, d=4, nobj=2, x.to.obj = x.to.obj, gridtype=gridtype)
##' integ.pts <- integcontrol$integ.pts
##' expanded.indices <- integcontrol$expanded.indices
##'
##' # Kriging models
##' n.init <- 11
##' design <- integ.pts[sample.int(n=nrow(integ.pts), size=n.init, replace=FALSE),]
##' response <- t(apply(design, 1, fun))
##' mf1 <- km(~., design = design, response = response[,1], lower=c(.1,.1))
##' mf2 <- km(~., design = design, response = response[,2], lower=c(.1,.1))
##' model <- list(mf1, mf2)
##'
##' # Conditional simulations
##' Simu <- t(Reduce(rbind, lapply(model, simulate, nsim=10, newdata=integ.pts, cond=TRUE,
##'                                    checkNames=FALSE, nugget.sim = 10^-8)))
##'
##' # Useful precalculations
##' library(KrigInv)
##' precalc.data <- lapply(model, FUN=KrigInv:::precomputeUpdateData, integration.points=integ.pts)
##'
##' # Compute criterion for all points on the grid
##' crit_grid <- lapply(X=1:prod(n.s), FUN=crit_SUR_Eq, model=model,
##'                     integcontrol=integcontrol, equilibrium = "NE",
##'                     Simu=Simu, precalc.data=precalc.data, n.ynew=10, IS=FALSE, cross=FALSE)
##' crit_grid <- unlist(crit_grid)
##'
##' # Draw contour of the criterion
##' filled.contour(seq(0, 1, length.out = n.s[1]), seq(0, 1, length.out = n.s[2]),
##'                matrix(pmax(0, crit_grid), n.s[1], n.s[2]), main = "SUR criterion",
##'                xlab = expression(x[1]), ylab = expression(x[2]), color = terrain.colors,
##'                plot.axes = {axis(1); axis(2);
##'                             points(design[,1], design[,2], pch = 21, bg = "white")
##'                            }
##' )
##' }
##'
crit_SUR_Eq <- function(idx, model, integcontrol, Simu, precalc.data=NULL, equilibrium,
                        n.ynew=NULL, cross=FALSE, IS=FALSE, plot=FALSE, kweights = NULL, Nadir = NULL){

  # if (is.null(cand.pts)) cand.pts <- integ.pts

  expanded.indices <- integcontrol$expanded.indices
  integ.pts <- integcontrol$integ.pts
  n.s <- integcontrol$n.s

  xnew <- matrix(as.numeric(integ.pts[idx,]), nrow=1)

  # if (!is.null(idx))  {
  #   xnew <- matrix(as.numeric(integ.pts[idx,]), nrow=1)
  # } else if (!is.null(x)) {
  #   xnew <- matrix(as.numeric(x), nrow=1)
  #   idx <- get_alternatives(s, i, expanded.indices)
  # } else {
  #   cat("At least one of the inputs idx and x should be provided \n")
  #   return(NA)
  # }

  nobj <- length(model)
  nsim <- ncol(Simu)/nobj
  nsimpts <- nrow(Simu)

  if (is.null(n.ynew)) n.ynew <- nsim

  if (!checkPredict(x=xnew, model=model)){

    # Precalculations if not provided
    if (is.null(precalc.data)) {
      precalc.data <- lapply(model, FUN=precomputeUpdateData, integration.points=integ.pts)
    }

    # Precalculs magiques pour la mise a jour rapide des trajectoires
    lambda <- matrix(NA, nsimpts, nobj)
    Ynew1  <- matrix(NA, nsim, nobj)
    Ynew2  <- matrix(NA, n.ynew, nobj)

    for (u in 1:nobj){
      prednew <- predict_nobias_km(model[[u]], newdata=data.frame(x=xnew), "UK", checkNames=FALSE)

      kn <- computeQuickKrigcov(model=model[[u]], integration.points=as.data.frame(x=integ.pts), X.new=data.frame(x=xnew),
                                precalc.data=precalc.data[[u]], F.newdata=prednew$F.newdata, c.newdata=prednew$c)
      kn_plus <- predict(model[[u]], data.frame(x=xnew), "UK", checkNames=FALSE)$sd^2
      lambda[,u] <- kn/kn_plus
      if (IS) {
        Ynew2[,u]  <- qnorm(p=sample(seq(1/(2*n.ynew),1, 1/n.ynew)), mean=prednew$mean, sd=prednew$sd)
      } else {
        Ynew2[,u]  <- rnorm(n=n.ynew, mean=prednew$mean, sd=prednew$sd)
      }
    }
    # Ynew1 is a [nsim x nobj] matrix - one line for each simulation
    # if (!is.null(J)) Ynew1  <- matrix(Simu[idx,], nsim, nobj)
    # else
    Ynew1  <- matrix(Simu[idx,], nsim, nobj)

    sorted <- !is.unsorted(expanded.indices[,nobj])
    if (plot) plot(NA, xlim=c(-200, 100), ylim=c(-50,-10))

    Gamma <- apply(Ynew2, 1, computeGamma, Simu=Simu, lambda=lambda, Ynew=Ynew1, n.s=n.s, kweights = kweights, Nadir=Nadir,
                   expanded.indices=expanded.indices, cross=cross, sorted=sorted, equilibrium = equilibrium, plot=plot)

    return(mean(Gamma, na.rm =TRUE))
  } else {
    return(NA)
  }
}
#----------------------------------------------------------------
## ' Computes Gamma hat for a given ynew and a set of trajectories
## ' @title Gamma calculations
## ' @param ynew is a [nobj] vector
## ' @param Simu is a list of length nobj containing matrices of size [nsim x npts]
## ' @param lambda is a [npts x nobj] matrix
## ' @param Ynew is a [nsim x nobj] matrix
## ' @param equilibrium equilibrium type
## ' @param n.s scalar of vector. If scalar, total nb of strategies (to be divided equally among players), otherwise nb of strategies per player.
## ' @param expanded.indices matrix containing the indices of the integ.pts on the grid
## ' @param cross if TRUE, all the combinations of trajectories are used
## ' @param sorted Boolean; if TRUE, the last column of expanded.indices is assumed to be sorted in increasing order. This provides a substantial efficiency gain.
## ' @param plot if TRUE, draws equilibria samples (should always be turned off)
## ' @param kweights kriging weights for \code{CKS} (TESTING)
computeGamma <- function(ynew, Simu, lambda, equilibrium, Ynew, n.s = NULL, kweights = NULL, Nadir = NULL,
                         expanded.indices = NULL, cross = FALSE, sorted = NULL, plot=FALSE){

  if (is.null(sorted)) {
    if (is.null(expanded.indices)) sorted <- FALSE
    else                           sorted <- !is.unsorted(expanded.indices[,nobj])
  }

  nobj <- length(ynew)
  nsim <- ncol(Simu)/nobj

  for (u in 1:nobj) {
    Simu[,(1:nsim)+(u-1)*nsim] <- Simu[,(1:nsim)+(u-1)*nsim] + tcrossprod(lambda[,u], rep(ynew[u], nsim) - Ynew[,u])
  }

  NE_simu_new <- getEquilibrium(Simu, equilibrium = equilibrium, nobj=nobj, n.s=n.s, expanded.indices=expanded.indices,
                                sorted=sorted, cross=cross, kweights = kweights, Nadir=Nadir)

  # Remove simulations without equilibrium
  NE_simu_new <- NE_simu_new[which(!is.na(NE_simu_new[,1])),, drop = FALSE]

  if (plot) {
    points(NE_simu_new[,1], NE_simu_new[,2], pch=sample.int(25,1), col=sample(rainbow(25),1))
  }

  if(is.null(dim(NE_simu_new)))  return(NA)
  else{
    result <- determinant(cov(NE_simu_new), logarithm = TRUE)[[1]][1]
    if(is.na(result)) return(NA)
    ## If rank of cov(NE_simu_new) inferior to nobj
    if(result == -Inf){
      result <- eigen(cov(NE_simu_new), only.values = T)
      result <- sum(log(result$values[which(result$values >0)]))
    }
    return(result)
  }
}
