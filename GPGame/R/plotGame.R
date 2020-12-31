##' Plot equilibrium for 2 objectives test problems with evaluations on a grid. The number of variables is not limited.
##' @title Visualisation of equilibrium solution in input/output space
##' @param fun name of the function considered
##' @param domain optional matrix for the bounds of the domain (for now [0,1]^d only), (two columns matrix with min and max)
##' @param n.grid number of divisions of the grid in each dimension (must correspond to \code{n.s} for Nash equilibriums)
##' @param graphs either \code{"design"}, \code{"objective"} or \code{"both"} (default) for which graph to display
##' @param x.to.obj,integcontrol see \code{\link[GPGame]{solve_game}} (for Nash equilibrium only)
##' @param equilibrium either "\code{NE}" for Nash, "\code{KSE}" for Kalai-Smoridinsky and "\code{NKSE}" for Nash-Kalai-Smoridinsky
##' @param fun.grid optional matrix containing the values of \code{fun} at \code{integ.pts}. Computed if not provided.
##' @param Nadir optional vector of size \code{nobj}. Replaces the nadir point for \code{KSE}. If only a subset of values needs to be defined, 
##' the other coordinates can be set to \code{Inf}.
##' @param ... further arguments to \code{fun}
##' @return list returned by invisible() with elements:
##' \itemize{
## '  \item \code{trueEq} equilibrium (list)
##'  \item \code{trueEqdesign} design corresponding to equilibrium value \code{trueEq}
##'  \item \code{trueEqPoff} corresponding values of the objective
##'  \item \code{trueParetoFront} Pareto front
##'  \item \code{response.grid}
##'  \item \code{integ.pts, expanded.indices}
##' }
##' @importFrom graphics lines
##' @export
## ' @details
## ' Options to plot Shadow and nadir points?
##' @examples
##' \dontrun{
##' library(GPareto)
##'
##' ## 2 variables
##' dom <- matrix(c(0,0,1,1),2)
##'
##' plotGameGrid("P1", domain = dom, n.grid = 51, equilibrium = "NE")
##' plotGameGrid("P1", domain = dom, n.grid = rep(31,2), equilibrium = "NE") ## As in the tests
##' plotGameGrid("P1", domain = dom, n.grid = 51, equilibrium = "KSE")
##' plotGameGrid("P1", domain = dom, n.grid = rep(31,2), equilibrium = "NKSE")
##' plotGameGrid("P1", graphs = "design", domain = dom, n.grid = rep(31,2), equilibrium = "NKSE")
##'
##' ## 4 variables
##' dom <- matrix(rep(c(0,1), each = 4), 4)
##' plotGameGrid("ZDT3", domain = dom, n.grid = 25, equilibrium = "NE", x.to.obj = c(1,1,2,2))
##'
##' }
##'
plotGameGrid <- function(fun=NULL, domain=NULL, n.grid, graphs = c("both", "design", "objective"), x.to.obj = NULL,
                         integcontrol = NULL, equilibrium = c("NE", "KSE", "CKSE", "NKSE"), fun.grid = NULL, Nadir = NULL, ...){

  if (is.null(fun) && is.null(fun.grid)) {
    cat("Either fun or fun.grid must be provided \n")
    return(NA)
  }

  expanded.indices <- integcontrol$expanded.indices
  integ.pts <- integcontrol$integ.pts

  equilibrium <- match.arg(equilibrium)
  graphs <- match.arg(graphs)

  if (is.null(domain) && is.null(integ.pts)) {
    cat("At least one of the following inputs must be provided: domain, integcontrol$integ.pts")
  }
  if (!is.null(integ.pts))    d <- ncol(integ.pts)
  else if (!is.null(domain))  d <- nrow(domain)

  if(equilibrium == "KSE" && is.null(integ.pts)){
    integ.pts <- expand.grid(sapply(1:d, function(i){seq(0, 1, length.out = n.grid)}, simplify = F))
  } else {
    if (is.null(integ.pts) || is.null(expanded.indices)){
      nobj <- max(length(unique(x.to.obj)), 2)
      if (length(n.grid)==1) n.grid <- rep(n.grid, nobj)
      res <- generate_integ_pts(n.s=n.grid, d = d, nobj = nobj, x.to.obj=x.to.obj, gridtype="cartesian")
      integ.pts <- res[[1]]
      expanded.indices <- res[[2]]
    }
  }
  n.integ.pts <- nrow(integ.pts)

  if (!is.null(fun.grid)) {
    if (nrow(fun.grid) != n.integ.pts) {
      cat("fun.grid inconsistent with either integ.pts or n.grid \n")
      return(NA)
    }
  }

  # Actual function values and Nash equilibrium
  if (is.null(fun.grid)) fun.grid <- t(apply(integ.pts, 1, fun, ... = ...))
  nobj <- ncol(fun.grid)
  trueEq <- getEquilibrium(Z=fun.grid, equilibrium = equilibrium, nobj=nobj, n.s=n.grid, return.design=TRUE,
                           expanded.indices = expanded.indices, sorted=TRUE, Nadir = Nadir)

  trueEqPoff <- trueEq[[1]]
  trueEqdesign <- integ.pts[trueEq[[2]],]
  # I.nd <- is_dominated(t(fun.grid))
  I.nd <- nonDom(fun.grid, return.idx=TRUE)
  trueParetoFront <- fun.grid[!I.nd,]

  # Plot actual problem and solution
  cols <- rep("black", n.integ.pts)
  lwds <- rep(1, n.integ.pts)
  pchs <- rep(1, n.integ.pts)

  ## Add Nash equilibrium for NKS
  if(equilibrium == "NKSE"){
    trueNEq <- getNashEquilibrium(Z=fun.grid, nobj=nobj, n.s=n.grid, return.design=TRUE,
                                  expanded.indices = expanded.indices, sorted=TRUE, cross=FALSE)
    cols[trueNEq[[2]]] <- "orange"
    lwds[trueNEq[[2]]] <- 2
    pchs[trueNEq[[2]]] <- 4
    Nadir <- apply(trueParetoFront, 2, min)
    Shadow <- apply(trueNEq[[1]], 2, max)
  }

  cols[!I.nd] <- "red"
  lwds[!I.nd] <- 2
  pchs[!I.nd] <- 2
  cols[trueEq[[2]]] <- "green"
  lwds[trueEq[[2]]] <- 3
  pchs[trueEq[[2]]] <- 3

  if(equilibrium == "KSE"){
    Shadow <- apply(trueParetoFront, 2, max)
    Nadir <- apply(trueParetoFront, 2, min)
  }

  if (nobj == 2 && graphs != "design") {
    # pdf(paste0(pbname, "_obj.pdf"), width=5, height=5)
    # par(mfrow=c(1,1), mar=c(3,3,2,1), mgp=c(2,1,0))
    plot(fun.grid[,1], fun.grid[,2], col=cols, lwd=lwds, pch=pchs,
         xlab=expression(y[1]), ylab=expression(y[2]), main="Objective space")

    if(equilibrium == "KSE" || equilibrium == "NKSE"){
      lines(x = c(Nadir[1], Shadow[1], Shadow[1], Nadir[1], Nadir[1], Shadow[1]),
            y = c(Nadir[2], Nadir[2], Shadow[2], Shadow[2], Nadir[2], Shadow[2]), lwd = 2, lty = 2, col = "cyan")
    }

    # dev.off()
  }

  if (d==2 && graphs != "objective") {
    # pdf(paste0(pbname, "_design.pdf"), width=5, height=5)
    # par(mfrow=c(1,1), mar=c(3,3,2,1), mgp=c(2,1,0))
    plot(integ.pts[,1], integ.pts[,2], col=cols, lwd=lwds, pch=pchs, xlim=c(0,1), ylim=c(0,1),
         xlab=expression(x[1]), ylab=expression(x[2]), main="Design space")
    # dev.off()
  } else {
    # pdf(paste0(pbname, "_design.pdf"), width=10, height=10)
    # par(mfrow=c(1,1), mar=c(3,3,2,1), mgp=c(2,1,0))
    if(graphs != "objective")
      pairs(integ.pts, col=cols, lwd=lwds, pch=pchs, xlim=c(0,1), ylim=c(0,1))
    # dev.off()
  }

  invisible(list(trueEqdesign = trueEqdesign, trueParetoFront = trueParetoFront, trueEqPoff = trueEqPoff, #trueEq = trueEq,
                 response.grid = fun.grid, integ.pts = integ.pts, expanded.indices = expanded.indices))


}

##' Plot equilibrium search result (2-objectives only)
##' @param res list returned by \code{\link[GPGame]{solve_game}}
##' @param equilibrium either "\code{NE}" for Nash, "\code{KSE}" for Kalai-Smoridinsky and "\code{NKSE}" for Nash-Kalai-Smoridinsky
##' @param add logical; if \code{TRUE} adds the first graphical output to an already existing plot; if \code{FALSE}, (default) starts a new plot
##' @param UQ_eq logical; should simulations of the equilibrium be displayed?
##' @param simus optional matrix of conditional simulation if \code{UQ_Eq} is \code{TRUE}
##' @param integcontrol list with \code{n.s} element (maybe n.s should be returned by solve_game). See \code{\link[GPGame]{solve_game}}.
##' @param simucontrol optional list for handling conditional simulations. See \code{\link[GPGame]{solve_game}}.
##' @param Nadir optional vector of size \code{nobj}. Replaces the nadir point for \code{KSE}. If only a subset of values needs to be defined, 
##' the other coordinates can be set to \code{Inf}.
##' @param ncores number of CPU available (> 1 makes mean parallel \code{TRUE})
##' @export
##' @examples
##' \dontrun{
##' library(GPareto)
##' library(parallel)
##'
##' # Turn off on Windows
##' parallel <- FALSE # TRUE
##' ncores <- 1
##' if(parallel) ncores <- detectCores()
##' cov.reestim <- TRUE
##' n.sim <- 20
##' n.ynew <- 20
##' IS <- TRUE
##' set.seed(1)
##'
##' pb <- "P1" # 'P1' 'PDE' 'Diff'
##' fun <- P1
##'
##' equilibrium = "NE"
##'
##' d <- 2
##' nobj <- 2
##' n.init <- 20
##' n.ite <- 4
##' model.trend <- ~1
##' n.s <- rep(31, 2) #31
##' x.to.obj   <- c(1,2)
##' gridtype <- 'cartesian'
##' nsimPoints <- 800
##' ncandPoints <- 200
##' sur_window_filter <- NULL
##' sur_pnash_filter  <- NULL
##' Pnash_only_filter <- NULL
##' res <- solve_game(fun, equilibrium = equilibrium, crit = "sur", model = NULL, n.init=n.init,
##'   n.ite = n.ite, nobj=nobj, x.to.obj = x.to.obj, integcontrol=list(n.s=n.s, gridtype=gridtype),
##'   simucontrol=list(n.ynew=n.ynew, n.sim=n.sim, IS=IS), ncores = ncores, d = d,
##'   filtercontrol=list(filter=sur_window_filter, nsimPoints=nsimPoints, ncandPoints=ncandPoints),
##'   kmcontrol=list(model.trend=model.trend), trace=3,
##'   seed=1)
##' plotGame(res, equilibrium = equilibrium)
##'
##' dom <- matrix(c(0,0,1,1),2)
##' plotGameGrid("P1", graphs = "objective", domain = dom, n.grid = 51, equilibrium = equilibrium)
##' plotGame(res, equilibrium = equilibrium, add = TRUE)
##'
##'
##' }
plotGame <- function(res, equilibrium = "NE", add = FALSE, UQ_eq = TRUE, simus = NULL, integcontrol = NULL, simucontrol = NULL, Nadir = NULL, ncores = 1){

  if (length(res$model)>2) {
    cat("plotGame works only for two players/objectives \n")
    return(NA)
  }

  if (UQ_eq) {
    if (is.null(simus)) {
      # Check simucontrol
      if (!is.null(simucontrol$n.sim)) n.sim <- simucontrol$n.sim else n.sim <- 10
      if (!is.null(simucontrol$IS)) IS <- simucontrol$IS else IS <- TRUE
      simus <- try(t(Reduce(rbind, mclapply(res$model, simulate, nsim = n.sim, newdata = res$integcontrol$integ.pts, cond=TRUE,
                                            checkNames = FALSE, nugget.sim = 10^-8, mc.cores = ncores))))
      if (typeof(simus) == "character") {
        cat("Conditional simulations failed - maybe there are too many integration points \n Correct or set UQ_eq = FALSE \n")
        Eq_simu <- NULL
        return(NA)
      }
    }
    if (is.null(integcontrol$n.s)) integcontrol$n.s <- apply(res$integcontrol$expanded.indices, 2, max)
    Eq_simu <- getEquilibrium(simus, equilibrium = equilibrium, nobj = length(res$model), n.s = integcontrol$n.s,
                              expanded.indices = res$integcontrol$expanded.indices, sorted = TRUE, cross = FALSE, Nadir = Nadir)
  } else {
    Eq_simu <- NULL
  }

  xlim <- range(c(res$model[[1]]@y, Eq_simu[,1], res$predEq[[length(res$predEq)]]$NEPoff[1]))
  ylim <- range(c(res$model[[2]]@y, Eq_simu[,2], res$predEq[[length(res$predEq)]]$NEPoff[2]))

  if(!add){
    plot(res$model[[1]]@y, res$model[[2]]@y, xlab = expression(f[1]), ylab = expression(f[2]), xlim=xlim, ylim=ylim, pch = 20, col = 'blue')
  }else{
    points(res$model[[1]]@y, res$model[[2]]@y, pch = 20, col = 'blue')
  }
  points(Eq_simu, col = "violet", pch = 3, lwd = 2)
  if(!is.null(res$predEq)) points(res$predEq[[length(res$predEq)]]$NEPoff, pch = 2, lwd = 2, col="green")
}
