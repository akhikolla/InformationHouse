#' Main function to solve games.
##' @title Main solver
##' @details If \code{noise.var="given_by_fn"}, \code{fn} returns a list of two vectors, the first being the objective functions and the second
##' the corresponding noise variances.
##'
##' \code{integcontrol} controls the way the design space is discretized. One can directly provide a set of points \code{integ.pts} with
##' corresponding indices \code{expanded.indices} (for \code{NE}). Otherwise, the points are generated according to the number of strategies \code{n.s}.
##' If \code{n.s} is a scalar, it corresponds to the total number of strategies (to be divided equally among players),
##' otherwise it corresponds to the nb of strategies per player. In addition, one may choose the type of discretization with \code{gridtype}.
##' Options are '\code{lhs}' or '\code{cartesian}'. Finally, \code{lb} and \code{ub} are vectors specifying the bounds for the design variables.
##' By default the design space is \code{[0,1]^d}.
##' A \code{renew} slot is available, if \code{TRUE}, then \code{integ.pts} are changed at each iteration. Available only for \code{KSE} and \code{CKSE}.
##' For \code{CKSE}, setting the slot \code{kweights=TRUE} allows to increase the number of integration points, with \code{nsamp} 
##' (default to \code{1e4}) virtual simulation points.
##'
##' \code{simucontrol} controls options on conditional GP simulations. Options are \code{IS}: if \code{TRUE}, importance sampling is used for \code{ynew};
##' \code{n.ynew} number of samples of \code{Y(x_{n+1})} and \code{n.sim} number of sample path generated.
##'
##' \code{filtercontrol} controls filtering options. \code{filter} sets how to select a subset of simulation and candidate points,
##' either either a single value or a vector of two to use different filters for simulation and candidate points.
##' Possible values are '\code{window}', '\code{Pnash}' (for \code{NE}), '\code{PND}' (probability of non domination), '\code{none}'.
##' \code{nsimPoints} and \code{ncandPoints} set the maximum number of simulation/candidate points wanted
##' (use with filter '\code{Pnash}' for now). Default values are \code{800} and \code{200}, resp.
##' \code{randomFilter} (\code{TRUE} by default except for filter \code{window}) sets whereas the filter acts randomly or deterministically.
##' For more than 3 objectives, \code{PND} is estimated by sampling; the number of samples is controled by \code{nsamp} (default to \code{max(20, 5 * nobj)}).
##'
##' \code{kmcontrol} Options for handling \code{nobj} \code{\link[DiceKriging]{km}} models.
##' \code{cov.reestim} (Boolean, \code{TRUE} by default) specifies if the kriging hyperparameters
##' should be re-estimated at each iteration,
##'
##' \code{returncontrol} sets options for the last iterations and what is returned by the algorithm.
##' \code{track.Eq} allows to estimate the equilibrium at each iteration; options are '\code{none}' to do nothing,
##' "\code{mean}" (default) to compute the equilibrium of the prediction mean (all candidates),
##'  "\code{empirical}" (for \code{KSE}) and "\code{pex}"/"\code{psim}" (\code{NE} only)
##' for using \code{Pnash} estimate (along with mean estimate, on integ.pts only, NOT reestimated if \code{filter.simu} or \code{crit} is \code{Pnash}).
##' The boolean \code{force.exploit.last} (default to \code{TRUE}) allows to evaluate the equilibrium on the predictive mean - if not already evaluated -
##' instead of using \code{crit} (i.e., \code{sur}) for \code{KSE} and \code{CKSE}.
##'
##' @param fun fonction with vectorial output
##' @param equilibrium either '\code{NE}', '\code{KSE}', '\code{CKSE}' or '\code{NKSE}' for Nash / Kalai-Smorodinsky / Copula-Kalai-Smorodinsky
##' / Nash-Kalai-Smorodinsky equilibria
##' @param crit '\code{sur}' (default) is available for all equilibria, '\code{psim}' and '\code{pex}' are available for Nash
##' @param model list of \code{\link[DiceKriging]{km}} models
##' @param n.init number of points of the initial design of experiments if no model is given
##' @param n.ite number of iterations of sequential optimization
##' @param d variable dimension
##' @param nobj number of objectives (players)
##' @param x.to.obj for \code{NE} and \code{NKSE}, which variables for which objective
##' @param noise.var noise variance. Either a scalar (same noise for all objectives), a vector (constant noise, different for each objective),
##' a function (type closure) with vectorial output (variable noise, different for each objective) or \code{"given_by_fn"}, see Details.
##' If not provided, \code{noise.var} is taken as the average of \code{model@noise.var}.
##' @param integcontrol optional list for handling integration points. See Details.
##' @param simucontrol optional list for handling conditional simulations. See Details.
##' @param filtercontrol optional list for handling filters. See Details.
##' @param kmcontrol optional list for handling \code{\link[DiceKriging]{km}} models. See Details.
##' @param returncontrol optional list for choosing return options. See Details.
##' @param ncores number of CPU available (> 1 makes mean parallel \code{TRUE})
##' @param trace controls the level of printing: \code{0} (no printing), \code{1} (minimal printing), \code{3} (detailed printing)
##' @param seed to fix the random variable generator
##' @param Nadir optional vector of size \code{nobj}. Replaces the nadir point for \code{KSE}. If only a subset of values needs to be defined, 
##' the other coordinates can be set to \code{Inf}.
##' @param ... additional parameter to be passed to \code{fun}
##' @return
##' A list with components:
##' \itemize{
##' \item{\code{model}}{: a list of objects of class \code{\link[DiceKriging]{km}} corresponding to the last kriging models fitted.}
##' \item{\code{Jplus}}{: recorded values of the acquisition function maximizer}
##' \item{\code{integ.pts} and  \code{expanded.indices}}{: the discrete space used,}
##' \item{\code{predEq}}{: a list containing the recorded values of the estimated best solution,}
##' \item{\code{Eq.design, Eq.poff}}{: estimated equilibrium and corresponding pay-off}
##' }
##'
##' Note: with CKSE, kweights are not used when the mean on integ.pts is used
##'
##' @export
## ' @importFrom grDevices dev.off pdf rainbow
## ' @importFrom graphics axis pairs par points title
##' @importFrom stats pnorm qnorm rnorm dnorm runif
##' @importFrom methods slot
## ' @importFrom graphics filled.contour
##' @import DiceKriging DiceDesign parallel
##' @importFrom MASS mvrnorm
## ' @importFrom grDevices terrain.colors
## ' @importFrom graphics legend
##' @useDynLib GPGame, .registration = TRUE
##' @references
##' V. Picheny, M. Binois, A. Habbal (2016+), A Bayesian optimization approach to find Nash equilibria,
##' \emph{https://arxiv.org/abs/1611.02440}.
##' @examples
##' \dontrun{
##'
##' ################################################################
##' # Example 1: Nash equilibrium, 2 variables, 2 players, no filter
##' ################################################################
##' # Define objective function (R^2 -> R^2)
##' fun1 <- function (x)
##' {
##'   if (is.null(dim(x)))    x <- matrix(x, nrow = 1)
##'   b1 <- 15 * x[, 1] - 5
##'   b2 <- 15 * x[, 2]
##'   return(cbind((b2 - 5.1*(b1/(2*pi))^2 + 5/pi*b1 - 6)^2 + 10*((1 - 1/(8*pi)) * cos(b1) + 1),
##'                -sqrt((10.5 - b1)*(b1 + 5.5)*(b2 + 0.5)) - 1/30*(b2 - 5.1*(b1/(2*pi))^2 - 6)^2-
##'                 1/3 * ((1 - 1/(8 * pi)) * cos(b1) + 1)))
##' }
##'
##' # To use parallel computation (turn off on Windows)
##' library(parallel)
##' parallel <- FALSE #TRUE #
##' if(parallel) ncores <- detectCores() else ncores <- 1
##'
##' # Simple configuration: no filter, discretization is a 21x21 grid
##'
##' # Grid definition
##' n.s <- rep(21, 2)
##' x.to.obj   <- c(1,2)
##' gridtype <- 'cartesian'
##'
##' # Run solver with 6 initial points, 4 iterations
##' # Increase n.ite to at least 10 for better results
##' res <- solve_game(fun1, equilibrium = "NE", crit = "sur", n.init=6, n.ite=4,
##'                   d = 2, nobj=2, x.to.obj = x.to.obj,
##'                   integcontrol=list(n.s=n.s, gridtype=gridtype),
##'                   ncores = ncores, trace=1, seed=1)
##'
##' # Get estimated equilibrium and corresponding pay-off
##' NE <- res$Eq.design
##' Poff <- res$Eq.poff
##'
##' # Draw results
##' plotGame(res)
##'
##'################################################################
##' # Example 2: same example, KS equilibrium with given Nadir
##'################################################################
##' # Run solver with 6 initial points, 4 iterations
##' # Increase n.ite to at least 10 for better results
##' res <- solve_game(fun1, equilibrium = "KSE", crit = "sur", n.init=6, n.ite=4,
##'                   d = 2, nobj=2, x.to.obj = x.to.obj,
##'                   integcontrol=list(n.s=400, gridtype="lhs"),
##'                   ncores = ncores, trace=1, seed=1, Nadir=c(Inf, -20))
##'
##' # Get estimated equilibrium and corresponding pay-off
##' NE <- res$Eq.design
##' Poff <- res$Eq.poff
##'
##' # Draw results
##' plotGame(res, equilibrium = "KSE", Nadir=c(Inf, -20))
##'
##' ################################################################
##' # Example 3: Nash equilibrium, 4 variables, 2 players, filtering
##' ################################################################
##' fun2 <- function(x, nobj = 2){
##'   if (is.null(dim(x)))     x <- matrix(x, 1)
##'   y <- matrix(x[, 1:(nobj - 1)], nrow(x))
##'   z <- matrix(x[, nobj:ncol(x)], nrow(x))
##'   g <- rowSums((z - 0.5)^2)
##'   tmp <- t(apply(cos(y * pi/2), 1, cumprod))
##'   tmp <- cbind(t(apply(tmp, 1, rev)), 1)
##'   tmp2 <- cbind(1, t(apply(sin(y * pi/2), 1, rev)))
##'   return(tmp * tmp2 * (1 + g))
##' }
##'
##' # Grid definition: player 1 plays x1 and x2, player 2 x3 and x4
##' # The grid is a lattice made of two LHS designs of different sizes
##' n.s <- c(44, 43)
##' x.to.obj   <- c(1,1,2,2)
##' gridtype <- 'lhs'
##'
##' # Set filtercontrol: window filter applied for integration and candidate points
##' # 500 simulation and 200 candidate points are retained.
##' filtercontrol <- list(nsimPoints=500, ncandPoints=200,
##'                    filter=c("window", "window"))
##'
##' # Set km control: lower bound is specified for the covariance range
##' # Covariance type and model trend are specified
##' kmcontrol <- list(lb=rep(.2,4), model.trend=~1, covtype="matern3_2")
##'
##' # Run solver with 20 initial points, 4 iterations
##' # Increase n.ite to at least 20 for better results
##' res <- solve_game(fun2, equilibrium = "NE", crit = "psim", n.init=20, n.ite=2,
##'                   d = 4, nobj=2, x.to.obj = x.to.obj,
##'                   integcontrol=list(n.s=n.s, gridtype=gridtype),
##'                   filtercontrol=filtercontrol,
##'                   kmcontrol=kmcontrol,
##'                   ncores = 1, trace=1, seed=1)
##'
##' # Get estimated equilibrium and corresponding pay-off
##' NE <- res$Eq.design
##' Poff <- res$Eq.poff
##'
##' # Draw results
##' plotGame(res)
##'
##' ################################################################
##' # Example 4: same example, KS equilibrium
##' ################################################################
##'
##' # Grid definition: simple lhs
##' integcontrol=list(n.s=1e4, gridtype='lhs')
##'
##' # Run solver with 20 initial points, 4 iterations
##' # Increase n.ite to at least 20 for better results
##' res <- solve_game(fun2, equilibrium = "KSE", crit = "sur", n.init=20, n.ite=2,
##'                   d = 4, nobj=2,
##'                   integcontrol=integcontrol,
##'                   filtercontrol=filtercontrol,
##'                   kmcontrol=kmcontrol,
##'                   ncores = 1, trace=1, seed=1)
##'
##' # Get estimated equilibrium and corresponding pay-off
##' NE <- res$Eq.design
##' Poff <- res$Eq.poff
##'
##' # Draw results
##' plotGame(res, equilibrium = "KSE")
##' }

solve_game <- function(
  fun, ..., equilibrium="NE", crit="sur", model=NULL, n.init=NULL, n.ite, d, nobj, x.to.obj=NULL, noise.var = NULL,
  Nadir=NULL, integcontrol=NULL, simucontrol=NULL, filtercontrol=NULL, kmcontrol=NULL, returncontrol=NULL,
  ncores=1, trace=1, seed=NULL) {
  
  t1 <- Sys.time()
  set.seed(seed)
  if(ncores > 1) parallel <- TRUE
  
  ####################################################################################################
  #### CHECK INPUTS ##################################################################################
  ####################################################################################################
  # Check critcontrols
  if (crit == 'pex') PnashMethod <- 'exact' else PnashMethod <- 'simu'
  if (crit %in% c('psim', 'pex') && equilibrium!="NE"){
    cat("pex and psim available only for Nash equilibria; crit switched to sur \n")
    crit <- "sur"
  }
  
  # Check integcontrol
  integ.pts <- integcontrol$integ.pts
  expanded.indices <- integcontrol$expanded.indices
  n.s <- integcontrol$n.s
  gridtype <- switch(1+is.null(integcontrol$gridtype), integcontrol$gridtype, "lhs")
  lb <- switch(1+is.null(integcontrol$lb), integcontrol$lb, rep(0, d))
  ub <- switch(1+is.null(integcontrol$ub), integcontrol$ub, rep(1, d))
  integcontrol$renew <- switch(1+is.null(integcontrol$renew), integcontrol$renew, equilibrium %in% c("KSE", "CKSE"))
  integcontrol$kweights <- switch(1+is.null(integcontrol$kweights), integcontrol$kweights, FALSE)
  integcontrol$nsamp <- switch(1+is.null(integcontrol$nsamp), integcontrol$nsamp, 1e4)
  if (equilibrium %in% c("NE", "NKSE") && integcontrol$renew) {
    cat("integcontrol$renew=TRUE only available for KSE and CKSE \n")
    integcontrol$renew <- FALSE
  }
  
  # Check simucontrol
  n.ynew <- switch(1+is.null(simucontrol$n.ynew), simucontrol$n.ynew, 10)
  n.sim <- switch(1+is.null(simucontrol$n.sim), simucontrol$n.sim, 10)
  IS <- switch(1+is.null(simucontrol$IS), simucontrol$IS, TRUE)
  cross <- FALSE
  
  # Check filtercontrol
  filter <- switch(1+is.null(filtercontrol$filter), filtercontrol$filter, 'none')
  if (length(filter)==1) filter <- rep(filter, 2)
  randomFilter <- switch(1+is.null(filtercontrol$randomFilter), filtercontrol$randomFilter, FALSE)
  if (length(randomFilter)==1) randomFilter <- rep(randomFilter, 2)
  nsimPoints <- switch(1+is.null(filtercontrol$nsimPoints), filtercontrol$nsimPoints, 800)
  ncandPoints <- switch(1+is.null(filtercontrol$ncandPoints), filtercontrol$ncandPoints, 200)
  filtercontrol$nsamp <- switch(1+is.null(filtercontrol$nsamp), filtercontrol$nsamp, max(20, 5* nobj))
  
  if (equilibrium %in% c("KSE", "CKSE") && "Pnash" %in% filter) {
    cat("Pnash filter only available for NE; switching to PND \n")
    filter[which(filter=="Pnash")] <- 'PND'
  }
  simu.filter <- filter[1]
  cand.filter <- filter[2]
  
  # Check kmcontrol
  cov.reestim <- switch(1+is.null(kmcontrol$cov.reestim), kmcontrol$cov.reestim, TRUE) 
  model.trend <- switch(1+is.null(kmcontrol$model.trend), kmcontrol$model.trend, ~1)
  kmlb <- switch(1+is.null(kmcontrol$lb), kmcontrol$lb, rep(.1,d))
  kmub <- switch(1+is.null(kmcontrol$ub), kmcontrol$ub, rep(1,d))
  kmnugget <- switch(1+is.null(kmcontrol$nugget), kmcontrol$nugget, 1e-8)
  control <- switch(1+is.null(kmcontrol$control), kmcontrol$control, list(trace=FALSE))
  covtype <- switch(1+is.null(kmcontrol$covtype), kmcontrol$covtype, "matern5_2")
  
  # Check returncontrol and track.Eq
  if (!is.null(returncontrol$track.Eq)) track.Eq <- returncontrol$track.Eq else track.Eq <- 'none'
  if(is.null(returncontrol$force.exploit.last)) returncontrol$force.exploit.last <- TRUE
  if (equilibrium %in% c("KSE", "CKSE") && track.Eq %in% c("psim", "pex")) {
    cat("track.Eq must be set to none, mean or empirical for KSE - switching to empirical \n")
    track.Eq <- "empirical"
  } 
  if (equilibrium %in% c("NE", "NKSE") && track.Eq == "empirical") {
    # Pnash case
    cat("track.Eq= empirical option only available for KSE and CKSE; switching to psim \n")
    track.Eq <- "psim"
  }
  
  # Check Noise
  if (!is.null(noise.var)) {
    if (typeof(noise.var) == "closure") noise.var <- match.fun(noise.var)
    else if (typeof(noise.var) == "double" && length(noise.var==1)) noise.var <- rep(noise.var, nobj)
    kmnugget <- NULL
  }
  
  if (trace>0) cat("--------------------------\n Starting", equilibrium, "search with:", "\n",
                   crit, "strategy,","\n",
                   simu.filter, "simulation point filter,", "\n",
                   cand.filter, "candidate point filter,", "\n",
                   "among (", n.s, ") strategies \n --------------------------\n " )
  
  ####################################################################################################
  #### INITIALIZE VARIABLES AND MODELS ###############################################################
  ####################################################################################################
  
  ### Initialize variables
  window <- NULL
  all.Jplus <- c()
  Eq.design.estimate <- Eq.poff.estimate <- trackPnash <- NULL
  predEq <- list()
  include.obs <- FALSE
  
  #### Initial design and models ####################
  if (is.null(model)){
    if (is.null(n.init))     n.init <- 5*d
    design <- lhsDesign(n.init, d, seed = seed)$design
    response <- t(apply(design, 1, fun, ... = ...))
    
    if (!is.null(noise.var)) {
      if (typeof(noise.var) == "closure") {
        newnoise.var <- apply(design, 1, noise.var, ...)
      } else if (typeof(noise.var) == "double") {
        newnoise.var <- matrix(noise.var, nrow=n.init, ncol=ncol(response), byrow=TRUE)
      } else {#noise.var ="given_by_fn"
        # newnoise.var <- response[[2]]
        # ynew <- response[[1]]
        tmp <- newnoise.var <- NULL # initialization
        for (i in 1:length(response)) {
          tmp <- rbind(tmp, response[[i]][[1]])
          newnoise.var <- rbind(newnoise.var, response[[i]][[2]])
        }
        response <- tmp
      }
    } else {
      newnoise.var <- NULL
    }
    nobj <- ncol(response)
    
    my.km <- function(i) {
      km(model.trend, design = design, response = response[, i], covtype=covtype,
         control=control, lower=kmlb, upper=kmub, nugget=kmnugget, noise.var=newnoise.var[,i])
    }
    model <- mclapply(1:nobj, my.km, mc.cores=ncores)
  }
  
  #### Integration points ###########################
  if (is.null(integcontrol$integ.pts) || is.null(integcontrol$expanded.indices)) {
    if (equilibrium %in% c("KSE", "CKSE")) include.obs <- TRUE else include.obs <- FALSE
    res <- generate_integ_pts(n.s=n.s, d=d, nobj=nobj, x.to.obj=x.to.obj, equilibrium=equilibrium,
                              gridtype=gridtype, lb = lb, ub = ub, include.obs=include.obs, model=model)
    integcontrol$integ.pts <- integ.pts <- res$integ.pts
    integcontrol$expanded.indices <- expanded.indices <- res$expanded.indices
  }
  n.integ.pts <- nrow(integ.pts)
  cand.pts <- integ.pts
  
  if (is.null(expanded.indices)) {
    sorted <- FALSE
  } else {
    sorted <- !is.unsorted(expanded.indices[, ncol(expanded.indices)])
  }
  t2 <- Sys.time() - t1
  if (trace>2) cat("Time for initialization: ", t2, units(t2), "\n")
  t1 <- Sys.time()
  
  ####################################################################################################
  #### MAIN LOOP STARTS HERE #########################################################################
  ####################################################################################################
  Jnres <- rep(NA, n.ite + 1)
  
  for (ii in 1:(n.ite+1)){
    if (ii < (n.ite+1) && trace>0){cat("--- Iteration #", ii, " ---\n")}else{
      if(trace>0){cat}("--- Post-processing ---\n")
    }
    t0 <- t1 <- Sys.time()
    
    
    ##################################################################
    # RENEW INTEGRATION POINTS
    ##################################################################
    if (ii > 1 && integcontrol$renew) {
      include.obs <- equilibrium %in% c("KSE", "CKSE")
      res <- generate_integ_pts(n.s=n.s, d=d, nobj=nobj, x.to.obj=x.to.obj, equilibrium=equilibrium,
                                gridtype=gridtype, lb = lb, ub = ub, include.obs=include.obs, model=model, seed=ii)
      integcontrol$integ.pts <- integ.pts <- res$integ.pts
      integcontrol$expanded.indices <- expanded.indices <- res$expanded.indices
    }
    ##################################################################
    # MODELS UPDATE
    ##################################################################
    if (ii > 1) {
      
      # if (ii == n.ite){
      #   crit <- finalcrit
      # }
      if (ii == (n.ite+1) && !(equilibrium %in% c("KSE", "CKSE"))){
        cand.filter <- 'none'
      }
      
      newmodel <- model
      X.new <- matrix(xnew,nrow=1)
      
      my.update <- function(u) {
        try(update(object = model[[u]], newX = X.new, newy=ynew[u], newX.alreadyExist=FALSE, newnoise.var = newnoise.var[u],
                   cov.reestim = cov.reestim, kmcontrol = list(control = list(trace = FALSE))), silent = TRUE)
      }
      newmodel <- mclapply(1:nobj, my.update, mc.cores=ncores)
      
      for (u in 1:nobj){
        if (typeof(newmodel[[u]]) == "character" && cov.reestim) {
          cat("Error in hyperparameter estimation - old hyperparameter values used instead for model ", u, "\n")
          newmodel[[u]] <- try(update(object = model[[u]], newX = X.new, newy=ynew[u], newnoise.var = newnoise.var[u],
                                      newX.alreadyExist=FALSE, cov.reestim = FALSE), silent = TRUE)
        }
        if (typeof(newmodel[[u]]) == "character") {
          cat("Unable to udpate kriging model ", u, " at iteration", ii-1, "- optimization stopped \n")
          cat("lastmodel ", u, " is the model at iteration", ii-1, "\n")
          cat("par and values contain the ",ii , "th observation \n \n")
          # if (ii > 1) allX.new <- rbind(model[[u]]@X[(model[[u]]@n + 1):(model[[u]]@n+ii-1),, drop=FALSE], X.new)
          return(list(
            par    = X.new,
            values = ynew,
            nsteps = ii,
            model = model,
            Jplus = Jplus, integcontrol=integcontrol))
        } else {
          model[[u]] <- newmodel[[u]]
        }
      }
    }
    
    t2 <- Sys.time() - t1
    if (trace>2) cat("Time for models updates: ", t2, units(t2), "\n")
    t1 <- Sys.time()
    observations <- Reduce(cbind, lapply(model, slot, "y"))
    
    #--------------------------------------------------------#
    # Regular loop : find infill point
    #--------------------------------------------------------#
    # if (ii < (n.ite+1) || return.Eq){
    # if (ii==(n.ite+1)) include.obs <- TRUE
    
    pred <- mclapply(model, FUN=predict, newdata = integ.pts, checkNames = FALSE, type = "UK", light.return = TRUE, mc.cores=ncores)
    
    t2 <- Sys.time() - t1
    if (trace>2) cat("Time for predictions: ", t2, units(t2), "\n")
    t1 <- Sys.time()
    
    # if (equilibrium == "KSE" && is.null(noise.var)){
    #   # PFobs <- observations[nonDomInd(observations),, drop = FALSE]
    #   PFobs <- nonDom(observations)
    #   NSobs <- PFobs[unique(c(apply(PFobs, 2, which.min), # Shadow
    #                           apply(PFobs, 2, which.max) # Nadir
    #   )),, drop = FALSE]
    # } else {
    #   NSobs <- NULL
    # }
    
    # if (!is.null(nadir)) 
    
    ##################################################################
    # FILTERING INTEGRATION POINTS
    ##################################################################
    if (simu.filter == "window") {
      
      if (is.null(window) || crit != "sur") {
        # Special case for initialization - otherwise do not change the old "window" value
        predmean <- Reduce(cbind, lapply(pred, function(alist) alist$mean))
        Eq_simu <- getEquilibrium(Z = predmean,  equilibrium = equilibrium, nobj = nobj, expanded.indices=expanded.indices,
                                  n.s=n.s, sorted = sorted, cross = cross, kweights = NULL, Nadir=Nadir)#, NSobs = NSobs)
        if(nrow(Eq_simu) < 1 || all(is.na(Eq_simu))){
          ## if no window has been defined in previous iterations, use the mean, otherwise keep the old one
          if(is.null(window)){
            window <- colMeans(observations)
          }
        } else if(nrow(Eq_simu) == 1) {
          window <- as.vector(Eq_simu)
        } else {
          window <- colMeans(Eq_simu, na.rm = TRUE)
        }
      }
      options <- list(window = window)
    } else if (simu.filter == "Pnash") {
      options <- list(method=PnashMethod, nsim=100)
    }
    if (trace > 1 && simu.filter == "window") cat("window for integ.pts", window, "\n")
    
    # Reduce the number of integration points first
    if (simu.filter != 'none') {
      # include.obs <- equilibrium %in% c("KSE", "CKSE") || ii==(n.ite+1)
      
      filt <- filter_for_Game(n.s.target = pmin(n.s, round(nsimPoints^(1/length(n.s)))), integcontrol=integcontrol,
                              model = model, predictions=pred, type=simu.filter, options = options, random=randomFilter[1],
                              include.obs=TRUE, ncores = ncores, equilibrium=equilibrium, nsamp=filtercontrol$nsamp, Nadir=Nadir)
      I <- filt$I
    } else {
      I <- 1:nrow(integ.pts)
    }
    
    # ## Add potential minimas of each objective based on EI, and EI x Pdom
    # if (equilibrium == "KSE" && simu.filter != 'none'){
    #   IKS <- NULL # points of interest for KS (i.e., related to KS and Nadir)
    #   
    #   ## EI on each objective
    #   for (jj in 1:nobj) {
    #     # Based on DiceOptim
    #     xcr <-  (min(model[[jj]]@y) - pred[[jj]]$mean)/pred[[jj]]$sd
    #     test <- (min(model[[jj]]@y) - pred[[jj]]$mean)*pnorm(xcr) + pred[[jj]]$sd * dnorm(xcr)
    #     test[which(pred[[jj]]$sd/sqrt(model[[jj]]@covariance@sd2) < 1e-06)] <- NA
    #     IKS <- c(IKS, which.max(test))
    #   }
    #   
    #   ## EI x Pdom on each objective (to find potential Nadir points)
    #   PNDom <- prob.of.non.domination(paretoFront = PFobs,
    #                                   model = model, integration.points = integ.pts, predictions = pred,
    #                                   nsamp = filtercontrol$nsamp)
    #   for (jj in 1:nobj) {
    #     # Based on DiceOptim
    #     xcr <-  -(max(PFobs[,jj]) - pred[[jj]]$mean)/pred[[jj]]$sd
    #     test <- (-max(PFobs[,jj]) + pred[[jj]]$mean)*pnorm(xcr) + pred[[jj]]$sd * dnorm(xcr)
    #     test[which(pred[[jj]]$sd/sqrt(model[[jj]]@covariance@sd2) < 1e-06)] <- NA
    #     test <- test * PNDom
    #     IKS <- c(IKS, which.max(test))
    #   }
    #   
    #   I <- unique(c(I, IKS))
    # } else {
    #   IKS <- NULL
    # }
    
    my.integ.pts <- integ.pts[I,]
    my.expanded.indices <- expanded.indices[I,]
    my.pred <- lapply(pred, function(alist) list(mean=alist$mean[I], sd=alist$sd[I]))
    
    if (equilibrium %in% c("KSE", "CKSE") && length(n.s)==1) {
      my.n.s <- length(I)
    } else {
      my.n.s <- apply(my.expanded.indices, 2, function(x) length(unique(x)))
    }
    
    if (equilibrium == "CKSE" && integcontrol$kweights){
      Xsamp <- matrix(runif(d * integcontrol$nsamp), integcontrol$nsamp)
      kweights <- list()
      for(u in 1:nobj){
        kn <- covMat1Mat2(model[[u]]@covariance, Xsamp, my.integ.pts)
        Knn <- try(chol2inv(chol(covMat1Mat2(model[[u]]@covariance, my.integ.pts, my.integ.pts) + diag(1e-6, nrow(my.integ.pts)))))
        if (typeof(Knn)== "character") {
          cat("Unable to compute kweights, last model returned \n")
          return(list(model=model, predEq=predEq, integcontrol=integcontrol))
        }
        kweights <- c(kweights, list(list(kn = kn, Knn = Knn)))
      }
    } else {kweights = NULL}
    
    if (trace>1) cat("Number of strategies after simu filtering: (", my.n.s, ")\n")
    
    my.expanded.indices2 <- my.expanded.indices
    for (i in 1:ncol(my.expanded.indices)) {
      unik_i <- unique(my.expanded.indices[,i])
      for (j in 1:length(unik_i)) {
        my.expanded.indices2[which(my.expanded.indices[,i]==unik_i[j]),i] <- j
      }
    }
    my.expanded.indices <- my.expanded.indices2
    
    my.integcontrol <- list(expanded.indices=my.expanded.indices, integ.pts=my.integ.pts, n.s=my.n.s)
    
    t2 <- Sys.time() - t1
    if (trace>2) cat("Time for filter 1 (simu): ", t2, units(t2), "\n")
    t1 <- Sys.time()
    
    ##################################################################
    # Precalculations and simulations (if not performed already)
    ##################################################################
    precalc.data <- mclapply(model, FUN=precomputeUpdateData, integration.points=my.integ.pts, mc.cores=ncores)
    K <- my.pred[[1]]$sd/model[[1]]@covariance@sd2 < 1e-12
    index.obs <- which(K)
    index.other <- which(!K)
    
    if (crit == 'sur' || cand.filter=="window"){
      my.Simu <- matrix(NA, nrow=nrow(my.integ.pts), ncol=n.sim*nobj)
      current.simu <- try(t(Reduce(rbind, mclapply(model, simulate, nsim=n.sim, newdata=my.integ.pts[index.other,,drop=FALSE], 
                                                   cond=TRUE, checkNames=FALSE, nugget.sim = 10^-6, mc.cores=ncores))))
      
      if (typeof(current.simu)== "character") {
        cat("Unable to compute GP draws, last model returned \n")
        return(list(model=model, predEq=predEq, integcontrol=integcontrol))
      } else {
        my.Simu[index.other,] <- current.simu
      }
      
      if (length(index.obs)>0){
        my.obs <- t(Reduce(rbind, lapply(my.pred, function(alist) alist$mean[index.obs])))
        my.Simu[index.obs,] <- matrix(do.call("rbind", rep(list( my.obs), n.sim)), nrow=nrow( my.obs))
      }
    }
    
    t2 <- Sys.time() - t1
    if (trace>2) cat("Time for precalc + simu: ", t2, units(t2), "\n")
    t1 <- Sys.time()
    
    ##################################################################
    # FILTERING CANDIDATE POINTS
    ##################################################################
    cand.pts <- my.integ.pts
    if ((simu.filter == "window" && crit == 'sur') || cand.filter == "window") {
      Eq_simu <- getEquilibrium(my.Simu, equilibrium = equilibrium, nobj = nobj, expanded.indices=my.expanded.indices,
                                n.s=my.n.s, sorted = sorted, cross = cross, kweights = kweights, Nadir=Nadir)
      # ,NSobs = matrix(do.call("rbind", rep(list(NSobs), n.sim)), nrow=nrow(NSobs)))
      Eq_simu <- Eq_simu[which(!is.na(Eq_simu[,1])),, drop = FALSE]
      Jnres[ii] <- determinant(cov(Eq_simu), logarithm = TRUE)[[1]][1]
      temp <- apply(Eq_simu, 2, range, na.rm = TRUE)
      
      if (!any(is.na(temp))) window <- temp
      options <- list(window = window)
    } else if (cand.filter == "Pnash") {
      options <- list(method=PnashMethod, nsim = 100)
    }
    
    if (trace>1 && cand.filter == "window") cat("window for cand.pts", window, "\n")
    
    J <- 1:nrow(cand.pts) # No filter case
    if (cand.filter != 'none') {
      if (cand.filter == "Pnash" && simu.filter == "Pnash") {
        # Special case if Pnash used twice: do not redo calculations
        J <- I[1:min(length(I), ncandPoints)] # deterministic filter
        if (randomFilter[2]) {
          J <- sample.int(length(I), size=min(length(I), ncandPoints), replace=FALSE, prob=filt$crit[I])
        }
      } else {
        # Regular case
        # include.obs <- ii==(n.ite+1)
        J <- filter_for_Game(n.s.target = pmin(my.n.s, round(ncandPoints^(1/length(my.n.s)))), integcontrol=my.integcontrol,
                             model = model, predictions=my.pred, type=cand.filter, options = options,
                             random=randomFilter[2], include.obs=FALSE, ncores = ncores, 
                             equilibrium=equilibrium, nsamp=filtercontrol$nsamp, Nadir=Nadir)$I
      }
    }
    
    # ## Add extremal points for KSE / CKSE
    # if(!is.null(IKS)){
    #   J <- unique(c(J, which(I %in% IKS)))
    # }
    
    ## Remove already evaluated designs from candidates
    if((ii > 1) && (ii<(n.ite+1))){
      my.pred.s <- my.pred[[1]]$sd[J]
      J <- J[which(my.pred.s >= sqrt(model[[1]]@covariance@sd2)/1e6)]
    }
    my.cand.pts <- cand.pts[J,, drop=FALSE]
    cand.s <- my.expanded.indices[J,, drop=FALSE]
    
    t2 <- Sys.time() - t1
    if (trace>2) cat("Time for filter 2 (cand): ", t2, units(t2), "\n")
    t1 <- Sys.time()
    
    if (trace>0) cat("Nb of integration / candidate pts:", length(I), "/", length(J), "\n")
    
    ##################################################################
    # CRITERION MINIMIZATION
    ##################################################################
    my.integ.pts <- data.frame(my.integ.pts)
    if (ii < n.ite+1) {
      
      ## Special treatment to force exploration on last iteration
      FailToExploit <- TRUE
      if (ii == n.ite && returncontrol$force.exploit.last) {
        if (equilibrium %in% c("KSE", "CKSE")) {
          predmean <- Reduce(cbind, lapply(pred, function(alist) alist$mean))
          currentEq <- try(getEquilibrium(Z = predmean,  equilibrium = equilibrium, nobj = nobj,
                                          expanded.indices=expanded.indices, n.s=n.s, kweights = NULL,
                                          sorted = sorted, cross = cross, return.design = TRUE, Nadir=Nadir)) #, NSobs = NSobs)
          
          if (typeof(currentEq)== "character") {
            cat("Unable to compute exploitation step, last model returned \n")
            return(list(model=model, predEq=predEq, integcontrol=integcontrol))
          }
          
          i <- currentEq$NE
          xnew <- integ.pts[i,]
          FailToExploit <- FALSE
          if (duplicated(rbind(integ.pts[currentEq$NE,], model[[1]]@X), fromLast = TRUE)[1]) FailToExploit <- TRUE
        } else {
          crit <- ifelse(crit=="sur", "pex", crit)
        }
      }
      
      ## Regular case
      if (FailToExploit) {
        if (crit == "sur") {
          Jplus <- try(unlist(mclapply(X=J, FUN=crit_SUR_Eq, model=model, integcontrol=my.integcontrol,
                                       Simu=my.Simu, equilibrium = equilibrium, precalc.data=precalc.data,
                                       n.ynew=n.ynew, IS=IS, cross=cross, kweights = kweights,
                                       mc.cores = ncores, Nadir=Nadir)))
        } else if (crit == 'pex') {
          Jplus <- try(-crit_PNash(idx=J, integcontrol=my.integcontrol,
                                   type = 'exact', model = model, ncores = ncores))
        } else if (crit == "psim") {
          Jplus <- try(-crit_PNash(idx=J, integcontrol=my.integcontrol,
                                   type = 'simu', model = model, ncores = ncores, control = list(nsim = 100)))
        }
        if (typeof(Jplus)== "character") {
          cat("Unable to optimize criterion, last model returned \n")
          return(list(model=model, predEq=predEq, integcontrol=integcontrol))
        }
        if (all(is.na(Jplus))){
          cat("No equilibrium found - last model returned \n")
          return(list(model = model, predEq=predEq, Jplus = Jplus, integcontrol=integcontrol))
        }
        Jplus[is.na(Jplus)] <- max(Jplus, na.rm = TRUE)
        
        ## Get best solution
        i <- which.min(Jplus)
        xnew <- my.cand.pts[i,]
      }
      
      t2 <- Sys.time() - t1
      if (trace>2) cat("Time for optim: ", t2, units(t2), "\n")
      
      if (trace > 0 && crit=="sur") cat("Jplus:", min(Jplus), '\n')
      if (trace > 0 && crit!="sur") cat("Pmax:", -min(Jplus), '\n')
      
      ## Compute new observation
      ynew <- try(fun(xnew, ...))
      
      if (typeof(ynew) == "character" ) {
        cat("Unable to compute objective function at iteration ", i, "- optimization stopped \n")
        cat("Problem occured for the design: ", xnew, "\n")
        cat("Last model and problematic design (xnew) returned \n")
        return(list(model=model, predEq=predEq, integcontrol=integcontrol, xnew=xnew))
      }
      
      if (!is.null(noise.var)) {
        if (typeof(noise.var) == "closure") {
          newnoise.var <- noise.var(xnew)
        } else if (typeof(noise.var) == "double") {
          newnoise.var <- noise.var
        } else {#noise.var ="given_by_fn"
          newnoise.var <- ynew[[2]]
          ynew <- ynew[[1]]
        }
      } else {
        newnoise.var <- NULL
      }
      
      all.Jplus <- c(all.Jplus, min(Jplus))
      
      if (trace>1) cat("New design evaluated: \n")
      if (trace>1) cat(xnew, "\n")
      if (trace>1) cat("Corresponding new observation: \n")
      if (trace>1) cat(ynew, "\n")
      
      t2 <- Sys.time() - t0
      if (trace>0) cat("Total iteration time: ", t2, units(t2), "\n")
      t1 <- Sys.time()
    }
    ##################################################################
    # Tracking of estimated best solution
    ##################################################################
    
    if (track.Eq != 'none' || (ii == n.ite+1)) {
      
      if (track.Eq == "empirical") {
        observations <- Reduce(cbind, lapply(model, slot, "y"))
        currentEq <- try(getEquilibrium(Z=observations, equilibrium = equilibrium, nobj=nobj,
                                        return.design=TRUE, cross=cross, kweights = kweights, Nadir=Nadir))
        
      } else if (track.Eq %in% c("psim", "pex")) {
        
        if (simu.filter == "Pnash") {
          estPnash <- filt$crit
        } else {
          if (track.Eq == 'pex') mytype <- "exact"
          else                   mytype <- "simu"
          estPnash <- try(crit_PNash(idx = 1:nrow(my.integcontrol$expanded.indices), integcontrol=my.integcontrol,
                                     type = mytype, model = model, ncores = ncores,
                                     control = list(nsim = 100)))
        }
        if (typeof(ynew) == "character" ) {
          currentEq <- NULL
        } else {
          ImaxPnash <- which.max(estPnash)
          maxPnash <- max(estPnash)
          currentEq <- list(predEqPoff = unlist(lapply(my.pred, function(x) x$mean[ImaxPnash])),
                            Eq = my.integ.pts[ImaxPnash,, drop = F], Pnash = maxPnash)
          if (trace>2) cat("Maximum estimated Pnash:", maxPnash," \n")
        }
      } else {
        ## Eq of the GP predictive means
        predmean <- Reduce(cbind, lapply(pred, function(alist) alist$mean))
        currentEq <- try(getEquilibrium(Z = predmean,  equilibrium = equilibrium, nobj = nobj,
                                        expanded.indices=expanded.indices, n.s=n.s, kweights = NULL,
                                        sorted = sorted, cross = cross, return.design = T, Nadir=Nadir)) #, NSobs = NSobs)
      }
      predEq <- c(predEq, list(currentEq))
      
      t2 <- Sys.time() - t1
      if (trace>0) cat("Time for estimating best candidate: ", t2, units(t2), "\n")
      t1 <- Sys.time()
    }
  }
  ####################################################################################################
  #### MAIN LOOP ENDS HERE ###########################################################################
  ####################################################################################################
  
  #--------------------------------------------------------#
  # When optimization is over: find best Eq estimate
  #--------------------------------------------------------#
  Eq.design.estimate <- integ.pts[predEq[[length(predEq)]]$NE,]
  Eq.poff.estimate <- predEq[[length(predEq)]]$NEPoff
  
  return(list(model = model, Eq.design=Eq.design.estimate, Eq.poff=Eq.poff.estimate,
              Jplus = all.Jplus, integcontrol=integcontrol, predEq = predEq))
  
}