#----------------------------------------------------------------
##' Select candidate points for conditional simulations or for criterion evaluation, based on a "window" or a probability related to the equilibrium at hand.
##' @title All-purpose filter
##' @param n.s.target scalar or vector of number of strategies (one value per player) to select. For \code{NE}, if \code{n.s.target} is a scalar
##' then each player will have \code{round(n.s.target^(1/nobj)} strategies.
##' @param integcontrol is a list containing: \code{integ.pts}, a [\code{npts x dim}] matrix defining the grid,
##' \code{expanded.indices} a matrix containing the indices of the integ.pts on the grid and \code{n.s},
##' a \code{nobj} vector containting the number of strategies per player
##' @param model is a list of \code{nobj} \code{nobj} \code{\link[DiceKriging]{km}} objects
##' @param predictions is a list of size \code{nobj}
##' @param type either "\code{window}", "\code{PND}" or "\code{Pnash}", see details
##' @param equilibrium either '\code{NE}', '\code{KSE}' or '\code{NKSE}' for Nash/Kalai-Smoridinsky/Nash-Kalai-Smoridinsky equilibria
##' @param options a list containing either the window (matrix or target) or the parameters for Pnash: method
##' ("\code{simu}" or "\code{exact}") and \code{nsim}
##' @param ncores \code{\link[parallel]{mclapply}} is used if \code{> 1} for parallel evaluation
##' @param random Boolean. If \code{FALSE}, the best points according to the filter criterion are chosen,
##' otherwise the points are chosen by random sampling with weights proportional to the criterion.
##' @param include.obs Boolean. If \code{TRUE}, the observations are included to the filtered set.
##' @param min.crit Minimal value for the criterion, useful if \code{random = TRUE}.
##' @param nsamp number of samples to estimate the probability of non-domination, useful when \code{type=PND} and \code{nobj}>3.
##' @param Nadir optional vector of size \code{nobj}. Replaces the nadir point for \code{KSE}. If only a subset of values needs to be defined, 
##' the other coordinates can be set to \code{Inf}.
##' @return List with two elements: \code{I} indices selected and \code{crit} the filter metric at all candidate points
##' @details If \code{type == "windows"}, points are ranked based on their distance to \code{option$window} (when it is a target vector),
##' or based on the probability that the response belongs to \code{option$window}.
##' The other options, "\code{PND}" (probability of non-domination, i.e., of not being dominated by the current Pareto front)
##' and "\code{Pnash}" (probability of realizing a Nash equilibrium) base the ranking of points on the associated probability.
##' @export
filter_for_Game <- function(n.s.target, model=NULL, predictions=NULL, type="window", equilibrium="NE",
                            integcontrol, options = NULL, ncores = 1, random=TRUE, include.obs=FALSE,
                            min.crit = 1e-12, nsamp = NULL, Nadir=NULL) {
  
  expanded.indices <- integcontrol$expanded.indices
  integ.pts <- integcontrol$integ.pts
  nobj <- ncol(expanded.indices)
  if (is.null(Nadir)) Nadir <- rep(Inf, nobj)
  
  if (is.null(options) && type=="Pnash") options <- list(method = 'simu', nsim = 100)
  if(type == "Pnash" && is.null(options$method)){
    options$method <- 'simu'
    options$nsim <- 100
  }
  if ((length(n.s.target) == 1) && (equilibrium %in% c("NE", "NKSE")))
    n.s.target <- rep(round(n.s.target^-nobj), nobj)
  
  if (equilibrium %in% c("KSE", "CKSE") && type=="Pnash") {
    cat("Pnash filter not available for KSE; switching to PND \n")
    type <- "PND"
  }
  
  if (type == "PND") {
    # Proba of non-domination
    crit <- prob.of.non.domination(model=model, integration.points=integ.pts, predictions=predictions, nsamp=nsamp)
    #---------------------------------
  } else if (type == "window") {
    crit <- rep(1, nrow(integ.pts))
    if (is.null(nrow(options$window))) {
      # Density at target
      for (u in 1:length(model)) {
        crit <- crit * dnorm( (options$window[u] - predictions[[u]]$mean)/predictions[[u]]$sd)
      }
    } else {
      # Proba of being in a box
      for (u in 1:length(model)) {
        crit <- crit * ( pnorm( (options$window[2,u] - predictions[[u]]$mean)/predictions[[u]]$sd) -
                           pnorm( (options$window[1,u] - predictions[[u]]$mean)/predictions[[u]]$sd) )
      }
    }
    #---------------------------------
  } else if (type == "Pnash") {
    crit <- crit_PNash(idx = 1:nrow(expanded.indices), integcontrol=integcontrol, model = model,
                       type = options$method, ncores = ncores, control = options)
  }
  crit2 <- crit[which(is.na(crit))] <- 0 #crit2 useful for KSE
  crit <- pmax(crit, min.crit)
  
  # More indices than needed in case of replications -> Fixed, selection is now with mat_s
  if (random) {
    if (nrow(integ.pts) < 1e5) {
      idx <- sample.int(nrow(integ.pts), replace = FALSE, prob = crit)
    } else {
      idx <- sample.int(nrow(integ.pts), size = 10*max(n.s.target), replace = FALSE, prob = crit)
    }
  } else {
    idx <- order(crit, decreasing = TRUE)
  }
  
  if (include.obs) {
    tmp <- duplicated(rbind(model[[1]]@X, integ.pts))
    J <- which(tmp[model[[1]]@n + (1:nrow(integ.pts))])
    idx <- unique(c(J, idx))
  }
  
  if (equilibrium %in% c("KSE", "CKSE") && type != "PND") {
    #---------------------------------
    ## Add potential minimas of each objective based on EI (utopia), and EI x Pdom (nadir)
    IKS <- NULL # points of interest for KS (i.e., related to KS and Nadir)
    # if (type == "PND") {
    #   PNDom <- crit2
    # } else {
    if (max(Nadir) == Inf) {
      PNDom <- prob.of.non.domination(model=model, integration.points=integ.pts, predictions=predictions, nsamp =nsamp)
    }
    # }
    PFobs <- nonDom(Reduce(cbind, lapply(model, slot, "y")))
    for (jj in 1:length(model)) {
      discard <- which(predictions[[jj]]$sd/sqrt(model[[jj]]@covariance@sd2) < 1e-06)
      # EI(min) on each objective (to find potential Utopia points)
      xcr <-  (min(model[[jj]]@y) - predictions[[jj]]$mean)/predictions[[jj]]$sd
      test <- (min(model[[jj]]@y) - predictions[[jj]]$mean)*pnorm(xcr) + predictions[[jj]]$sd * dnorm(xcr)
      test[discard] <- NA
      IKS <- c(IKS, which.max(test))
      # EI(max) x Pdom on each objective (to find potential Nadir points) unless Nadir is provided
      if (Nadir[jj] == Inf) {
        xcr <-  -(max(PFobs[,jj]) - predictions[[jj]]$mean)/predictions[[jj]]$sd
        test <- (-max(PFobs[,jj]) + predictions[[jj]]$mean)*pnorm(xcr) + predictions[[jj]]$sd * dnorm(xcr)
        test[discard] <- NA
        test <- test * PNDom
        IKS <- c(IKS, which.max(test))
      }
    }
    idx <- unique(c(IKS, idx))
    #---------------------------------
    I <- idx[1:n.s.target]
  } else {
    mat_s <- matrix(NA, max(n.s.target), length(n.s.target))
    for(i in 1:length(n.s.target)){
      mat_s[,i] <- unique(expanded.indices[idx, i])[1:n.s.target[i]]
    }
    
    I <- rep(TRUE, nrow(expanded.indices))
    for(i in 1:ncol(expanded.indices)){
      I <- I & (expanded.indices[,i] %in% mat_s[,i])
    }
    I <- which(I)
  }
  
  return(list(I=I, crit=crit))
}
