#' Calculate MSE over single dimension
#' 
#' Calculated using grid of integration points.
#' Can be calculated exactly, but not much reason in 1D.
#'
#' @param xl Vector of points in 1D
#' @param theta Correlation parameters
#' @param CorrMat Function that gives correlation matrix for vectors of 1D points.
#'
#' @return MSE value
#' @export
#'
#' @examples
#' CGGP_internal_calcMSE(xl=c(0,.5,.9), theta=c(1,2,3),
#'          CorrMat=CGGP_internal_CorrMatCauchySQT)
CGGP_internal_calcMSE <- function(xl, theta, CorrMat) {
  S = CorrMat(xl, xl, theta)
  xp = seq(-10^(-4),1+10^(-4),l=401)
  Cp = CorrMat(xp,xl,theta)
  n = length(xl)
  cholS = chol(S)
  CiCp = backsolve(cholS,backsolve(cholS,t(Cp), transpose = TRUE))
  MSE_MAPal = mean(1 - rowSums(t(CiCp)*Cp))
  
  MSE_MAPal
}


#' Calculate MSE over blocks
#' 
#' Delta of adding block is product over i=1..d of IMSE(i,j-1) - IMSE(i,j)
#'
#' @param valsinds Block levels to calculate MSEs for
#' @param MSE_MAP Matrix of MSE values
#'
#' @return All MSE values
#' @export
#'
#' @examples
#' SG <- CGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SG <- CGGPfit(SG, Y=y)
#' MSE_MAP <- outer(1:SG$d, 1:8, 
#'      Vectorize(function(dimlcv, lcv1) {
#'         CGGP_internal_calcMSE(SG$xb[1:SG$sizest[dimlcv]],
#'         theta=SG$thetaMAP[(dimlcv-1)*SG$numpara+1:SG$numpara],
#'         CorrMat=SG$CorrMat)
#'  }))
#' CGGP_internal_calcMSEde(SG$po[1:SG$poCOUNT, ], MSE_MAP)
CGGP_internal_calcMSEde <- function(valsinds, MSE_MAP) {
  maxparam <- -Inf # Was set to -10 and ruined it.
  if(is.matrix(valsinds)){
    MSE_de = rep(0, dim(valsinds)[1])
    
    for (levellcv2 in 1:dim(valsinds)[1]) {
      MSE_de[levellcv2] = 0
      for (levellcv in 1:dim(valsinds)[2]) {
        if (valsinds[levellcv2, levellcv] > 1.5) {
          MSE_de[levellcv2] = MSE_de[levellcv2] + max(log(-MSE_MAP[levellcv, valsinds[levellcv2, levellcv]] + 
                                                            MSE_MAP[levellcv, valsinds[levellcv2, levellcv] - 1]),maxparam)
          
        } else {
          # This is when no ancestor block, 1 comes from when there is no data. 
          # 1 is correlation times integrated value over range.
          # This depends on correlation function.
          MSE_de[levellcv2] = MSE_de[levellcv2] + max(log(-MSE_MAP[levellcv, valsinds[levellcv2, levellcv]] + 1),maxparam)
          
        }
      }
    }
  } else {
    MSE_de = 0
    
    for (levellcv in 1:length(valsinds)) {
      if (valsinds[levellcv] > 1.5) {
        MSE_de = MSE_de + max(log(-MSE_MAP[levellcv, valsinds[levellcv]] + MSE_MAP[levellcv, valsinds[levellcv] -1]),maxparam)
        
      } else {
        MSE_de = MSE_de + max(log(-MSE_MAP[levellcv, valsinds[levellcv]] + 1),maxparam)
        
      }
    }
  }
  
  MSE_de = exp(MSE_de)
  
  return(MSE_de)
}



#' Add points to CGGP
#' 
#' Add `batchsize` points to `SG` using `theta`.
#'
#' @param CGGP Sparse grid object
#' @param batchsize Number of points to add
#' @param selectionmethod How points will be selected: one of `UCB`, `TS`,
#' `MAP`, `Oldest`, `Random`, or `Lowest`.
#' `UCB` uses Upper Confidence Bound estimates for the parameters.
#' `TS` uses Thompson sampling, a random sample from the posterior.
#' `MAP` uses maximum a posteriori parameter estimates.
#' `Oldest` adds the block that has been available the longest.
#' `Random` adds a random block.
#' `Lowest` adds the block with the lowest sum of index levels.
#' `UCB` and `TS` are based on bandit algorithms and account for uncertainty
#' in the parameter estimates, but are the slowest.
#' `MAP` is fast but doesn't account for parameter uncertainty.
#' The other three are naive methods that are not adaptive and won't
#' perform well.
#' @importFrom stats quantile sd var
#'
#' @return SG with new points added.
#' @export
#' @family CGGP core functions
#'
#' @examples
#' SG <- CGGPcreate(d=3, batchsize=100)
#' y <- apply(SG$design, 1, function(x){x[1]+x[2]^2})
#' SG <- CGGPfit(SG, Y=y)
#' SG <- CGGPappend(CGGP=SG, batchsize=20, selectionmethod="MAP")
CGGPappend <- function(CGGP,batchsize, selectionmethod = "MAP"){
  # ===== Check inputs =====
  if (!(selectionmethod %in% c("UCB", "TS", "MAP", "Oldest", "Random", "Lowest"))) {
    stop("selectionmethod in CGGPappend must be one of UCB, TS, MAP, Oldest, Random, or Lowest")
  }
  
  if (!is.null(CGGP$design_unevaluated)) {
    stop("Can't append if CGGP has unevaluated design points.")
  }
  
  # Track how many design points there currently are in $design
  n_before <- if (is.null(CGGP[["design"]]) || length(CGGP$design)==0) {
    0
  } else {
    nrow(CGGP$design)
  }
  max_polevels = apply(CGGP$po[1:CGGP$poCOUNT, ,drop=FALSE], 2, max)
  
  separateoutputparameterdimensions <- is.matrix(CGGP$thetaMAP)
  # nopd is numberofoutputparameterdimensions
  nopd <- if (separateoutputparameterdimensions) {
    if (length(CGGP$y)>0) {ncol(CGGP$y)} else {ncol(CGGP$ys)}
  } else {
    1
  }
  
  
  # ==============================.
  # ====    Calculate IMSE    ====
  # ==============================.
  
  # Calculate integrated mean squared error (IMSE) values for the given method
  if(selectionmethod=="MAP"){
    # Set up blank array to store MSE values
    MSE_MAP = array(0, dim=c(CGGP$d, CGGP$maxlevel,nopd))
    
    # Loop over dimensions and design refinements
    for (opdlcv in 1:nopd) {
      thetaMAP.thisloop <- if (nopd==1) CGGP$thetaMAP else CGGP$thetaMAP[, opdlcv]
      for (dimlcv in 1:CGGP$d) {
        for (levellcv in 1:max_polevels[dimlcv]) {
          # Calculate some sort of MSE from above, not sure what it's doing
          MSE_MAP[dimlcv, levellcv, opdlcv] = 
            max(0, 
                abs(
                  CGGP_internal_calcMSE(
                    CGGP$xb[1:CGGP$sizest[levellcv]],
                    thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
                    CGGP$CorrMat
                  )
                )
            )
          if (levellcv > 1.5) { # If past 1st level, it is as good as one below
            MSE_MAP[dimlcv, levellcv, opdlcv] = 
              min(MSE_MAP[dimlcv, levellcv, opdlcv], MSE_MAP[dimlcv, levellcv - 1, opdlcv])
          }
        }
      }
    }
    
    # Integrated MSE
    IMES_MAP = rep(0, CGGP$ML)
    
    # For all possible blocks, calculate MSE_MAP, need to apply it over nopd
    IMES_MAP_beforemean = apply(MSE_MAP, 3,
                                function(x) {
                                  CGGP_internal_calcMSEde(
                                    CGGP$po[1:CGGP$poCOUNT, , drop=F],
                                    x)
                                })
    if (CGGP$poCOUNT==1) {
      IMES_MAP_beforemean <- matrix(IMES_MAP_beforemean, nrow=1)
    }
    if (!is.matrix(IMES_MAP_beforemean)) {stop("Need a matrix here 0923859")}
    # Need as.matrix in case of single value
    #   i.e. when only supp data and only po is initial point
    # If multiple output but single opd, need to take mean
    sigma2MAP.thisloop <- if (nopd==1) {
      mean(CGGP$sigma2MAP)
    } else {
      CGGP$sigma2MAP
    }
    IMES_MAP[1:CGGP$poCOUNT] = rowMeans(
      sweep(IMES_MAP_beforemean, 2,
            sigma2MAP.thisloop, "*")
    )
    
    # Clean up to avoid silly errors
    rm(opdlcv, thetaMAP.thisloop, sigma2MAP.thisloop)
  } else if (selectionmethod %in% c("UCB", "TS")) { # selectionmethod is UCB or TS
    MSE_PostSamples = array(0, c(CGGP$d, CGGP$maxlevel,CGGP$numPostSamples, nopd))
    # Dimensions can be considered independently
    # Loop over dimensions and design refinements
    for (opdlcv in 1:nopd) { # Loop over output parameter dimensions
      thetaPostSamples.thisloop <- if (nopd==1) {
        CGGP$thetaPostSamples
      } else {
        CGGP$thetaPostSamples[ , , opdlcv]
      }
      for (dimlcv in 1:CGGP$d) { # Loop over each input dimension
        for (levellcv in 1:max_polevels[dimlcv]) {
          for(samplelcv in 1:CGGP$numPostSamples){
            # Calculate some sort of MSE from above, not sure what it's doing
            MSE_PostSamples[dimlcv, levellcv,samplelcv, opdlcv] = 
              max(0, 
                  abs(
                    CGGP_internal_calcMSE(
                      CGGP$xb[1:CGGP$sizest[levellcv]],
                      thetaPostSamples.thisloop[(dimlcv-1)*CGGP$numpara +
                                                  1:CGGP$numpara,
                                                samplelcv],
                      CGGP$CorrMat)
                  )
              )
            if (levellcv > 1.5) { # If past first level, it is as good as one below it
              MSE_PostSamples[dimlcv, levellcv,samplelcv, opdlcv] = 
                min(MSE_PostSamples[dimlcv, levellcv,samplelcv, opdlcv],
                    MSE_PostSamples[dimlcv, levellcv - 1,samplelcv, opdlcv])
            }
          }
        }
      }
    }
    rm(opdlcv, dimlcv, levellcv, samplelcv) # Avoid dumb mistakes
    IMES_PostSamples = matrix(0, CGGP$ML,CGGP$numPostSamples)
    
    # Calculate sigma2 for all samples if needed
    if (is.null(CGGP$sigma2_samples)) {
      CGGP$sigma2_samples <- CGGP_internal_calc_sigma2_samples(CGGP)
    }
    sigma2.allsamples.alloutputs <- CGGP$sigma2_samples
    
    for(samplelcv in 1:CGGP$numPostSamples){
      if (nopd == 1) { # Will be a matrix
        # Multiply by sigma2. If multiple output dimensions with
        #  shared parameters, take mean.
        #  Needed because each thetasample will have a different sigma2.
        sigma2.thistime <- mean(sigma2.allsamples.alloutputs[samplelcv,])
        IMES_PostSamples[1:CGGP$poCOUNT,samplelcv] = sigma2.thistime *
          CGGP_internal_calcMSEde(CGGP$po[1:CGGP$poCOUNT,],
                                  MSE_PostSamples[,,samplelcv,])
        rm(sigma2.thistime) # Avoid mistakes
      } else { # Is a 3d array, need to use an apply and then apply again with mean
        IMES_PostSamples_beforemean <- 
          apply(MSE_PostSamples[,,samplelcv,], 3,
                function(x){
                  CGGP_internal_calcMSEde(CGGP$po[1:CGGP$poCOUNT,,drop=F], x)
                })
        if (!is.matrix(IMES_PostSamples_beforemean)) {
          # Happens when CGGP$poCOUNT is 1, when only initial block avail
          if (CGGP$poCOUNT!=1) {stop("Something is wrong here #279287522")}
          IMES_PostSamples_beforemean <- matrix(IMES_PostSamples_beforemean, nrow=1)
        }
        # Need sigma2 for this theta sample, already calculated in sigma2.allsamples.alloutputs
        IMES_PostSamples[1:CGGP$poCOUNT,samplelcv] <-
          apply(IMES_PostSamples_beforemean, 1,
                function(x) {
                  # Weight by sigma2 samples
                  mean(sigma2.allsamples.alloutputs[samplelcv,] *
                         x)
                })
      }
    }; rm(samplelcv)
    # Get UCB IMES using 90% upper conf bound
    IMES_UCB = numeric(CGGP$ML)
    IMES_UCB[1:CGGP$poCOUNT] = apply(IMES_PostSamples[1:CGGP$poCOUNT,, drop=F],1,quantile, probs=0.9)
  } else {
    # Can be Oldest or Random or Lowest
  }
  
  
  
  # =============================.
  # ====    Append points    ====
  # =============================.
  
  # Append points to design until limit until reaching max_design_points
  max_design_points = CGGP$ss + batchsize
  while (max_design_points > CGGP$ss + min(CGGP$pogsize[1:CGGP$poCOUNT]) - .5) {
    if(selectionmethod=="MAP"){
      IMES = IMES_MAP
    } else if(selectionmethod=="UCB"){
      IMES = IMES_UCB
    } else if(selectionmethod=="TS"){
      IMES = IMES_PostSamples[,sample(1:CGGP$numPostSamples,1)]
    } else if(selectionmethod=="Oldest"){
      IMES = seq.int(from=CGGP$poCOUNT, to=1, by=-1)
      # Multiply by size so it gets undone below
      IMES <- IMES * CGGP$pogsize[1:CGGP$poCOUNT]
    } else if(selectionmethod=="Random"){
      IMES = rep(1,CGGP$poCOUNT)
      # Multiply by size so it gets undone below
      IMES <- IMES * CGGP$pogsize[1:CGGP$poCOUNT]
    } else if(selectionmethod=="Lowest"){
      IMES = rowSums(CGGP$po[1:CGGP$poCOUNT,])
      # Make the lowest the highest value
      IMES <- max(IMES) + 1 - IMES
      # Multiply by size so it gets undone below
      IMES <- IMES * CGGP$pogsize[1:CGGP$poCOUNT]
    } else {
      stop("Selection method not acceptable")
    }
    CGGP$uoCOUNT = CGGP$uoCOUNT + 1 #increment used count
    
    # Find which blocks are still valid for selecting
    stillpossible <- which(CGGP$pogsize[1:CGGP$poCOUNT] < 
                             (max_design_points - CGGP$ss + 0.5))
    
    # Pick block with max IMES per point in the block
    metric <- IMES[1:CGGP$poCOUNT] / CGGP$pogsize[1:CGGP$poCOUNT]
    # Find the best one that still fits
    M_comp = max(metric[stillpossible])
    # Find which ones are close to M_comp and pick randomly among them
    possibleO = stillpossible[metric[stillpossible] >= 0.99*M_comp]
    
    
    # If more than one is possible and near the best, randomly pick among them.
    if(length(possibleO)>1.5){
      pstar = sample(possibleO,1)
    } else{
      pstar = possibleO
    }
    
    l0 =  CGGP$po[pstar, ] # Selected block
    # Need to make sure there is still an open row in uo to set with new values
    if (CGGP$uoCOUNT > nrow(CGGP$uo)) {
      CGGP <- CGGP_internal_addrows(CGGP)
    }
    CGGP$uo[CGGP$uoCOUNT,] = l0 # Save selected block
    CGGP$ss =  CGGP$ss + CGGP$pogsize[pstar] # Update selected size
    
    
    # ================================.
    # ====    Update ancestors    ====
    # ================================.
    
    # Protect against initial block which has no ancestors
    if (CGGP$pilaCOUNT[pstar] > 0) { # Protect for initial block
      new_an = CGGP$pila[pstar, 1:CGGP$pilaCOUNT[pstar]]
      total_an = new_an
      for (anlcv in 1:length(total_an)) { # Loop over ancestors
        if (total_an[anlcv] > 1.5) { # If there's more than 1, do this
          total_an = unique(
            c(total_an,
              CGGP$uala[total_an[anlcv], 1:CGGP$ualaCOUNT[total_an[anlcv]]])
          )
        }
      }
      
      CGGP$ualaCOUNT[CGGP$uoCOUNT]  = length(total_an)
      CGGP$uala[CGGP$uoCOUNT, 1:length(total_an)] = total_an
      
      # Loop over all ancestors, update weight
      for (anlcv in 1:length(total_an)) {
        lo = CGGP$uo[total_an[anlcv],]
        if (max(abs(lo - l0)) < 1.5) {
          CGGP$w[total_an[anlcv]] = CGGP$w[total_an[anlcv]] + (-1)^abs(round(sum(l0-lo)))
        }
      }
    }
    CGGP$w[CGGP$uoCOUNT] = CGGP$w[CGGP$uoCOUNT] + 1
    
    
    # Update data. Remove selected item, move rest up.
    # First get correct indices to change. Protect when selecting initial point
    new_indices <- if (CGGP$poCOUNT>1) {1:(CGGP$poCOUNT - 1)} else {numeric(0)}
    old_indices <- setdiff(seq.int(1, CGGP$poCOUNT, 1), pstar)
    # Then change the data
    CGGP$po[new_indices,] = CGGP$po[old_indices,]
    CGGP$pila[new_indices,] = CGGP$pila[old_indices,]
    CGGP$pilaCOUNT[new_indices] = CGGP$pilaCOUNT[old_indices]
    CGGP$pogsize[new_indices] = CGGP$pogsize[old_indices]
    if(selectionmethod=="MAP"){
      IMES_MAP[new_indices] = IMES_MAP[old_indices]
    }
    if(selectionmethod=="UCB"){
      IMES_UCB[new_indices] = IMES_UCB[old_indices]
    }
    if(selectionmethod=="TS"){
      IMES_PostSamples[new_indices,] = IMES_PostSamples[old_indices,]
    }
    # And reduce number of available blocks by one.
    CGGP$poCOUNT = CGGP$poCOUNT - 1
    
    # ==========================================.
    # ====    Update new possible blocks    ====
    # ==========================================.
    
    # Loop over possible descendents of selected block, add them if possible    
    for (dimlcv in 1:CGGP$d) {
      lp = l0
      
      lp[dimlcv] = lp[dimlcv] + 1
      
      if (max(lp) <= CGGP$maxlevel && CGGP$poCOUNT < 4 * CGGP$ML) {
        kvals = which(lp > 1.5) # Dimensions above base level
        
        canuse = 1
        ap = rep(0, CGGP$d)
        nap = 0
        for (activedimlcv in 1:length(kvals)) {
          lpp = lp
          lpp[kvals[activedimlcv]] = lpp[kvals[activedimlcv]] - 1
          
          ismem = rep(1, CGGP$uoCOUNT)
          for (dimdimlcv in 1:CGGP$d) {
            ismem  = ismem *
              (CGGP$uo[1:CGGP$uoCOUNT, dimdimlcv] == lpp[dimdimlcv])
          }
          
          if (max(ismem) > 0.5) {
            ap[activedimlcv] = which(ismem > 0.5)
            nap = nap + 1
          } else{
            canuse = 0
          }
        }
        if (canuse > 0.5) { # If it can be used, add to possible blocks
          CGGP$poCOUNT = CGGP$poCOUNT + 1
          CGGP$po[CGGP$poCOUNT,] = lp
          CGGP$pogsize[CGGP$poCOUNT] = prod(CGGP$sizes[lp])
          CGGP$pila[CGGP$poCOUNT, 1:nap] = ap[1:nap]
          CGGP$pilaCOUNT[CGGP$poCOUNT] = nap
          
          max_polevels_old = max_polevels
          max_polevels = apply(CGGP$po[1:CGGP$poCOUNT, ,drop=F], 2, max)
          
          if(selectionmethod=="MAP"){
            for (opdlcv in 1:nopd) { # Loop over output parameter dimensions
              thetaMAP.thisloop <- if (nopd==1) CGGP$thetaMAP else CGGP$thetaMAP[, opdlcv]
              for (dimlcv in 1:CGGP$d) {
                if((max_polevels_old[dimlcv]+0.5)<max_polevels[dimlcv]){
                  levellcv = max_polevels[dimlcv]
                  MSE_MAP[dimlcv, levellcv,
                          opdlcv] = max(0, abs(CGGP_internal_calcMSE(CGGP$xb[1:CGGP$sizest[levellcv]],
                                                                     thetaMAP.thisloop[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara],
                                                                     CGGP$CorrMat)))
                  if (levellcv > 1.5) { # If past first level, it is as good as one below it. Why isn't this a result of calculation?
                    MSE_MAP[dimlcv, levellcv, opdlcv] = min(MSE_MAP[dimlcv, levellcv, opdlcv], MSE_MAP[dimlcv, levellcv - 1, opdlcv])
                  }
                }
              }
            }
            # Clean up
            rm(thetaMAP.thisloop, opdlcv)
          } else if (selectionmethod %in% c("UCB", "TS")){ # selection method is UCB or TS
            for (opdlcv in 1:nopd) {
              thetaPostSamples.thisloop <- if (nopd==1) CGGP$thetaPostSamples else CGGP$thetaPostSamples[, , opdlcv]
              for (dimlcv_2 in 1:CGGP$d) { # dimlcv is already used for which descendent to add
                if((max_polevels_old[dimlcv_2]+0.5)<max_polevels[dimlcv_2]){
                  levellcv = max_polevels[dimlcv_2]
                  for(samplelcv in 1:CGGP$numPostSamples){
                    # Calculate some sort of MSE from above, not sure what it's doing
                    MSE_PostSamples[dimlcv_2, levellcv,
                                    samplelcv, opdlcv] = max(0,
                                                             abs(CGGP_internal_calcMSE(
                                                               CGGP$xb[1:CGGP$sizest[levellcv]],
                                                               thetaPostSamples.thisloop[(dimlcv_2-1)*CGGP$numpara+1:CGGP$numpara,
                                                                                         samplelcv],
                                                               CGGP$CorrMat)))
                    if (levellcv > 1.5) { # If past first level, it is as good as one below it. Why isn't this a result of calculation?
                      MSE_PostSamples[dimlcv_2, levellcv,
                                      samplelcv, opdlcv] = min(MSE_PostSamples[dimlcv_2, levellcv,samplelcv, opdlcv],
                                                               MSE_PostSamples[dimlcv_2, levellcv - 1,samplelcv, opdlcv])
                    }
                  }; rm(samplelcv)
                }
              }; rm(dimlcv_2)
            }
            # Clean up
            rm(thetaPostSamples.thisloop, opdlcv)
          } else {
            # Can be Oldest or Random or Lowest
          }
          
          if(selectionmethod=="MAP"){
            # IMES_MAP[CGGP$poCOUNT] = CGGP_internal_calcMSEde(as.vector(CGGP$po[CGGP$poCOUNT, ]), MSE_MAP)
            # Need to apply first
            IMES_MAP_beforemeannewpoint <- apply(MSE_MAP, 3,
                                                 function(x) {CGGP_internal_calcMSEde(as.vector(CGGP$po[CGGP$poCOUNT, ]), x)})
            # Take weighted mean over dimensions
            IMES_MAP[CGGP$poCOUNT] <- mean(CGGP$sigma2MAP * IMES_MAP_beforemeannewpoint)
          } else if (selectionmethod=="UCB" || selectionmethod=="TS"){
            for(samplelcv in 1:CGGP$numPostSamples){
              if (nopd == 1) { # is a matrix
                # Each sample has different sigma2, so use. If multiple output
                #  parameter dimensions, take mean over sigma2.
                sigma2.thistime <- mean(sigma2.allsamples.alloutputs[samplelcv,])
                IMES_PostSamples[CGGP$poCOUNT,samplelcv] = sigma2.thistime *
                  CGGP_internal_calcMSEde(as.vector(CGGP$po[CGGP$poCOUNT, ]),
                                          MSE_PostSamples[,,samplelcv,])
                rm(sigma2.thistime)
              } else { # is an array, need to apply
                IMES_PostSamples_beforemeannewpoint = apply(MSE_PostSamples[,,samplelcv,],
                                                            3, # 3rd dim since samplelcv removes 3rd
                                                            function(x) {
                                                              CGGP_internal_calcMSEde(as.vector(CGGP$po[CGGP$poCOUNT, ]), x)
                                                            }
                )
                IMES_PostSamples[CGGP$poCOUNT,samplelcv] <- mean(sigma2.allsamples.alloutputs[samplelcv,] * 
                                                                   IMES_PostSamples_beforemeannewpoint)
              }
            }; rm(samplelcv)
            IMES_UCB[CGGP$poCOUNT] = quantile(IMES_PostSamples[CGGP$poCOUNT,],probs=0.9)
          } else if (selectionmethod %in% c("Oldest", "Random", "Lowest")) {
            # nothing needed
          } else {stop("Not possible #9235058")}
        }
      }
    }
  }
  
  
  # Get design and other attributes updated
  CGGP <- CGGP_internal_getdesignfromCGGP(CGGP)
  
  
  # Check if none were added, return warning/error
  if (n_before == nrow(CGGP$design)) {
    warning("No points could be added. You may need a larger batch size.")
  } else {
    # Save design_unevaluated to make it easy to know which ones to add
    CGGP$design_unevaluated <- CGGP$design[(n_before+1):nrow(CGGP$design),]
  }
  
  return(CGGP)
}
