#' Create sparse grid GP
#'
#' @param d Input dimension
# @param xmin Min x values, vector. Must be rep(0,d).
# @param xmax Max x values, vector. Must be rep(1,d).
#' @param batchsize Number added to design each batch
# @param nugget Nugget term added to diagonal of correlation matrix,
#' for now only on predictions
#' @param corr Name of correlation function to use. Must be one of "CauchySQT",
#' "CauchySQ", "Cauchy", "Gaussian", "PowerExp", "Matern32", "Matern52".
#' @param grid_sizes Size of grid refinements.
#' @param Xs Supplemental X data
#' @param Ys Supplemental Y data
#' @param supp_args Arguments used to fit if Xs and Ys are given
#' @param HandlingSuppData How should supplementary data be handled?
#' * Correct: full likelihood with grid and supplemental data
#' * Only: only use supplemental data
#' * Ignore: ignore supplemental data
#' 
#' @importFrom stats rbeta
#'
#' @return CGGP
#' @export
#' @family CGGP core functions
#'
#' @examples
#' CGGPcreate(d=8,200)
CGGPcreate <- function(d, batchsize, corr="PowerExponential",
                       grid_sizes=c(1,2,4,4,8,12,20,28,32),
                       Xs=NULL, Ys=NULL,
                       HandlingSuppData="Correct",
                       supp_args=list()
) {
  if (d < 2) {stop("d must be at least 2")}
  
  # ==================================.
  # ====    Create CGGP object    ====
  # ==================================.
  
  # This is list representing our GP object
  CGGP <- list()
  class(CGGP) <- c("CGGP", "list") # Give it class CGGP
  CGGP$d <- d
  CGGP$numPostSamples <- 100
  CGGP$HandlingSuppData <- HandlingSuppData
  CGGP <- CGGP_internal_set_corr(CGGP, corr)
  CGGP$nugget <- 0
  
  
  # Partial matching is very bad! Keep these as length 0 instead of NULL,
  #  otherwise CGGP$Y can return CGGP$Ys
  CGGP$Y <- numeric(0)
  CGGP$y <- numeric(0)
  
  # ====================================================.
  # ==== If supplemental data is given, fit it here ====
  # ====================================================.
  if (!is.null(Xs) && !is.null(Ys)) {
    if (!is.null(supp_args) && length(supp_args) > 0 &&
        is.null(names(supp_args))) {
      stop("Give names for supp_args")
    }
    supp_args$CGGP <- CGGP
    supp_args$Xs <- Xs
    supp_args$Ys <- Ys
    CGGP <- do.call(CGGP_internal_fitwithonlysupp, supp_args)
  }
  
  # ===========================================.
  # ====    Start setting up CGGP stuff    ====
  # ===========================================.
  
  # Levels are blocks. Level is like eta from paper.
  CGGP$ML = min(choose(CGGP$d + 6, CGGP$d), 10000) #max levels
  
  # Track evaluated blocks, aka used levels
  CGGP$uo = matrix(0, nrow = CGGP$ML, ncol = CGGP$d) # blocks that have been selected
  CGGP$uoCOUNT = 0 # number of selected blocks
  
  # Track the blocks that are allowed to be evaluated
  CGGP$po = matrix(0, nrow = 4 * CGGP$ML, ncol = CGGP$d) #proposed levels tracker
  # Only option at first is initial block (1,1,...,1)
  CGGP$po[1, ] <- rep(1, CGGP$d)
  CGGP$poCOUNT <- 1
  
  # Ancestors are blocks one level down in any dimension.
  CGGP$maxgridsize = 400
  CGGP$pila = matrix(0, nrow = CGGP$ML, ncol =CGGP$maxgridsize ) #proposed immediate level ancestors
  CGGP$pala = matrix(0, nrow = CGGP$ML, ncol =CGGP$maxgridsize ) #proposedal all level ancestors
  CGGP$uala = matrix(0, nrow = CGGP$ML, ncol =CGGP$maxgridsize ) #used all level ancestors
  CGGP$pilaCOUNT = rep(0, CGGP$ML) #count of number of pila
  CGGP$palaCOUNT = rep(0, CGGP$ML) #count of number of pala
  CGGP$ualaCOUNT = rep(0, CGGP$ML) #count of number of uala
  
  # Initial block (1,1,...,1) has no ancestors
  CGGP$pilaCOUNT[1] <- 0
  CGGP$pila[1, 1] <- 0
  
  # CGGP$sizes = c(1,2,4,4,8,12,32) # Num of points added to 1D design as you go further in any dimension
  CGGP$sizes <- grid_sizes
  CGGP$maxlevel = length(CGGP$sizes)
  # Proposed grid size
  CGGP$pogsize = rep(0, 4 * CGGP$ML)
  CGGP$pogsize[1:CGGP$poCOUNT] = apply(matrix(CGGP$sizes[CGGP$po[1:CGGP$poCOUNT, ]], CGGP$poCOUNT, CGGP$d), 1, prod)
  # Selected sample size
  CGGP$ss = 0
  
  CGGP$w = rep(0, CGGP$ML) #keep track of + and - for prediction
  CGGP$uoCOUNT = 0 ###1 # Number of used levels
  
  
  # =========================.
  # ====   Add Blocks    ====
  # =========================.
  
  # sample has unexpected behavior, e.g., sample(34,1), see help file for sample
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  
  # While number selected + min sample size <= batch size, i.e.,
  #  still have enough spots for a block, keep adding blocks
  while (batchsize > (CGGP$ss + min(CGGP$pogsize[1:CGGP$poCOUNT]) - 0.5)) {
    CGGP$uoCOUNT = CGGP$uoCOUNT + 1 #increment used count
    
    
    if (CGGP$uoCOUNT < 1.5) { # Nothing picked yet, so take base block (1,1,...,1)
      pstar <- 1
    } else if (CGGP$uoCOUNT < (CGGP$d + 1.5)) {
      # Next d iterations pick the (2,1,1,1,1),(1,2,1,1,1) blocks b/c we need
      #  info on each dimension before going adaptive
      pstar = 1
    } else{ # Next d iterations randomly pick from boxes w/ min # of pts
        criteriahere = rowSums(CGGP$po[1:CGGP$poCOUNT,])
        A1 = (CGGP$pogsize[1:CGGP$poCOUNT] < min(batchsize - CGGP$ss + 0.5,CGGP$maxgridsize)) 
        MCN = min(criteriahere[A1])
        A2 = criteriahere <= 0.5 + MCN
        pstar <- resample(which(A1 & A2), 1)
    }
    
    l0 =  CGGP$po[pstar, ] # Selected block e.g. (2,1,1,2)
    # Need to make sure there is still an open row in uo to set with new values
    if (CGGP$uoCOUNT > nrow(CGGP$uo)) {
      CGGP <- CGGP_internal_addrows(CGGP)
    }
    CGGP$uo[CGGP$uoCOUNT, ] = l0 # Store new block
    CGGP$ss =  CGGP$ss + CGGP$pogsize[pstar] # Update selected sample size
    
    # Ancestors of block just selected
    # Need to give possibility for initial block, has no ancestors, and 1:0 is bad
    new_an = CGGP$pila[pstar,
                       if (CGGP$pilaCOUNT[pstar]>.5) {1:CGGP$pilaCOUNT[pstar]}
                       else {numeric(0)}]
    total_an = new_an
    
    # Loop over ancestors of block just selected
    if (length(total_an) > .5) { # Initial block has no total_an , need this to avoid 1:0
      for (anlcv in 1:length(total_an)) {
        # If more than one ancestor, update with unique ones.
        if (total_an[anlcv] > 1.5) {
          total_an = unique(c(total_an,
                              CGGP$uala[total_an[anlcv],
                                        1:CGGP$ualaCOUNT[total_an[anlcv]]]))
        }
      }
      # Update storage of ancestors
      CGGP$ualaCOUNT[CGGP$uoCOUNT]  = length(total_an)
      CGGP$uala[CGGP$uoCOUNT, 1:length(total_an)] = total_an
    }
    
    # Loop over ancestors, update coefficient
    if (length(total_an) > .5) { # Initial block has no total_an , need this to avoid 1:0
      for (anlcv in 1:length(total_an)) {
        lo = CGGP$uo[total_an[anlcv], ]
        if (max(abs(lo - l0)) < 1.5) {
          CGGP$w[total_an[anlcv]] = CGGP$w[total_an[anlcv]] + (-1) ^ abs(round(sum(l0-lo)))
          
        }
      }
    }
    CGGP$w[CGGP$uoCOUNT] = CGGP$w[CGGP$uoCOUNT] + 1
    
    CGGP$po[pstar,] <- 0 # Clear just used row
    if (CGGP$poCOUNT > 1.5) { # Move up other options if there are others left
      po_rows_to_move <- (1:CGGP$poCOUNT)[-pstar] # Moving all po rows except selected
      CGGP$po[1:(CGGP$poCOUNT - 1), ] = CGGP$po[po_rows_to_move, ]
      CGGP$pila[1:(CGGP$poCOUNT - 1), ] = CGGP$pila[po_rows_to_move, ]
      CGGP$pilaCOUNT[1:(CGGP$poCOUNT - 1)] = CGGP$pilaCOUNT[po_rows_to_move]
      CGGP$pogsize[1:(CGGP$poCOUNT - 1)] = CGGP$pogsize[po_rows_to_move]
    }
    
    # One less option now
    CGGP$poCOUNT = CGGP$poCOUNT - 1
    
    # Loop over dimensions to add new possible blocks
    for (dimlcv in 1:CGGP$d) {
      # The block e.g. (1,2,1,1,3) just selected
      lp = l0
      # Increase single dimension by 1, will see if it is possible
      lp[dimlcv] = lp[dimlcv] + 1
      
      # Check if within bounds
      if (max(lp) <= CGGP$maxlevel && CGGP$poCOUNT < 4*CGGP$ML) {
        # Dimensions which are past first design level
        kvals = which(lp > 1.5)
        
        canuse = 1 # Can this block be used? Will be set to 0 below if not.
        ap = rep(0, CGGP$d) # Ancestors
        nap = 0 # Number ancestors
        
        # Loop over dims at 2+
        for (activedimlcv in 1:length(kvals)) {
          lpp = lp # The block selected with 1 dim incremented
          lpp[kvals[activedimlcv]] = lpp[kvals[activedimlcv]] - 1
          
          ismem = rep(1, CGGP$uoCOUNT) # Boolean
          # Loop over dimensions
          for (dimdimlcv in 1:CGGP$d) {
            # Set to 0 or 1 if all points already selected have same value
            ismem  = ismem * (CGGP$uo[1:CGGP$uoCOUNT,
                                      dimdimlcv] == lpp[dimdimlcv])
          }
          # If any are still 1,
          if (max(ismem) > 0.5) {
            ap[activedimlcv] = which(ismem > 0.5)
            nap = nap + 1 # Count number that are >=1
          } else{ # All are 0, so can't use
            canuse = 0
          }
        }
        # If it can be used, add to possible blocks
        if (canuse > 0.5) {
          CGGP$poCOUNT = CGGP$poCOUNT + 1
          CGGP$po[CGGP$poCOUNT, ] = lp
          CGGP$pogsize[CGGP$poCOUNT] = prod(CGGP$sizes[lp])
          CGGP$pila[CGGP$poCOUNT, 1:nap] = ap[1:nap]
          CGGP$pilaCOUNT[CGGP$poCOUNT] = nap
        }
      }
    }
  }
  
  # Create points for design
  #  These are distances from the center 0.5.
  xb = rep(
    c(
      3 / 8, # 1/8, 7/8
      1 / 2, # 0, 1
      1 / 4, # 1/4, 3/4
      1 / 8, # 3/8, 5/8
      15 / 32, # etc
      7 / 16,
      3 / 16,
      5 / 16,
      7 / 32,
      11 / 32,
      3 / 32,
      13 / 32,
      9 / 32,
      5 / 32,
      1 / 32,
      1 / 16,
      seq(31,1,-4)/64,
      seq(29,1,-4)/64,
      seq(63,1,-4)/128,
      seq(61,1,-4)/128
    ),
    "each" = 2
  )
  CGGP$xb = 0.5 + c(0, xb * rep(c(-1, 1), length(xb) / 2))
  CGGP$xindex = 1:length(xb) # Why not length(CGGP$xb), which is one longer than xb?
  # After this xb is
  #  [1] 0.50000 0.12500 0.87500 0.25000 0.75000 0.37500 0.62500 0.28125
  #      0.71875 0.31250 0.68750 0.00000 1.00000 0.18750 0.81250
  # [16] 0.06250 0.93750 0.43750 0.56250 0.40625 0.59375 0.09375 0.90625
  #      0.21875 0.78125 0.34375 0.65625 0.46875 0.53125 0.15625
  # [31] 0.84375 0.03125 0.96875
  CGGP$sizest = cumsum(CGGP$sizes) # Total # of points in 1D design along axis
  
  
  # ======================================.
  # ==== Get design and return object ====
  # ======================================.
  # This is all to create design from uo.
  # If only supp data is given, don't run it.
  if (CGGP$uoCOUNT > 0) {
    # Get design from uo and other data
    CGGP <- CGGP_internal_getdesignfromCGGP(CGGP)
    CGGP$design_unevaluated <- CGGP$design
  }
  
  return(CGGP)
}