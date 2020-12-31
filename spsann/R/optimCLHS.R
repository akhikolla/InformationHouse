#' Optimization of sample configurations for spatial trend identification and estimation (IV)
#'
#' Optimize a sample configuration for spatial trend identification and estimation using the method proposed 
#' by Minasny and McBratney (2006), known as the conditioned Latin hypercube sampling. An utility function 
#' _U_ is defined so that the sample reproduces the marginal distribution and correlation matrix of the
#' numeric covariates, and the class proportions of the factor covariates (__CLHS__). The utility function 
#' is obtained aggregating three objective functions: __O1__, __O2__, and __O3__.
#' 
# @inheritParams spJitter
#' @template spSANN_doc
#' @inheritParams optimACDC
#' @template spJitter_doc
#' 
#' @param clhs.version (Optional) Character value setting the CLHS version that should be used. Available
#' options are: `"paper"`, for the formulations of __O1__, __O2__, and __O3__ as presented in the original 
#' paper by Minasny and McBratney (2006); `"fortran"`, for the formulations of __O1__ and __O3__ that include
#' a scaling factor as implemented in the late Fortran code by Budiman Minasny (ca. 2015); and `"update"`, for
#' formulations of __O1__, __O2__, and __O3__ that include the modifications proposed the authors of this
#' package in 2018 (see below). Defaults to `clhs.version = "paper"`.
#' 
#' @details 
#' \subsection{Marginal sampling strata}{
#' Reproducing the marginal distribution of the numeric covariates depends upon the definition of marginal 
#' sampling strata. _Equal-area_ marginal sampling strata are defined using the sample quantiles estimated 
#' with \code{\link[stats]{quantile}} using a continuous function (`type = 7`), that is, a function that
#' interpolates between existing covariate values to estimate the sample quantiles. This is the procedure
#' implemented in the original method of Minasny and McBratney (2006), which creates breakpoints that do not
#' occur in the population of existing covariate values. Depending on the level of discretization of the
#' covariate values, that is, how many significant digits they have, this can create repeated breakpoints, 
#' resulting in empty marginal sampling strata. The number of empty marginal sampling strata will ultimately
#' depend on the frequency distribution of the covariate and on the number of sampling points. The effect of
#' these features on the spatial modelling outcome still is poorly understood.
#' }
#' \subsection{Correlation between numeric covariates}{
#' The _correlation_ between two numeric covariates is measured using the sample Pearson's _r_, a descriptive
#' statistic that ranges from -1 to +1. This statistic is also known as the sample linear correlation 
#' coefficient. The effect of ignoring the correlation among factor covariates and between factor and numeric
#' covariates on the spatial modelling outcome still is poorly understood.
#' }
#' \subsection{Multi-objective combinatorial optimization}{
#' A method of solving a multi-objective combinatorial optimization problem (MOCOP) is to aggregate the 
#' objective functions into a single utility function _U_. In the __spsann__ package, as in the original
#' implementation of the CLHS by Minasny and McBratney (2006), the aggregation is performed using the 
#' __weighted sum method__, which uses weights to incorporate the __a priori__ preferences of the user about 
#' the relative importance of each objective function. When the user has no preference, the objective functions 
#' receive equal weights.
#' 
#' The weighted sum method is affected by the relative magnitude of the different objective function values. 
#' The objective functions implemented in `optimCLHS` have different units and orders of magnitude. The 
#' consequence is that the objective function with the largest values, generally __O1__, may have a numerical
#' dominance during the optimization. In other words, the weights may not express the true preferences of the
#' user, resulting that the meaning of the utility function becomes unclear because the optimization will 
#' likely favour the objective function which is numerically dominant.
#' 
#' An efficient solution to avoid numerical dominance is to scale the objective functions so that they are 
#' constrained to the same approximate range of values, at least in the end of the optimization. In the 
#' original implementation of the CLHS by Minasny and McBratney (2006), `clhs.version = "paper"`, `optimCLHS`
#' uses the naive aggregation method, which ignores that the three objective functions have different units 
#' and orders of magnitude. In a 2015 Fortran implementation of the CLHS, `clhs.version = "fortran"`, scaling
#' factors were included to make the values of the three objective function more comparable. The effect of
#' ignoring the need to scale the objective functions, or using arbitrary scaling factors, on the spatial
#' modelling outcome still is poorly understood. Thus, an updated version of __O1__, __O2__, and __O3__ has
#' been implemented in the __spsann__ package. The need formulation aim at making the values returned by the 
#' objective functions more comparable among themselves without having to resort to arbitrary scaling factors.
#' The effect of using these new formulations have not been tested yet.
#' }
#' 
#' @return
#' `optimCLHS` returns an object of class `OptimizedSampleConfiguration`: the optimized sample configuration
#' with details about the optimization.
#' 
#' `objCLHS` returns a numeric value: the energy state of the sample configuration -- the objective function
#' value.
#' 
#' @references
#' Minasny, B.; McBratney, A. B. A conditioned Latin hypercube method for sampling in the presence of 
#' ancillary information. _Computers & Geosciences_, v. 32, p. 1378-1388, 2006.
#'
#' Minasny, B.; McBratney, A. B. Conditioned Latin Hypercube Sampling for calibrating soil sensor data to 
#' soil properties. Chapter 9. Viscarra Rossel, R. A.; McBratney, A. B.; Minasny, B. (Eds.) _Proximal Soil
#' Sensing_. Amsterdam: Springer, p. 111-119, 2010.
#'
#' Roudier, P.; Beaudette, D.; Hewitt, A. A conditioned Latin hypercube sampling algorithm incorporating
#' operational constraints. _5th Global Workshop on Digital Soil Mapping_. Sydney, p. 227-231, 2012.
#'
#' @note
#' The (only?) difference of `optimCLHS` to the original Fortran implementation of Minasny and McBratney
#' (2006), and to the `clhs` function implemented in the former 
#' __[clhs](https://CRAN.R-project.org/package=clhs)__ package by Pierre Roudier, is
#' the annealing schedule.
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[spsann]{optimACDC}}
#' @concept spatial trend
#' @aliases optimCLHS objCLHS CLHS
#' @export
#' @examples
#' data(meuse.grid, package = "sp")
#' candi <- meuse.grid[1:1000, 1:2]
#' covars <- meuse.grid[1:1000, 5]
#' schedule <- scheduleSPSANN(
#'   chains = 1, initial.temperature = 20, x.max = 1540, y.max = 2060, 
#'   x.min = 0, y.min = 0, cellsize = 40)
#' set.seed(2001)
#' res <- optimCLHS(
#'   points = 10, candi = candi, covars = covars, use.coords = TRUE,
#'   clhs.version = "fortran", weights = list(O1 = 0.5, O3 = 0.5), schedule = schedule)
#' objSPSANN(res) - objCLHS(
#'   points = res, candi = candi, covars = covars, use.coords = TRUE, 
#'   clhs.version = "fortran", weights = list(O1 = 0.5, O3 = 0.5))
# MAIN FUNCTION ###############################################################################################
optimCLHS <-
  function (points, candi,
            # O1, O2, and O3
            covars, use.coords = FALSE, clhs.version = c("paper", "fortran", "update"),
            # SPSANN
            schedule = scheduleSPSANN(), plotit = FALSE, track = FALSE,
            boundary, progress = "txt", verbose = FALSE,
            # MOOP
            weights) {
            # weights = list(O1 = 1/3, O2 = 1/3, O3 = 1/3)) {
    
    # Objective function name
    objective <- "CLHS"
    
    # Check spsann arguments
    eval(.check_spsann_arguments())
    
    # Check other arguments
    check <- .optimCLHScheck(candi = candi, covars = covars, use.coords = use.coords)
    if (!is.null(check)) { stop (check, call. = FALSE) }
    
    # Set plotting options
    eval(.plotting_options())
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare for jittering
    eval(.prepare_jittering())
    
    # Prepare 'covars' and base data
    eval(.prepare_clhs_covars())
    
    # Identify CLHS version
    clhs.version <- match.arg(clhs.version)

    # Compute initial energy state
    energy0 <- .objCLHS(
      sm = sm, breaks = breaks, id_num = id_num, pcm = pcm, id_fac = id_fac, n_pts = n_pts + n_fixed_pts, 
      pop_count = pop_count, n_candi = n_candi, weights = weights, covars_type = covars_type, 
      clhs.version = clhs.version)
    
    # Other settings for the simulated annealing algorithm
    old_sm <- sm
    new_sm <- sm
    best_sm <- sm
    old_energy <- energy0
    best_energy <- .bestEnergyCLHS(covars_type = covars_type)
    actual_temp <- schedule$initial.temperature
    k <- 0 # count the number of jitters
   
    # Set progress bar
    eval(.set_progress())
    
    # Initiate the annealing schedule
    for (i in 1:schedule$chains) {
      n_accept <- 0
      
      for (j in 1:schedule$chain.length) { # Initiate one chain
        
        for (wp in 1:n_pts) { # Initiate loop through points
          k <- k + 1
          
          # Plotting and jittering
          eval(.plot_and_jitter())
          
          # Update sample matrix and compute the new energy state
          new_sm[wp, ] <- covars[new_conf[wp, 1], ]
          new_energy <- .objCLHS(
            sm = new_sm, breaks = breaks, id_num = id_num, pcm = pcm, id_fac = id_fac, 
            n_pts = n_pts + n_fixed_pts, pop_count = pop_count, n_candi = n_candi, weights = weights, 
            covars_type = covars_type, clhs.version = clhs.version)
          
          # Evaluate the new system configuration
          accept <- .acceptSPSANN(old_energy[[1]], new_energy[[1]], actual_temp)
          if (accept) {
            old_conf <- new_conf
            old_energy <- new_energy
            old_sm <- new_sm
            n_accept <- n_accept + 1
          } else {
            new_energy <- old_energy
            new_conf <- old_conf
            new_sm <- old_sm
          }
          if (track) energies[k, ] <- new_energy
          
          # Record best energy state
          if (new_energy[[1]] < best_energy[[1]] / 1.0000001) {
            best_k <- k
            best_conf <- new_conf
            best_energy <- new_energy
            best_old_energy <- old_energy
            old_conf <- old_conf
            best_sm <- new_sm
            best_old_sm <- old_sm
          }
          
          # Update progress bar
          eval(.update_progress())
          
        } # End loop through points
        
      } # End the chain
      
      # Check the proportion of accepted jitters in the first chain
      eval(.check_first_chain())
      
      # Count the number of chains without any change in the objective function.
      # Restart with the previously best configuration if it exists.
      if (n_accept == 0) {
        no_change <- no_change + 1
        if (no_change > schedule$stopping) {
          # if (new_energy[[1]] > best_energy[[1]] * 1.000001) {
            # old_conf <- old_conf
            # new_conf <- best_conf
            # old_energy <- best_old_energy
            # new_energy <- best_energy
            # new_sm <- best_sm
            # old_sm <- best_old_sm
            # no_change <- 0
            # cat("\nrestarting with previously best configuration\n")
          # } else { 
            break 
          # }
        }
        if (verbose) {
          cat("\n", no_change, "chain(s) with no improvement... stops at", schedule$stopping, "\n")
        }
      } else {
        no_change <-  0
      }
      
      # Update control parameters
      actual_temp <- actual_temp * schedule$temperature.decrease
      x.max <- x_max0 - (i / schedule$chains) * (x_max0 - x.min) + cellsize[1]
      y.max <- y_max0 - (i / schedule$chains) * (y_max0 - y.min) + cellsize[2]
      
    } # End the annealing schedule
    
    # Prepare output
    eval(.prepare_output())
  }
# INTERNAL FUNCTION - CHECK ARGUMENTS #########################################################################
# candi: candidate locations
# covars: covariates
# use.coords: should the coordinates be used
.optimCLHScheck <-
  function (candi, covars, use.coords) {
    
    # covars
    if (is.vector(covars)) {
      if (use.coords == FALSE) {
        return ("'covars' must have two or more columns")
      }
      if (nrow(candi) != length(covars)) {
        return ("'candi' and 'covars' must have the same number of rows")
      }
    } else {
      if (nrow(candi) != nrow(covars)) {
        return ("'candi' and 'covars' must have the same number of rows")
      }
    }
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE ###########################################################
# This function is used to calculate the criterion value of CLHS.
# Aggregation is done using the weighted sum method.
.objCLHS <-
  function (sm, breaks, id_num, pcm, id_fac, n_pts, pop_count, n_candi, weights, covars_type, clhs.version) {
    
    # Objective functions
    if (any(covars_type == c("numeric", "both"))) {
      obj_O1 <- weights$O1 * .objO1(sm = sm, breaks = breaks, id_num = id_num, clhs.version = clhs.version)
      obj_O3 <- weights$O3 * .objO3(sm = sm, id_num = id_num, pcm = pcm, clhs.version = clhs.version)
    }
    if (any(covars_type == c("factor", "both"))) {
      obj_O2 <- weights$O2 * 
        # .objO2(sm = sm, id_fac = id_fac, n_pts = n_pts, pop_prop = pop_prop, clhs.version = clhs.version)
        .objO2(sm = sm, id_fac = id_fac, n_pts = n_pts, pop_count = pop_count, n_candi = n_candi, 
               clhs.version = clhs.version)
    }
    
    # Prepare output, a data.frame with the weighted sum in the first column followed by the values of the
    # constituent objective functions (IN ALPHABETICAL ORDER).
    if (covars_type == "both") {
      res <- data.frame(
        obj = obj_O1 + obj_O2 + obj_O3, 
        O1 = obj_O1, 
        O2 = obj_O2, 
        O3 = obj_O3)
    } else if (covars_type == "numeric") {
      res <- data.frame(
        obj = obj_O1 + obj_O3, 
        O1 = obj_O1, 
        O3 = obj_O3)
    } else {
      res <- data.frame(
        obj = obj_O2)
    }
    return (res)
    # } else {
    #   if (covars_type == "numeric") {
    #     return (data.frame(obj = obj_O1 + obj_O3, O1 = obj_O1, O3 = obj_O3))
    #   } else {
    #     return (data.frame(obj = obj_O2))
    #   }
    # }
  }
# CALCULATE OBJECTIVE FUNCTION VALUE ##########################################################################
#' @rdname optimCLHS
#' @export
objCLHS <-
  function (points, candi, covars, use.coords = FALSE, 
            clhs.version = c("paper", "fortran", "update"), 
            weights) {
            # weights = list(O1 = 1/3, O2 = 1/3, O3 = 1/3)) {
    
    # Check arguments
    check <- .optimCLHScheck(candi = candi, covars = covars, use.coords = use.coords)
    if (!is.null(check)) stop (check, call. = FALSE)
    
    # Prepare points and candi
    eval(.prepare_points())
    
    # Prepare 'covars' and and base data
    eval(.prepare_clhs_covars())
    
    # Identify CLHS version
    clhs.version <- match.arg(clhs.version)
    
    # Output energy state
    out <- .objCLHS(
      sm = sm, breaks = breaks, id_num = id_num, pcm = pcm, id_fac = id_fac, n_pts = n_pts, n_candi = n_candi,
      pop_count = pop_count, weights = weights, covars_type = covars_type, clhs.version = clhs.version)
    return(out)
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE (O1) ######################################################
# sm: sample matrix
# breaks: break points of the marginal sampling strata
# id_num: number of the column containing numeric covariates
# clhs.version: CLHS version
.objO1 <- 
  function (sm, breaks, id_num, clhs.version) {
    
    # Count the number of points per marginal sampling strata
    sm_count <- sapply(1:length(id_num), function (i) 
      graphics::hist(sm[id_num][, i], breaks[[i]], plot = FALSE)$counts)
    
    out <- switch (clhs.version,
      paper = {
        # Minasny and McBratney (2006)
        sum(abs(sm_count - 1))
      },
      fortran = {
        # The late FORTRAN code of Budiman Minasny -- ca. 2015 -- implements scaling factors so that values
        # are "more" comparable among objective functions. For O1, the scaling factor is defined as the number
        # of samples, nrow(sm), multiplied by the number of continuous variables, length(id_num), that is, the
        # total number of marginal sampling strata among all continuous variables.
        n <- nrow(sm) * length(id_num)
        sum(abs(sm_count - 1)) / n
      },
      update = {
        # Dick Brus (Jul 2018) proposes to compute O1 as the mean of the absolute deviations of marginal
        # stratum sample sizes. This should be the same as implemented in the FORTRAN code.
        mean(abs(sm_count - 1))
      })
    
    # Output
    # return (sum(abs(sm_count)) / n)
    return (out)
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE (O2) ######################################################
# sm: sample matrix
# n_pts: number of points
# id_fac: columns of sm containing factor covariates
# pop_prop: population class proportions (DEPRECATED)
# pop_count: population class counts
# n_candi: number of candidate locations (population)
# clhs.version: CLHS version
.objO2 <-
  # function (sm, id_fac, n_pts, pop_prop, clhs.version) {
  function (sm, id_fac, n_pts, pop_count, n_candi, clhs.version) {
    
    # Count the number of sample points per class
    sm_count <- lapply(sm[, id_fac], function(x) table(x))
    
    # Compute the sample proportions (DEPRECATED)
    # sm_prop <- lapply(sm[, id_fac], function(x) table(x) / n_pts)
    
    # Compare the sample and population proportions (DEPRECATED)
    # sm_prop <- sapply(1:length(id_fac), function (i)
      # sum(abs(sm_prop[[i]] - pop_prop[[i]])))
    
    out <- switch (clhs.version,
      paper = {
        # Minasny and McBratney (2006)
        sm_prop <- lapply(sm_count, function (x) x / n_pts)
        pop_prop <- lapply(pop_count, function (x) x / n_candi)
        sum(sapply(1:length(id_fac), function (i) sum(abs(sm_prop[[i]] - pop_prop[[i]]))))
      },
      fortran = {
        # Minasny and McBratney (2006)
        sm_prop <- lapply(sm_count, function (x) x / n_pts)
        pop_prop <- lapply(pop_count, function (x) x / n_candi)
        sum(sapply(1:length(id_fac), function (i) sum(abs(sm_prop[[i]] - pop_prop[[i]]))))
      }, 
      update = {
        # Dick Brus (Jul 2018) proposes to compute O2 as the mean of the absolute deviations of marginal
        # stratum sample sizes, defined just like O1 in terms of sample sizes. Defined in this alternative 
        # way O1 and O2 should be fully comparable.
        # mean(abs(n_realized - n_populational)
        mean(sapply(1:length(id_fac), function (i) abs(pop_count[[i]] - pop_count[[i]])))
      })
    
    # Output
    # return (sum(sm_prop))
    return(out)
  }
# INTERNAL FUNCTION - CALCULATE THE CRITERION VALUE (O3) ######################################################
# sm: sample matrix
# id_num: columns of sm containing numeric covariates
# pcm: population correlation matrix
# clhs.version: CLHS version
.objO3 <-
  function (sm, id_num, pcm, clhs.version) {
    
    # Calculate sample correlation matrix
    scm <- stats::cor(x = sm[, id_num], use = "complete.obs")
    
    out <- switch (clhs.version,
      paper = {
        # Minasny and McBratney (2006)
        sum(abs(pcm - scm))
      },
      fortran = {
        # The late FORTRAN code of Budiman Minasny -- ca. 2015 -- implements scaling factors so that values
        # are "more" comparable among objective functions. For O3, the scaling factor is defined as 
        # n * n / 2 + n. The rationale for this scaling factor is not clear.
        n <- length(id_num)
        n <- n * n / 2 + n
        sum(abs(pcm - scm)) / n
      },
      update = {
        # Dick Brus (Jul 2018) proposes to compute O3 as the mean of the off diagonal elements of the matrix
        # of absolute differences between sample and population correlation matrices. Defined in this 
        # alternative way, O3 should be fully comparable with O1 and O2.
        r_diff <- abs(pcm - scm)
        mean(r_diff[row(r_diff) != col(r_diff)])
      })
    
    # Output
    # return(sum(abs(pcm - scm)) / n)
    return(out)
  }
# INTERNAL FUNCTION - PREPARE OBJECT TO STORE THE BEST ENERGY STATE ###########################################
.bestEnergyCLHS <-
  function (covars_type) {
    if (covars_type == "both") {
      return (data.frame(obj = Inf, O1 = Inf, O2 = Inf, O3 = Inf))
    } else {
      if (covars_type == "numeric") {
        return (data.frame(obj = Inf, O1 = Inf, O3 = Inf))
      } else {
        return (data.frame(obj = Inf))
      }
    }
  }
