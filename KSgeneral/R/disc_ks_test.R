#################################################
## This function computes the p-value of one sample two-sided
## Kolmogorov-Smirnov test when the distribution under the
## null hypothesis is purely discrete

disc_ks_test <- function(x, y, ..., exact = NULL, tol = 1e-8, sim.size = 1000000, num.sim = 10)
{
  upper_rectangles_1 <- function(S, n, y, knots.y, tol)
  {
    # Rectangle for the uniform order statistics obtained from Gleser(1985)
    # The rectangles are found in a similar way as in the function ks.test
    # in package dgof by Arnold and Emerson (2011)

    eps <- min(tol, min(diff(knots.y)) * tol)
    eps2 <- min(tol, min(diff(y(knots.y))) * tol)

    a <- rep(0, n)
    b <- a
    f_a <- a


    for (i in 1:n){

      b[i] <- min(c(knots.y[which(y(knots.y) - S > (i-1)/n - eps2)[1]], Inf), na.rm = TRUE)

    }

    f_b <- y(b)

    return (f_b)
  }

  #################################################################################

  lower_rectangles_1 <- function(S, n, y, knots.y, tol)
  {
    # Rectangle for the uniform order statistics obtained from Gleser(1985)
    # The rectangles are found in a similar way as in the function ks.test
    # in package dgof by Arnold and Emerson (2011)

    eps <- min(tol, min(diff(knots.y)) * tol)
    eps2 <- min(tol, min(diff(y(knots.y))) * tol)

    a <- rep(0, n)
    b <- a
    f_a <- a

    for (i in 1:n){

      a[i] <- min(c(knots.y[which(y(knots.y) + S >= i/n + eps2)[1]], Inf), na.rm = TRUE)
      f_a[i] <- ifelse(!(a[i] %in% knots.y), y(a[i]), y(a[i] - eps))

    }

    return(f_a)
  }

  ######################################################################
  ######################################################################
  ## Wood and Altavela (1978) simulated approach to approximating the
  ## asymptotic distribution of the KS statistic when the underlying
  ## distribution is purely discrete

  WA_Single <- function(size, n, lambda, y)
  {

    z <- knots(y)

    lth <- length(z) - 1
    pmf_ <- y(z)[1 : lth]

    a <- rep(0, (lth)*(lth))

    for (i in 1 : lth){

      for (j in i : lth){

        if (j == i){
          a[lth * i - lth + j] <- pmf_[i] * (1 - pmf_[j]) / 2
        }
        else{
          a[lth * i - lth + j] <- pmf_[i] * (1 - pmf_[j])
        }

        j <- j + 1
      }

      i <- i + 1
    }

    matrix_1 <- matrix(a, nrow = lth, ncol = lth, byrow = TRUE)
    matrix_2 <- t(matrix_1)

    Cov_matrix <- matrix_1 + matrix_2

    RV <- MASS::mvrnorm(size, rep(0, lth), Cov_matrix)

    counting <- 0
    for(i in 1 : size){

      storing <- RV[i, ]

      indicator <- 1
      for (j in 1: lth){
        indicator <- indicator * (abs(storing[j]) < lambda)
        j <- j + 1
      }
      if (indicator){
        counting <- counting + 1
      }

      i <- i + 1
    }
    prob <- counting / size
    options("digits" = 8)

    return (1 - prob)
  }
  ######################################################################
  ## Wood and Altavela (1978) approach without continuity correction
  WA_Multi_No_Correction <- function(size, n, lambda, y, reps){

    result <- 0
    for (i in 1 : reps){
      set.seed(i)
      result <- result + WA_Single(size, n, lambda, y)

      i <- i + 1
    }

    return(result/reps)
  }
  ######################################################################
  ## Wood and Altavela (1978) approach with continuity correction
  WA_Multi_Correction <- function(size, n, lambda, y, reps){

    result <- 0
    for (i in 1 : reps){
      set.seed(i)
      result <- result + WA_Single(size, n, (lambda - 0.5*(n)^(-1/2)), y)


      i <- i + 1
    }
    return(result/reps)

  }

  ######################################################################

  DNAME <- deparse(substitute(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 1L){
    stop("not enough 'x' data")
  }
  PVAL <- NULL

  #########################################################

  if (is.stepfun(y)){
    z <- knots(y)

    METHOD <- "One-sample Kolmogorov-Smirnov test"
    dev <- c(0, ecdf(x)(z) - y(z))
    STATISTIC <- max(abs(dev))

    #####
    if(is.null(exact)){
      warning("Recommend use of 'exact = TRUE', since for q*n^(1/2) large, the Wood and Altavela (1978)'s MC-based, approximated asymptotic results may not be accurate")
      exact <- (n <= 100000)
    }

    if (exact){

    #####

      upper_rect <- upper_rectangles_1(STATISTIC, n, y, z, tol)
      lower_rect <- lower_rectangles_1(STATISTIC, n, y, z, tol)

      df <- data.frame(rbind(upper_rect, lower_rect))
      write.table(df,"Boundary_Crossing_Time.txt", sep = ", ", row.names = FALSE, col.names = FALSE)

      PVAL <- KSgeneral::ks_c_cdf_Rcpp(n)

      file.remove("Boundary_Crossing_Time.txt")
    }
    else {
      if ((n > 100000) && (z > 15)){
        warning("For large sample sizes, use Wood and Altavela (1978)'s MC-based, approximated asymptotic results")
        PVAL <- WA_Multi_No_Correction(sim.size, n, (STATISTIC*n^(1/2)), y, num.sim)
      }
      else if ((n > 100000) && (z <= 15)){
        warning("For large sample sizes, use Wood and Altavela (1978)'s MC-based, approximated asymptotic results")
        PVAL <- WA_Multi_Correction(sim.size, n, (STATISTIC*n^(1/2)), y, num.sim)
      }
      else if ((n <= 100000) && (z > 15)){
        warning("For small sample sizes, use of asymptotic results may not be efficient or accurate enough, especially when the sample size is smaller than 1000. Recommend to put 'exact = TRUE'")
        PVAL <- WA_Multi_No_Correction(sim.size, n, (STATISTIC*n^(1/2)), y, num.sim)
      }
      else if ((n <= 100000) && (z <= 15)){
        warning("For small sample sizes, use of asymptotic results may not be efficient or accurate enough, especially when the sample size is smaller than 1000. Recommend to put 'exact = TRUE'")
        PVAL <- WA_Multi_Correction(sim.size, n, (STATISTIC*n^(1/2)), y, num.sim)
      }
    }

    nm_alternative <- "two-sided"

  }
  else{
    if (is.character(y)){
      y <- get(y, mode = "function", envir = parent.frame())
    }
    if (!is.function(y)){
      stop("'y' must be a function or a string naming a valid function")
    }
    METHOD <- "One-sample Kolmogorov-Smirnov test"
    TIES <- FALSE

    if (length(unique(x)) < n){
      stop("ties should not be present for the continuous Kolmogorov-Smirnov test")
    }

    x <- y(sort(x), ...) - (0 : (n - 1))/n
    STATISTIC <- max(c(x, 1/n - x))

    PVAL <- KSgeneral::cont_ks_c_cdf(STATISTIC, n)

  # file.remove("Boundary_Crossing_Time.txt")

    nm_alternative <- "two-sided"
  }

  names(STATISTIC) <- "D"

  PVAL <- min(1.0, max(0.0, PVAL))
  RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative, method = METHOD, data.name = DNAME)

  class(RVAL) <- "htest"

  # file.remove("Boundary_Crossing_Time.txt")

  return(RVAL)

}
