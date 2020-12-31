#################################################
## This function computes the complementary cdf of one sample two-sided
## Kolmogorov-Smirnov statistic when the distribution under the
## null hypothesis is purely discrete

disc_ks_c_cdf <- function(q, n, y, ..., exact = NULL, tol = 1e-8, sim.size = 1000000, num.sim = 10)
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

  PVAL <- NULL

  ######################################################################
  if (n <= 0){
    stop("'n' must be larger than zero")
    #PVAL <- NULL
  }
  else{
    if (q <= 0){
      warning("For any 'q' smaller than or equal to 0, the complementary cdf is 1")
      PVAL <- 1
    }
    else if (q >= 1){
      warning("For any 'q' greater than or equal to 1, the complementary cdf is 0")
      PVAL <- 0
    }
    else{
      if (is.stepfun(y)){
        z <- knots(y)

        if(is.null(exact)){
          exact <- (n <= 100000)
        }

        if (exact){# Exact-KS-FFT method

          upper_rect <- upper_rectangles_1(q, n, y, z, tol)
          lower_rect <- lower_rectangles_1(q, n, y, z, tol)

          df <- data.frame(rbind(upper_rect, lower_rect))
          write.table(df,"Boundary_Crossing_Time.txt", sep = ", ", row.names = FALSE, col.names = FALSE)

          PVAL <- KSgeneral::ks_c_cdf_Rcpp(n)

          file.remove("Boundary_Crossing_Time.txt")
        }
        else {# Wood and Altavela (1978)'s MC-based method
          if ((n > 100000) && (z > 15)){
            warning("For large sample sizes, use Wood and Altavela (1978)'s MC-based, approximated asymptotic results")
            PVAL <- WA_Multi_No_Correction(sim.size, n, (q*n^(1/2)), y, num.sim)
          }
          else if ((n > 100000) && (z <= 15)){
            warning("For large sample sizes, use Wood and Altavela (1978)'s MC-based, approximated asymptotic results")
            PVAL <- WA_Multi_Correction(sim.size, n, (q*n^(1/2)), y, num.sim)
          }
          else if ((n <= 100000) && (z > 15)){
            warning("For small sample sizes, use of asymptotic results may not be efficient or accurate enough, especially when the sample size is smaller than 1000. Recommend to put 'exact = TRUE'")
            PVAL <- WA_Multi_No_Correction(sim.size, n, (q*n^(1/2)), y, num.sim)
          }
          else if ((n <= 100000) && (z <= 15)){
            warning("For small sample sizes, use of asymptotic results may not be efficient or accurate enough, especially when the sample size is smaller than 1000. Recommend to put 'exact = TRUE'")
            PVAL <- WA_Multi_Correction(sim.size, n, (q*n^(1/2)), y, num.sim)
          }
        }

      }
      else{
        stop("'y' must be a step function")
      }
    }
  }

  PVAL <- min(1.0, max(0.0, PVAL))

  # file.remove("Boundary_Crossing_Time.txt")

  return(PVAL)

}
