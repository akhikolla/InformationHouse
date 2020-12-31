#################################################
## This function computes the cdf of one sample two-sided Kolmogorov-Smirnov
## statistic when the distribution under the null hypothesis is continuous

cont_ks_cdf <- function(q, n)
{

#  options("digits" = 20)

  PVAL <- NULL

  ################################################################
  if (n <= 0){
    stop("'n' must be larger than zero")
  }
  else{
    if (q <= 0){
      warning("For any 'q' smaller than or equal to 0, the cdf is 0")
      PVAL <- 0
    }
    else if (q >= 1){
      warning("For any 'q' greater than or equal to 1, the cdf is 1")
      PVAL <- 1
    }
    else{

      #################################
      ## Please see Dimitrova et al. (2017) for more details

      if ((n <= 140) && (n*(q)^2 > 12)){
        PVAL <- 1
      }
      else if ((n > 140) && (n <= 100000) && (n*(q)^(3/2) >= 1.4) && (n*(q)^2 >= 11)){
        PVAL <- 1
      }
      else if (n*(q)^2 >= 18){
        PVAL <- 1
      }
      else{

        #################################
        f_a <- rep(NA, n)
        f_b <- rep(NA, n)

        for (i in 1:n){

          f_b[i] <- (i-1)/n + q
          f_a[i] <- (i+0)/n - q

          if (f_b[i] > 1){
            f_b[i] <- 1
          }
          if (f_a[i] < 0){
            f_a[i] <- 0
          }

          i <- i + 1
        }

        df <- data.frame(rbind(f_b, f_a))
        write.table(df,"Boundary_Crossing_Time.txt", sep = ", ", row.names = FALSE, col.names = FALSE)
        #################################

        PVAL <- 1 - KSgeneral::ks_c_cdf_Rcpp(n)

        file.remove("Boundary_Crossing_Time.txt")

      }

    }
  }

  PVAL <- min(1.0, max(0.0, PVAL))

  return(PVAL)

}
####################################################
####################################################
## This function computes the complementary cdf of one sample two-sided
## Kolmogorov-Smirnov statistic when the distribution under the null
## hypothesis is continuous

cont_ks_c_cdf <- function(q, n)
{

  PVAL <- 1 - cont_ks_cdf(q, n)

  return(PVAL)

}
