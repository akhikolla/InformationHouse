#################################################
#################################################
## This function computes the complementary cdf of one sample two-sided
## Kolmogorov-Smirnov statistic when the distribution under the
## null hypothesis is mixed

mixed_ks_c_cdf <- function(q, n, jump_points, Mixed_dist, ..., tol = 1e-10)
{

  options("digits" = 16)


  if (n <= 0){
    stop("'n' must be larger than zero")
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
      Vec_Mixed_cdf <- function(x)
      {
        return (sapply(x, Mixed_dist))
      }

      a <- rep(NA,n)
      f_a <- a
      b <- a
      f_b <- a

      temp_a <- a
      temp_b <- a

      temp_length <- length(jump_points) - 1

      temp <- rep(0, (2*temp_length))

      #################################################
      #################################################

      for (i in 1:temp_length){

        temp[2*i-1] <- Vec_Mixed_cdf(jump_points[i])
        temp[2*i] <- Vec_Mixed_cdf(jump_points[i+1] - (1e-13))

        i <- i + 1
      }

      for (i in 1:n){

        # Rectangle for the uniform order statistics obtained from Gleser(1985)
        # The rectangles are found in a similar way as in the function ks.test
        # in package dgof by Arnold and Emerson (2011)

        a[i] <- min(c(jump_points[which(Vec_Mixed_cdf(jump_points) + q >= i/n + tol)[1]], Inf), na.rm = TRUE)
        b[i] <- min(c(jump_points[which(Vec_Mixed_cdf(jump_points) - q > (i-1)/n - tol)[1]], Inf), na.rm = TRUE)

        f_a[i] <- ifelse(!(a[i] %in% jump_points), Vec_Mixed_cdf(a[i]), Vec_Mixed_cdf(a[i] - tol))

        f_b[i] <- Vec_Mixed_cdf(b[i])

        temp_a[i] <- max(i/n - q + tol, 0)
        temp_b[i] <- min((i-1)/n + q - tol, 1)

        for (j in 1:temp_length){

          if (((temp_a[i] > temp[2*j-1]) && (temp_a[i] <= temp[2*j])) || ((temp_a[i] <= Vec_Mixed_cdf(jump_points[1] - (1e-13)))) || ((temp_a[i] > Vec_Mixed_cdf(jump_points[length(jump_points)])))){

            f_a[i] <- i/n - q
          }

          if (((temp_b[i] >= temp[2*j-1]) && (temp_b[i] < temp[2*j])) || ((temp_b[i] < Vec_Mixed_cdf(jump_points[1] - (1e-13)))) || ((temp_b[i] >= Vec_Mixed_cdf(jump_points[length(jump_points)])))){

            f_b[i] <- (i-1)/n + q
          }

          j <- j + 1
        }

        if (f_a[i] < 0){
          f_a[i] <- 0
        }
        if (f_b[i] > 1){
          f_b[i] <- 1
        }

        i <- i + 1

      }

      df <- data.frame(rbind(f_b, f_a))
      write.table(df,"Boundary_Crossing_Time.txt", sep = ", ", row.names = FALSE, col.names = FALSE)

      PVAL <- KSgeneral::ks_c_cdf_Rcpp(n)

      file.remove("Boundary_Crossing_Time.txt")

    }
  }

  PVAL <- min(1.0, max(0.0, PVAL))

  return(PVAL)

}
