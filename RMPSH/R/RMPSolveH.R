#' Recursive Modified Direct Search on Hyper-rectangle
#'
#' `RMPSolveH` can be Used to Minimize any Non-Convex Blackbox Function where Each Parameter
#' has an Upper Bound and Lower Bound.
#'
#'
#' @param x0 Vector of Initial Guess provided by User.
#' @param func The Function to be Optimized, should be provided by the User.
#' @param lb Vector of Lower Bounds, of same Dimension as 'x0'.
#' @param ub Vector of Upper Bound, of same Dimension as 'x0'
#' @param rho_1 'Step Decay Rate' for the First Run Only (Default is 2).
#' @param rho_2 'Step Decay Rate' for Second Run Onwards (Default is 2).
#' @param phi Lower Bound for 'Global Step Size'. Default value is \eqn{10^{-6}}.
#' @param no_runs Max Number of 'Runs'. Default Value is 1000.
#' @param max_iter Max Number of Iterations in each 'Run'. Default Value is 10000.
#' @param s_init Initial  'Global Step Size'. Default Value is 2. It must be set Less than or Equal to 2.
#' @param tol_fun Termination Tolerance on when to decrease the 'Global Step Size'. Default Value is \eqn{10^{-6}}. For more accuracy, user may set it to a Smaller Value
#' e.g., \eqn{10^{-20}}. However, for Expensive Objective Functions, for Faster Computation, User should set it to a Larger Value e.g, \eqn{10^{-3}}.
#' @param tol_fun_2 Termination Tolerance on the Difference of Norms of solution points in two Consecutive Runs. Default Value is \eqn{10^{-20}}.
#' However, for Expensive Objective Functions, for Faster Computation, user should set it to a Larger Value e.g, \eqn{10^{-6}}.
#' @param max_time Time Alloted (In Seconds) for Execution of RMPSH. Default is 36000 secs (10 Hours).
#' @param print_output Binary Command to Print Optimized Value of Objective Function after Each Iteration. Default is set as FALSE.
#'
#' @return The Optimal Solution Point.
#'
#' @examples
#'
#' g <- function(y)
#'  return(-20 * exp(-0.2 * sqrt(0.5 * (y[1] ^ 2 + y[2] ^ 2))) -
#'  exp(0.5 * (cos(2 * pi * y[1]) + cos(2 * pi * y[2]))) + exp(1) + 20)
#'
#' starting_point <- rep(1, 10)
#'
#' g(starting_point)
#'
#' solution <- RMPSolveH(starting_point, g, rep(-33, 10), rep(33, 10))
#'
#' g(solution)
#'
#' RMPSolveH(c(2, 4, 6, 2, 1), g, rep(-3, 5), rep(23, 5), print_output = TRUE)
#' # Will Print the Updates after Each Iteration
#'
#'
#' g <- function(y)
#'  return(sum(y ^ 2))
#' RMPSolveH(rep(2.3, 100),
#'           g,
#'           rep(-11, 100),
#'           rep(13, 100),
#'           max_time = 2,
#'           print = 1)
#' # Will Exit and Return Result after 2 Seconds
#'
#' @references
#' \itemize{
#'
#'   \item Das, Priyam \cr
#'    "Black-box optimization on hyper-rectangle using Recursive Modified Pattern Search and application to ROC-based Classification Problem" \cr
#'          (available at `arXiv \url{http://arxiv.org/abs/1604.08616}).
#' }
#'
#'@name RMPSolveH
NULL

Rcpp::sourceCpp('src/anti_transformation.cpp')
Rcpp::sourceCpp('src/transformation.cpp')
Rcpp::sourceCpp('src/update.cpp')


#' @rdname RMPSolveH
#' @export

RMPSolveH <-
  function(x0,
           func,
           lb,
           ub,
           rho_1 = 2,
           rho_2 = 2,
           phi = 1e-6,
           no_runs = 1e3,
           max_iter = 1e4,
           s_init = 2,
           tol_fun = 1e-6,
           tol_fun_2 = 1e-20,
           max_time = 36e3,
           print_output = FALSE)
  {
    M <- length(x0)
    start_value <- func(x0)
    start_time <- Sys.time()

    theta_array <- matrix(0, no_runs, M)
    each_run_solution <- array(0, no_runs)

    ill_condition <- 0
    for (iii in 1:no_runs)
    {
      epsilon <- s_init
      if (iii == 1)
      {
        rho <- rho_1
        theta <- anti_transformation(x0, lb, ub)
        if (min(ub - lb) < 0)
        {
          print("upper bound should be greater than lower bound")
          break
        }
        if (max(abs(theta)) > 1 ||
            max(abs(theta)) < 0)
        {
          print("Starting point is outside the domain")
          break
        }
      }
      else
      {
        theta <- as.array(theta_array[iii - 1,])
        rho <- rho_2
      }
      if (ill_condition == 1)
        break

      array_of_values <- array(0, max_iter)
      for (i in 1:max_iter)
      {
        current_lh <- func(transformation(theta, lb, ub))
        possible_x_coords <-
          update(theta, epsilon, rho, phi)
        total_lh <- rep(1, 2 * M)

        for (kk in 1:(2 * M))
        {
          candidate_theta <- theta
          coord_number <- ceiling(kk / 2)

          if (possible_x_coords[kk] == theta[coord_number])
            total_lh[kk] <- current_lh
          else
          {
            candidate_theta[coord_number] <- possible_x_coords[kk]
            total_lh[kk] <-
              func(transformation(candidate_theta, lb, ub))
          }
        }
        new_min <- min(total_lh)

        if (new_min < current_lh)
        {
          position_new_min <- which.min(total_lh)
          pos_of_theta <- ceiling(position_new_min / 2)
          theta[pos_of_theta] <- possible_x_coords[position_new_min]
        }

        array_of_values[i] <- min(new_min, current_lh)

        if (print_output)
        {
          cat(
            'Run no. ',
            iii,
            ',iteration no.',
            i,
            ', current fun value = ',
            array_of_values[i],
            '\n'
          )
        }

        if (i > 1)
          if (abs(array_of_values[i] - array_of_values[i - 1]) < tol_fun)
          {
            if (abs(epsilon) > phi)
              epsilon <- epsilon / rho
            else
              break
          }

        now_time <- Sys.time()
        now_time_spent <- as.numeric(now_time - start_time)

        if (now_time_spent > max_time)
        {
          cat("\n \n \n")
          paste0("Starting objective function value:  ", start_value)
          paste0("Final objective function value:  ", current_lh)
          ill_condition <- 1
          break
        }

      }

      theta_array[iii,] <- t(theta)
      each_run_solution[iii] <- func(transformation(theta, lb, ub))

      if (iii > 1)
      {
        old_soln <- theta_array[iii - 1,]
        new_soln <- theta_array[iii,]

        if (norm(as.matrix(new_soln - old_soln)) < tol_fun_2)
          break
      }
    }

    end_time <- Sys.time()
    time_spent <- as.numeric(end_time - start_time)
    final_value <- func(transformation(theta, lb, ub))

    paste0("Starting objective function value:  ", start_value)
    paste0("Final objective function value:  ", final_value)
    paste0("Total time required (in secs):  ", time_spent)
    paste0("Obtained minima point is :  ")

    transformation(theta, lb, ub)

    return(transformation(theta, lb, ub))

  }
