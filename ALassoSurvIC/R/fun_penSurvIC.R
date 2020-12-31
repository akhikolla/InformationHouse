fun_penSurvIC <- function(theta, tilde_b, arglist, tol) {

  n <- arglist$n
  npmle_set <- arglist$set
  u <- npmle_set[,2]
  if(missing(tol)) tol <- arglist$tol
  niter <- arglist$niter

  old_b <- tilde_b
  old_lambda <- arglist$initial_lambda

  all_b_info <- old_b
  all_lambda_info <- old_lambda
  all_distance <- NA

  distance <- tol + 1000
  iter <- 0

  while ((distance > tol) & (iter < niter)) {
    ew <- fun_ew(old_b, old_lambda, arglist)
    neg_gradient <- fun_dq(old_b, ew, arglist) * (-1)
    neg_hess <- fun_ddq(old_b, ew, arglist) * (-1)
    x <- chol(neg_hess)
    y <- solve(t(x)) %*% (as.numeric(neg_hess %*% old_b) - neg_gradient)

    new_b <- fun_shooting_algorithm(x, y, theta, tilde_b, arglist)
    new_lambda <- fun_updatelambda(new_b, ew, arglist)
    distance <- max(abs(c(new_b, new_lambda) - c(old_b, old_lambda)))

    all_b_info <- rbind(all_b_info, new_b)
    all_lambda_info <- rbind(all_lambda_info, new_lambda)
    all_distance <- c(all_distance, distance)

    old_b <- new_b
    old_lambda <- new_lambda
    iter <- iter + 1
  }

  results <- list()
  results$b <- new_b
  results$lambda <- new_lambda
  results$theta <- theta
  results$iteration <- iter
  results$convergence <- ifelse(iter < niter, TRUE, FALSE)
  results$distance <- distance
  results$tol <- tol
  results$niter <- niter

  return(results)

}
