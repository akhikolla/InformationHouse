fun_unpenSurvIC <- function(b, arglist) {

  n <- arglist$n
  npmle_set <- arglist$set
  tol <- arglist$tol
  niter <- arglist$niter

  old_b <- b
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

    new_b <- as.numeric(solve(t(x)%*%x)%*%t(x)%*%y)
    new_lambda <- fun_updatelambda(new_b, ew, arglist)

    distance <- max(abs(c(new_b, new_lambda) - c(old_b, old_lambda)))
    all_b_info <- rbind(all_b_info, new_b)
    all_lambda_info <- rbind(all_lambda_info, new_lambda)
    all_distance <- c(all_distance, distance)

    old_b <- new_b
    old_lambda <- new_lambda
    iter <- iter + 1
  }

  results <- NULL

  results <- list()
  results$b <- new_b
  results$lambda <- new_lambda
  results$iteration <- iter
  results$convergence <- ifelse(iter < niter, TRUE, FALSE)
  results$distance <- distance
  results$tol <- tol
  results$niter <- niter

  return(results)
}
