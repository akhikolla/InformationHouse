predict_clusterwise <-
function(object, newdata, ...) {
  if (!is(object, "disclapmixfit")) stop("object must be a disclapmixfit")
  probs <- rcpp_calculate_haplotype_probabilities_clusterwise(newdata, object$y, object$disclap_parameters, object$tau)
  return(probs)
}


theta_weir <- function(hapsum, ...) {
  if (!is(hapsum, "happrobsum_within_between")) stop("object must be a happrobsum_within_between")

  mi <- hapsum$match_within
  mij <- hapsum$match_between

  mW <- mean(mi)
  mA <- mean(mij[upper.tri(mij)])
  
  r <- length(mi)
  
  theta_approx <- (mW - mA) / (1 - mA)
  theta <- (((r-1)/r) * theta_approx) / (1 - (1/r)*theta_approx)
  
  thetas_approx <- (mi - mA) / (1 - mA)
  thetas <- (((r-1)/r) * thetas_approx) / (1 - (1/r)*thetas_approx)
  
  #return(weir_theta)
  return(list(theta = theta, theta_approx = theta_approx, thetas_subpops = thetas))
}


#  * New helper function: convert_to_compact_db
convert_to_compact_db <- function(x) { 
  #if (!is.matrix(x) || !is.data.frame(x)) {
  #  stop("x must be a matrix or data.frame")
  #}

  if (nrow(x) <= 1L) {
    ret <- data.frame(x, Ndb = 1L)
    ret$ind <- list("1" = 1L)
    return(ret)
  }

  #if (is.matrix(x) && (dim(x)[2L] == 1L)) {
  #  x <- as.vector(x) 
  #}
 
  x_ord <- do.call(order, as.data.frame(x))
   
  #if (is.vector(x)) {
  #  same_as_previous <- x[tail(x_ord, -1L)] == x[head(x_ord, -1L)]
  #} else {
    same_as_previous <- rowSums(x[tail(x_ord, -1L), , drop = FALSE] != x[head(x_ord, -1L), , drop = FALSE]) == 0L
  #}	
 
  indices <- split(x_ord, cumsum(c(TRUE, !same_as_previous)))
 
  #if (is.vector(x)) {
  #  x <- x[sapply(indices, function (x) x[[1L]]), drop = FALSE]
  #} else {
    x <- x[sapply(indices, function (x) x[[1L]]), , drop = FALSE]
  #}
  
  return(data.frame(x, Ndb = sapply (indices, length), ind = I(indices)))
}

#  * New helper function: find_haplotype_in_matrix
find_haplotype_in_matrix <- function(mat, haplotype) {
  if (!is.matrix(mat) || !is.integer(mat) || !is.integer(haplotype)) {
    stop("mat must be an integer matrix and haplotype an integer vector")
  }
  if (length(haplotype) != ncol(mat)) stop("Wrong dimensions")

  i <- rcpp_find_haplotype_in_matrix(mat, haplotype)

  if (i <= 0L) {
    i <- NULL
  }

  return(i)
}





#' Calculate haplotype diversity from a disclapmixfit
#' 
#' Calculate haplotype diversity from a \code{\link{disclapmixfit}} object. The
#' method is based on simulating a huge database that approximates the
#' population.
#' 
#' 
#' @param object a \code{\link{disclapmixfit}} object, usually from a result of
#' a call to \code{disclapmix}.
#' @param nsim number of haplotypes to generate for calculating the haplotype
#' diversity.
#' @return The calculated haplotype diversity.
#' @seealso \code{\link{disclapmix}} \code{\link{disclapmixfit}}
#' \code{\link{predict.disclapmixfit}} \code{\link{print.disclapmixfit}}
#' \code{\link{summary.disclapmixfit}} \code{\link{simulate.disclapmixfit}}
#' %\code{\link{haplotype_diversity}} \code{\link{clusterdist}}
#' @keywords print
haplotype_diversity <- function(object, nsim = 1e4L) {
  if (!is(object, "disclapmixfit")) stop("object must be a disclapmixfit")
  
  if (is.null(nsim) || length(nsim) != 1L || !is.integer(nsim) || nsim <= 0L) {
    stop("nsim must be >= 1L (note the L postfix for integer)")
  }
  
  get_db_counts <- function(x) { 
    order_x <- do.call(order, as.data.frame(x))
    equal.to.previous <- rowSums(x[tail(order_x, -1),] != x[head(order_x, -1),]) == 0
    indices <- split(order_x, cumsum(c(TRUE, !equal.to.previous)))
    Ns <- unlist(lapply(indices, length))
    return(Ns)
  }
  
  db <- simulate.disclapmixfit(object, nsim = nsim)
  Ns <- get_db_counts(db)
  freqs <- Ns / nsim
  #D <- 1 - sum(freqs^2)
  D <- (nsim / (nsim - 1)) * (1 - sum(freqs^2))
  
  return(D)
}





predict_with_variance <- function(fit, newdata, nsim = 1000L) {
  stopifnot(nsim >= 1L) 
  
  #new_betas <- mvtnorm::rmvnorm(nsim, fit$glm_coef, fit$covmat)
  new_betas <- MASS::mvrnorm(nsim, fit$glm_coef, fit$covmat)
  p_sim <- matrix(NA, nrow = nrow(newdata), ncol = nsim)
  
  for (iter in 1L:nsim) {
    new_beta <- new_betas[iter, , drop = FALSE]
    
    new_discps <- convert_coef_to_disclap_parameters_internal(new_beta, nrow(fit$y))
    
    new_wic <- rcpp_calculate_wic(fit$x, fit$y, new_discps, fit$tau)
    new_vic_matrix <- rcpp_calculate_vic(new_wic)
    
    new_tau_vector <- apply(new_vic_matrix, 2, sum) / nrow(fit$x)
    
    new_ys <- move_centers(fit$x, fit$y, new_vic_matrix)
    
    try({
      ps_new <- rcpp_calculate_haplotype_probabilities(newdata, new_ys, new_discps, new_tau_vector)
      p_sim[, iter] <- ps_new
    })
  }
  
  p_sim_mean <- apply(p_sim, 1, mean, na.rm = TRUE)
  p_sim_sd <- apply(p_sim, 1, sd, na.rm = TRUE)
  

  #p_sim_qs <- t(apply(p_sim, 1, quantile, c(0.005, 0.025, 0.05, 0.95, 0.975, 0.995)))
  #colnames(p_sim_qs) <- paste0("q", gsub("%", "perc", fixed = TRUE, colnames(p_sim_qs)))

  #return(data.frame(mean = p_sim_mean, sd = p_sim_sd, p_sim_qs))
  
  return(data.frame(mean = p_sim_mean, sd = p_sim_sd))
}




predict_with_variance_rnd_ys <- function(fit, newdata, nsim = 1000L) {
  stopifnot(nsim >= 1L) 
  
  #new_betas <- mvtnorm::rmvnorm(nsim, fit$glm_coef, fit$covmat)
  new_betas <- MASS::mvrnorm(nsim, fit$glm_coef, fit$covmat)
  p_sim <- matrix(NA, nrow = nrow(newdata), ncol = nsim)
  
  for (iter in 1L:nsim) {
    new_beta <- new_betas[iter, , drop = FALSE]
    
    new_discps <- convert_coef_to_disclap_parameters_internal(new_beta, nrow(fit$y))
    
    # new
    new_ys <- fit$y
    for (j in 1L:nrow(fit$y)) {
      for (k in 1L:ncol(fit$x)) {
        new_ys[j, k] <- disclap::rdisclap(1, new_discps[j, k]) + fit$y[j, k]
      }
    }
    
    new_wic <- rcpp_calculate_wic(fit$x, new_ys, new_discps, fit$tau)
    new_vic_matrix <- rcpp_calculate_vic(new_wic)
    
    new_tau_vector <- apply(new_vic_matrix, 2, sum) / nrow(fit$x)
    
    #new_ys <- move_centers(fit$x, fit$y, new_vic_matrix)
    
    try({
      ps_new <- rcpp_calculate_haplotype_probabilities(newdata, new_ys, new_discps, new_tau_vector)
      p_sim[, iter] <- ps_new
    })
  }
  
  p_sim_mean <- apply(p_sim, 1, mean, na.rm = TRUE)
  p_sim_sd <- apply(p_sim, 1, sd, na.rm = TRUE)
  

  #p_sim_qs <- t(apply(p_sim, 1, quantile, c(0.005, 0.025, 0.05, 0.95, 0.975, 0.995)))
  #colnames(p_sim_qs) <- paste0("q", gsub("%", "perc", fixed = TRUE, colnames(p_sim_qs)))

  #return(data.frame(mean = p_sim_mean, sd = p_sim_sd, p_sim_qs))
  
  return(data.frame(mean = p_sim_mean, sd = p_sim_sd))
}



