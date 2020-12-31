#' Simulate and Analyze LSBCLUST
#' 
#' Perform a single simulation run for the LSBCLUST model. Multiple data sets 
#' are generated for a single set of underlying parameters, 
#' 
#' @inheritParams rlsbclust
#' @param seed An optional seed to be set for the random number generator
#' @param nstart_T3 The number of random starts to use for \code{\link{T3Clusf}}
#' @param nstart_ak The number of random starts to use for \code{\link{akmeans}}
#' @param verbose Integer giving the number of iterations after which the loss values is printed.
#' @param parallel_data Logical indicating whether to parallelize over the data sets. If 
#' \code{FALSE}, parallelization is done over random starts (depending on \code{parallel}).
#' @param parallel Logical indicating whether to parallelize over random starts. 
#' Note that \code{parallel_data} has precedence over this
#' @param mc.cores The number of cores to use, passed to \code{\link{makeCluster}}
#' @param include_fits Logical indicating whether to include the model fits, or
#' or only the fit statistics
#' @param include_data Logical indicating whether to include the simulated data 
#' fitted on, or only the results
#' @param nstart From \code{\link{lsbclust}}
#' @param nstart.kmeans From \code{\link{lsbclust}}
# @param \dots Additional arguments passed to \code{\link{lsbclust}}
#' @export
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @examples 
#' set.seed(1)
#' res <- sim_lsbclust(ndata = 5, nobs = 100, size = c(10, 8), nclust = rep(5, 4), 
#'                     verbose = 0, nstart_T3 = 2, nstart_ak = 1, parallel_data = FALSE,
#'                     nstart = 2, nstart.kmeans = 5 )
#' 
sim_lsbclust <- function(ndata, nobs, size, nclust, clustsize = NULL, 
                         delta = rep(1L, 4L), ndim = 2L, alpha = 0.5, 
                         fixed = c("none", "rows", "columns"), 
                         err_sd = 1, svmins = 0.5, svmax = 5, 
                         seed = NULL, parallel = FALSE, parallel_data = TRUE,
                         verbose = 0, nstart_T3 = 20L, nstart_ak = 20L, 
                         mc.cores = detectCores() - 1, include_fits = FALSE, 
                         include_data = FALSE, ## Following are specific to lsbclust()
                         nstart, nstart.kmeans) {
  
  ## Start timing
  t0 <- proc.time()[3]
  
  ## Catch call
  cll <- match.call()
  
  ## Set seed if needed
  if (!is.null(seed))
    set.seed(seed)
  
  ## Expand nclust if needed
  if (length(nclust) == 1L)
    nclust <- rep(nclust, 4L)
  
  ## Expand svmins if needed
  if (length(svmins) == 1)
    svmins <- rep(svmins, ndim)
  
  ## Generate the data
  dat <- rlsbclust(ndata = ndata, nobs = nobs, size = size, nclust = nclust, 
                   clustsize = clustsize, delta = delta, ndim = ndim, 
                   alpha = alpha, fixed = fixed, err_sd = err_sd, 
                   svmins = svmins, svmax = svmax)
  
  ## Parallelize over data sets
  if (parallel_data) {
    
    ## Function to do all four analyses
    analyze <- function(dat, delta, nclust, 
                        ndim, alpha, fixed, 
                        parallel, verbose, nstart, 
                        nstart.kmeans) {
      
      ## Fit and evaluate lsbclust
      res_lsb <- lsbclust(data = dat$data, delta = delta, nclust = nclust, 
                               ndim = ndim, alpha = alpha, fixed = fixed, 
                               parallel = FALSE, verbose = verbose,  
                          nstart = nstart, nstart.kmeans = nstart.kmeans)
      eval_lsb <- cfsim(fitted = res_lsb, actual = dat)
      
      ## Double-centre array for T3Clusf
      dat$data <- carray(dat$data)
      
      ## Fit and evaluate T3clusf
      res_t3 <- T3Clusf(X = dat$data, Q = ndim, G = nclust[4L], 
                             parallel = FALSE, nstart = nstart_T3, verbose = verbose)
      eval_t3 <- cfsim(fitted = res_t3, actual = dat)
      
      ## Fit and evaluate akmeans
      res_ak <- akmeans(data = dat$data, centers = nclust[4L], ndim = ndim, 
                             nstart = nstart_ak)
      eval_ak <- cfsim(fitted = res_ak, actual = dat)
      
      ## Fit and evaluate akmeans_nodim
      res_ak_nodim <- akmeans(data = dat$data, centers = nclust[4L], ndim = NULL, 
                                   nstart = nstart_ak)
      eval_ak_nodim <- cfsim(fitted = res_ak_nodim, actual = dat)
      
      ## Return the results
      return(list(lsbclust_overall = eval_lsb$overall$stats, 
                  lsbclust_rows = eval_lsb$rows$stats,
                  lsbclust_columns = eval_lsb$columns$stats,
                  lsbclust_interactions = eval_lsb$interactions$stats, 
                  t3clus = eval_t3$interactions$stats, 
                  akmeans = eval_ak$interactions$stats,
                  akmeans_nodim = eval_ak_nodim$interactions$stats,
                  fit_lsbclust = res_lsb, fit_t3clus = res_t3, 
                  fit_akmeans = res_ak, fit_akmeans_nodim = res_ak_nodim))
      
    }
    
    ## Do fitting parallelized over data sets
    if (.Platform$OS.type == "windows") {
      cl <- makeCluster(mc.cores, type = "PSOCK")
    } else {
      cl <- makeCluster(mc.cores, type = "FORK")
    }
    registerDoParallel(cl)
    res <- foreach (i = seq_along(dat)) %dopar% analyze(dat[[i]], delta = delta, nclust = nclust, 
                                                        ndim = ndim, alpha = alpha, fixed = fixed, 
                                                        parallel = FALSE, verbose = verbose, nstart = nstart, 
                                                        nstart.kmeans = nstart.kmeans)
    stopCluster(cl)
    # res <- mclapply(X = dat, FUN = analyze, mc.cores = mc.cores)
    
    ## Combine into tables
    tab_lsb_ov <- t(sapply(res, "[[", "lsbclust_overall"))
    tab_lsb_rows <- t(sapply(res, "[[", "lsbclust_rows"))
    tab_lsb_col <- t(sapply(res, "[[", "lsbclust_columns"))
    tab_lsb_int <- t(sapply(res, "[[", "lsbclust_interactions"))
    tab_t3 <- t(sapply(res, "[[", "t3clus"))
    tab_ak <- t(sapply(res, "[[", "akmeans"))
    tab_ak_nodim <- t(sapply(res, "[[", "akmeans_nodim"))
    if (include_fits) {
      fit_lsbclust <- lapply(res, "[[", "fit_lsbclust")
      fit_t3clus <- lapply(res, "[[", "fit_t3clus")
      fit_akmeans <- lapply(res, "[[", "fit_akmeans")
      fit_akmeans_nodim <- lapply(res, "[[", "fit_akmeans_nodim")
    }
    
    ## Stop time
    t1 <- proc.time()[3]
    
    if (include_fits) {
      return(list(lsbclust_overall = tab_lsb_ov, 
                  lsbclust_rows = tab_lsb_rows,
                  lsbclust_columns = tab_lsb_col,
                  lsbclust_interactions = tab_lsb_int, 
                  t3clus = tab_t3, 
                  akmeans = tab_ak,
                  akmeans_nodim = tab_ak_nodim,
                  fit_lsbclust = fit_lsbclust,
                  fit_t3clus = fit_t3clus, 
                  fit_akmeans = fit_akmeans,
                  fit_akmeans_nodim = fit_akmeans_nodim, 
                  data = ifelse(include_data, dat, NA),
                  time = t1 - t0))
    } else {
      return(list(lsbclust_overall = tab_lsb_ov, 
                  lsbclust_rows = tab_lsb_rows,
                  lsbclust_columns = tab_lsb_col,
                  lsbclust_interactions = tab_lsb_int, 
                  t3clus = tab_t3, 
                  akmeans = tab_ak,
                  akmeans_nodim = tab_ak_nodim,
                  data =  ifelse(include_data, dat, NA),
                  time = t1 - t0))
    }
    
    
  } else {
    ## Parallelize over random starts
    ## Set progress bar
    # pb <- txtProgressBar(max = ndata, style = 3)
    res_lsb <- eval_lsb <- res_t3 <- eval_t3 <- 
      res_ak <- eval_ak <- res_ak_nodim <- eval_ak_nodim <- vector(length = ndata, mode = "list")
    
    ## Fit models and calculate fit
    for (i in seq_len(ndata)) {
      res_lsb[[i]] <- lsbclust(data = dat[[i]]$data, delta = delta, nclust = nclust, 
                               ndim = ndim, alpha = alpha, fixed = fixed, 
                               parallel = parallel, verbose = verbose, 
                               nstart = nstart, nstart.kmeans = nstart.kmeans)
      eval_lsb[[i]] <- cfsim(fitted = res_lsb[[i]], actual = dat[[i]])
      
      ## Double-centre array for T3Clusf
      dat[[i]]$data <- carray(dat[[i]]$data)
      res_t3[[i]] <- T3Clusf(X = dat[[i]]$data, Q = ndim, G = nclust[4L], 
                             parallel = parallel, nstart = nstart_T3, verbose = verbose)
      eval_t3[[i]] <- cfsim(fitted = res_t3[[i]], actual = dat[[i]])
      
      ## Do akmeans with ndim
      res_ak[[i]] <- akmeans(data = dat[[i]]$data, centers = nclust[4L], ndim = ndim, 
                             nstart = nstart_ak)
      eval_ak[[i]] <- cfsim(fitted = res_ak[[i]], actual = dat[[i]])
      
      ## Do akmeans without ndim
      res_ak_nodim[[i]] <- akmeans(data = dat[[i]]$data, centers = nclust[4L], ndim = NULL, 
                                   nstart = nstart_ak)
      eval_ak_nodim[[i]] <- cfsim(fitted = res_ak_nodim[[i]], actual = dat[[i]])
      
      # setTxtProgressBar(pb, value = i)
    }
    
    ## Convert results to table(s)
    ## LSBCLUST: interactions
    tab_lsb_int <- t(sapply(eval_lsb, function(x) x$interactions$stats))
    
    ## LSBCLUST: overall
    odelta <- delta[1] * delta[3] + delta[2] * delta[4] - delta[1] * delta[2]
    if (odelta) {
      tab_lsb_ov <- t(sapply(eval_lsb, function(x) x$overall$stats))    
    } else {
      tab_lsb_ov <- NULL    
    } 
    
    ## LSBCLUST: rows
    if (delta[2L]) {
      tab_lsb_rows <- t(sapply(eval_lsb, function(x) x$rows$stats))  
    } else {
      tab_lsb_rows <- NULL
    }
    
    ## LSBCLUST: columns
    if (delta[1L]) {
      tab_lsb_col <- t(sapply(eval_lsb, function(x) x$columns$stats))
    } else {
      tab_lsb_col <- NULL
    }
    
    ## T3Clusf
    tab_t3 <- t(sapply(eval_t3, function(x) x$interactions$stats))
    
    ## akmeans with ndim
    tab_ak <- t(sapply(eval_ak, function(x) x$interactions$stats))
    
    ## akmenas without ndim
    tab_ak_nodim <- t(sapply(eval_ak_nodim, function(x) x$interactions$stats))
    
    ## Stop time
    t1 <- proc.time()[3]
    
    if (include_fits) {
      return(list(lsbclust_overall = tab_lsb_ov, 
                  lsbclust_rows = tab_lsb_rows,
                  lsbclust_columns = tab_lsb_col,
                  lsbclust_interactions = tab_lsb_int, 
                  t3clus = tab_t3, 
                  akmeans = tab_ak,
                  akmeans_nodim = tab_ak_nodim,
                  fit_lsbclust = res_lsb,
                  fit_t3clus = res_t3, 
                  fit_akmeans = res_ak,
                  fit_akmeans_nodim = res_ak_nodim,
                  data = ifelse(include_data, dat, NA),
                  time = t1 - t0))
      
    } else {
      return(list(lsbclust_overall = tab_lsb_ov, 
                  lsbclust_rows = tab_lsb_rows,
                  lsbclust_columns = tab_lsb_col,
                  lsbclust_interactions = tab_lsb_int, 
                  t3clus = tab_t3, 
                  akmeans = tab_ak,
                  akmeans_nodim = tab_ak_nodim,
                  data = ifelse(include_data, dat, NA),
                  time = t1 - t0))
    }

  }
}