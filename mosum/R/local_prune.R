#' Multiscale MOSUM algorithm with localised pruning
#' 
#' Multiscale MOSUM procedure with (possibly) assymetric bandwidths and localised pruning based on Schwarz criterion.
#' @param x input data (a \code{numeric} vector or an object of classes \code{ts} and \code{timeSeries})
#' @param G a vector of bandwidths, given as either integers less than \code{length(x)/2}, 
#' or numbers between \code{0} and \code{0.5} describing the moving sum bandwidths relative to \code{length(x)}.
#' Asymmetric bandwidths obtained as the Cartesian product of the set \code{G} with itself are used for change-point analysis
#' @param max.unbalance a numeric value for the maximal ratio between maximal and minimal bandwidths to be used for candidate generation,
#' \code{1 <= max.unbalance <= Inf}
#' @param threshold string indicating which threshold should be used to determine significance.
#' By default, it is chosen from the asymptotic distribution at the significance level \code{alpha}.
#' Alternatively, it is possible to parse a user-defined function with \code{threshold.function}
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}. Use iff \code{threshold = "critical.value"}
#' @param threshold.function function object of form \code{function(G_l, G_r, length(x), alpha)}, to compute a
#' threshold of significance for different bandwidths \code{(G_l, G_r)}; use iff \code{threshold = "custom"}
#' @param criterion how to determine whether an exceeding point is a change-point; to be parsed to \link[mosum]{mosum}
#' @param eta,epsilon see \link[mosum]{mosum}
#' @param rule string for the choice of sorting criterion for change-point candidates in merging step. 
#' Possible values are: 
#' \itemize{
#' \item{\code{"pval"}}{smallest p-value}
#' \item{\code{"jump"}}{largest (rescaled) jump size}
#' }
#' @param penalty string specifying the type of penalty term to be used in Schwarz criterion; possible values are:
#' \itemize{
#' \item{\code{"log"}}{use \code{penalty = log(length(x))^pen.exp}}
#' \item{\code{"polynomial"}}{use \code{penalty = length(x)^pen.exp}}
#' }
#' @param pen.exp exponent for the penalty term (see \code{penalty});
#' @param do.confint flag indicating whether confidence intervals for change-points should be computed
#' @param level use iff \code{do.confint = TRUE}; a numeric value (\code{0 <= level <= 1}) with which
#' \code{100(1-level)\%} confidence interval is generated
#' @param N_reps use iff \code{do.confint = TRUE}; number of bootstrap replicates to be generated
#' @param ... further arguments to be parsed to \link[mosum]{mosum} calls
#' @return S3 object of class \code{multiscale.cpts}, which contains the following fields:
#'    \item{x}{input data}
#'    \item{cpts}{estimated change-points}
#'    \item{cpts.info}{data frame containing information about estimated change-points}
#'    \item{sc}{Schwarz criterion values of the estimated change-point set}
#'    \item{pooled.cpts}{set of change-point candidates that have been considered by the algorithm}
#'    \item{G}{input parameter}
#'    \item{threshold, alpha, threshold.function}{input parameters}
#'    \item{criterion, eta, epsilon}{input parameters}
#'    \item{rule, penalty, pen.exp}{input parameters}
#'    \item{do.confint}{input parameter}
#'    \item{ci}{object of class \code{cpts.ci} containing confidence intervals for change-points iff \code{do.confint = TRUE}}
#' @details See Algorithm 2 in the first referenced paper for a comprehensive
#' description of the procedure and further details.
#' @references A. Meier, C. Kirch and H. Cho (2019)
#' mosum: A Package for Moving Sums in Change-point Analysis. \emph{To appear in the Journal of Statistical Software}.
#' @references H. Cho and C. Kirch (2019)
#' Localised pruning for data segmentation based on multiscale change point procedures. \emph{arXiv preprint arXiv:1910.12486}.
#' @examples 
#' x <- testData(model = "mix", seed = 123)$x
#' mlp <- multiscale.localPrune(x, G = c(8, 15, 30, 70), do.confint = TRUE)
#' print(mlp)
#' summary(mlp)
#' par(mfcol=c(2, 1), mar = c(2, 4, 2, 2))
#' plot(mlp, display = "data", shaded = "none")
#' plot(mlp, display = "significance", shaded = "CI", CI = "unif")
#' @importFrom Rcpp evalCpp
#' @useDynLib mosum, .registration = TRUE
#' @export
multiscale.localPrune <- function(x, G=bandwidths.default(length(x)), max.unbalance = 4,
                            threshold=c('critical.value', 'custom')[1], alpha=.1, threshold.function = NULL,
                            criterion=c('eta', 'epsilon')[1], eta=0.4, epsilon=0.2,
                            rule=c('pval', 'jump')[1], penalty=c('log', 'polynomial')[1], pen.exp=1.01,
                            do.confint=FALSE, level=0.05, N_reps=1000, ...) {
  
  n <- length(x)
  
  if (class(G) == "integer" || class(G) == "numeric") {
    grid <- multiscale.grid(G, max.unbalance = max.unbalance)
  } else if(class(G) == 'multiscale.grid'){
    grid <- G
  } else stop('Expecting a vector of numbers')
  abs.bandwidth <- all(grid$grid>=1)
  
  stopifnot(max.unbalance >= 1)
  stopifnot(is.element(rule, c('pval', 'jump', 'lr', 'rl')))
  stopifnot(is.element(criterion, c('eta', 'epsilon')))
  stopifnot((criterion=='eta' & eta <= 1 & eta > 0) || (criterion=='epsilon' & epsilon <= 1 & epsilon > 0))
  stopifnot(!do.confint || N_reps>0)
  
  if (penalty == 'log') {
    log.penalty <- TRUE
  } else {
    if (penalty != 'polynomial') {
      stop('penalty has to set to log or polynomial')
    }
    log.penalty <- FALSE
  }
  
  if (threshold != 'critical.value' && threshold != 'custom') {
    stop('threshold must be either \'critical.value\' or \'custom\'')
  }
  
  all.cpts <- matrix(NA, ncol=6, nrow=0)
  for (i in seq_len(nrow(grid$grid))) {
    G1 <- grid$grid[[i,1]]
    G2 <- grid$grid[[i,2]]
    if (threshold == 'critical.value') {
      m <- mosum(x, G=G1, G.right=G2, ...,  
                        threshold='critical.value', alpha=alpha, 
                        criterion=criterion, eta=eta, epsilon=epsilon)
    } else{
      threshold_val <- threshold.function(G1, G2, n, alpha)
      m <- mosum(x, G=G1, G.right=G2, ..., 
                         threshold='custom', threshold.custom=threshold_val, alpha=alpha,
                         criterion=criterion, eta=eta, epsilon=epsilon)
    }

    if(length(m$cpts)>0){
      if (!abs.bandwidth) {
        G1 <- floor(G1*n)
        G2 <- floor(G2*n)
      }
      all.cpts <- rbind(all.cpts, 
                        cbind(m$cpts, G1, G2, G1+G2,
                              mosum.pValue(m$stat[m$cpts], n, G1, G2), m$stat[m$cpts]*sqrt(G1+G2)/sqrt(G1*G2)))
    }
  }

  all.cpts <- all.cpts[sort(all.cpts[, 1], decreasing=FALSE, index.return=TRUE)$ix,,drop=FALSE]
  all.cpts <- dup.merge(all.cpts) # if there are duplicates, only select one according to 'rule'
  ac <- nrow(all.cpts)
  if(ac > 0){
    lp <- local.prune(x, all.cpts, rule, log.penalty, pen.exp)
    est.cpts <- lp$est.cpts; est.cpts.ind <- lp$est.cpts.ind; min.cost <- lp$min.cost
  } else{
    est.cpts.ind <- est.cpts <- integer(0)
    min.cost <- sum(x^2) - n*mean(x)^2
  }
  
  est.cpts.info <- data.frame(cpts = all.cpts[est.cpts.ind, 1], 
                          G.left =  all.cpts[est.cpts.ind, 2], 
                          G.right =  all.cpts[est.cpts.ind, 3],
                          p.value = all.cpts[est.cpts.ind, 5],
                          jump = all.cpts[est.cpts.ind, 6])
  if (log.penalty) {
    penalty_term <- length(est.cpts)*log(n)^pen.exp
  } else {
    penalty_term <- length(est.cpts)*n^pen.exp
  }
  final.sc <- n/2*log(min.cost/n) + penalty_term
  if(!abs.bandwidth) G <- floor(n * sort(unique(c(grid$grid))))
  
  ret <- structure(list(x = x,
                        cpts = est.cpts, 
                        cpts.info = est.cpts.info,
                        sc = final.sc, 
                        pooled.cpts = all.cpts[, 1], 
                        G = G,
                        alpha = alpha,
                        threshold = threshold,
                        threshold.function = threshold.function,
                        criterion = criterion,
                        eta = eta,
                        epsilon = epsilon,
                        rule = rule,
                        penalty = penalty,
                        pen.exp = pen.exp,
                        do.confint = FALSE,
                        ci = NA), 
                   class = 'multiscale.cpts')
  if (do.confint) {
    ret$ci <- confint.multiscale.cpts(ret, level = level, N_reps = N_reps)
    ret$do.confint <- TRUE
  }
  return(ret)
}

#' Localised pruning algorithm
#' @keywords internal
local.prune <- function(x, all.cpts, rule, log.penalty, pen.exp){
  
  THRESH_MANUAL_MERGING <- 24
  
  n <- length(x)
  ac <- dim(all.cpts)[1]
  
  all.cpts <- cbind(all.cpts, 1:ac)
  cand_used <- rep(FALSE, ac)
  all.unique.cpts <- c(0, all.cpts[, 1], n)
  auc <- length(all.unique.cpts) - 2
  sums <- matrix(0, nrow = auc+1, ncol = 4) # calculated for efficient computation of rss
  for(j in 1:(auc + 1)){
    sums[j, 1:2] <- c(all.unique.cpts[j] + 1, all.unique.cpts[j + 1])
    sums[j, 3] <- sum(x[sums[j, 1]:sums[j, 2]])
    sums[j, 4] <- sum(x[sums[j, 1]:sums[j, 2]]^2)
  }
  min.cost <- sum(sums[, 4] - sums[, 3]^2/(sums[, 2] - sums[, 1]+1)) # min rss with all the candidates
  
  if(rule == 'pval'){
    u <- all.cpts[order(all.cpts[, 5], all.cpts[, 4], all.cpts[, 2], all.cpts[, 3]),, drop = FALSE]
    rule.seq <- u[, 7]; rm(u)
  }
  if(rule == 'jump'){
    u <- all.cpts[order(-all.cpts[, 6], all.cpts[, 4], all.cpts[, 2], all.cpts[, 3]),, drop = FALSE]
    rule.seq <- u[, 7]; rm(u)
  }
  if(rule == 'lr') rule.seq <- pool
  if(rule == 'rl') rule.seq <- rev(pool)
  
  current <- pool <- seq_len(ac); est.cpts.ind <- est.cpts <- integer(0)
  # current = C, pool = P, est.cpts.ind = B (index)
  while(length(pool)>0){
    # step 1
    j <- rule.seq[1]; adj <- 0
    
    # step 2
    le <- local.env(j, est.cpts.ind, all.cpts, current, ac)
    li <- le$li; li_final <- le$li_final
    ri <- le$ri; ri_final <- le$ri_final
    
    #step 3
    left <- li + 1
    right <- ri - 1
    cand.ind <- (left:right)[is.element(left:right, pool)]
    cand <- all.cpts[cand.ind, 1] # = D
    
    # left <- max(est.cpts.ind[est.cpts.ind < j]+1, j - all.cpts[j, 7])
    # right <- min(est.cpts.ind[est.cpts.ind > j]-1, j + all.cpts[j, 8])
    # if(sum(current < left) > 0) li <- max(current[current < left]) else li <- 0
    # if(sum(current > right) > 0) ri <- min(current[current > right]) else ri <- ac+1  
    
    ind_middl_tmp <- sums[(li + 1):(ri - 1), 2]
    ind_middl_tmp <- ind_middl_tmp[which(!cand_used[(li + 1):(ri - 1)])]
    ind_tmp <- c(sums[li + 1, 1] - 1, ind_middl_tmp, sums[ri, 2])
    sub.sums <- extract_sub(ind_tmp, x)
    
    doExhaustiveSearch <- TRUE
    # Too many candidates to do exhaustive search?
    if (length(cand) > THRESH_MANUAL_MERGING) {
      # Count neighbourhood size of neighbours
      
      cand.rule.seq <- rule.seq[is.element(rule.seq, cand.ind)]
      cand_size <- rep(NA, length(cand))
      cand_size[1] <- length(cand)
      for (i_tmp in 2:length(cand)) {
        jj <- cand.rule.seq[i_tmp]
        le_jj <- local.env(jj, est.cpts.ind, all.cpts, current, ac)
        left_jj <- le_jj$li + 1
        right_jj <- le_jj$ri - 1
        cand.ind_jj <- (left_jj:right_jj)[is.element(left_jj:right_jj, pool)]
        #cand_jj <- all.cpts[cand.ind_jj, 1] # = D
        cand_size[i_tmp] <- length(cand.ind_jj)
      }
      
      if (any(cand_size <= THRESH_MANUAL_MERGING)) {
        # Proceed with next candidate, for which exhaustive search IS possible
        rule_tmp <- cand.rule.seq[min(which(cand_size <= THRESH_MANUAL_MERGING))]
        ind_star <- which(rule.seq == rule_tmp)
        rule.seq[ind_star] <- rule.seq[1]; rule.seq[1] <- rule_tmp
        doExhaustiveSearch <- FALSE
      } else {
        # Count neighbourhood size of remaining candidates
        cand_size <- rep(NA, length(rule.seq))
        cand_size[1] <- length(cand)
        for (i_tmp in seq(from = 2, length.out = length(rule.seq) - 1)) {
          jj <- rule.seq[i_tmp]
          le_jj <- local.env(jj, est.cpts.ind, all.cpts, current, ac)
          left_jj <- le_jj$li + 1
          right_jj <- le_jj$ri - 1
          cand.ind_jj <- (left_jj:right_jj)[is.element(left_jj:right_jj, pool)]
          #cand_jj <- all.cpts[cand.ind_jj, 1] # = D
          cand_size[i_tmp] <- length(cand.ind_jj)
        }
        
        if (any(cand_size <= THRESH_MANUAL_MERGING)) {
          # Proceed with next candidate, for which exhaustive search IS possible
          ind_star <- min(which(cand_size <= THRESH_MANUAL_MERGING))
          rule_tmp <- rule.seq[ind_star]; rule.seq[ind_star] <- rule.seq[1]; rule.seq[1] <- rule_tmp
          doExhaustiveSearch <- FALSE
        } else {
          # No more exhaustive search possible at all
          # --> Do manual merging, until exhaustive search becomes possible
          while(length(cand) > THRESH_MANUAL_MERGING) {
            warn_msg <- paste0('Warning: ', length(cand), ' conflicting candidates, thinning manually')
            print(warn_msg) #; warning(warn_msg)
            k <- cand[which.min(diff(cand))]
            l <- which(sub.sums[, 2] == k)
            a <- sub.sums[l, ]; b <- sub.sums[l + 1, ]
            # as change-points are merged, the minimum rss in the local environment needs to be updated
            adj <- adj + (a[2] - a[1] + 1)*(b[2] - b[1] + 1)/(b[2] - a[1] + 1)*(a[3]/(a[2] - a[1] + 1) - b[3]/(b[2] - b[1] + 1))^2
            sub.sums[l + 1, 1] <- a[1]; sub.sums[l + 1, 3:4] <- sub.sums[l + 1, 3:4]+a[3:4]
            sub.sums <- sub.sums[-l,, drop = FALSE]
            cand <- setdiff(cand, k)
            k.ind <- which(all.cpts[, 1] == k)
            cand.ind <- setdiff(cand.ind, k.ind); pool <- setdiff(pool, k.ind); rule.seq <- setdiff(rule.seq, k.ind)
            cand_used[k.ind] <- TRUE
          }
        }
      }
    }
    
    if (doExhaustiveSearch) {
      # step 4
      # performs exhaustive search (Algorithm 2)
      out <- exhaust_sc(cand = cand, sub_sums = sub.sums, 
                         strength = pen.exp, log_penalty = log.penalty, 
                         n = n, auc = length(current), min_cost = min.cost)
      est.cpts <- c(est.cpts, out$est_cpts)
      current.est.cpts.ind <- all.cpts[all.cpts[, 1] %in% out$est_cpts, 7]
      est.cpts.ind <- c(est.cpts.ind, current.est.cpts.ind)
      
      # steps 5, 6
      # removal of candidates
      rm.set <- c(j, current.est.cpts.ind)
      if(length(current.est.cpts.ind) > 0){
        rm.set <- c(rm.set, cand.ind[cand.ind %in% min(current.est.cpts.ind):max(current.est.cpts.ind)])
        if(li_final) rm.set <- c(rm.set, cand.ind[cand.ind <= max(current.est.cpts.ind)])
        if(ri_final) rm.set <- c(rm.set, cand.ind[cand.ind >= min(current.est.cpts.ind)])
      }        
      
      # tmp <- (all.cpts[cand.ind, 1] - sub.sums[1, 1] >= all.cpts[cand.ind, 2]) & 
      #   (sub.sums[nrow(sub.sums), 2] - all.cpts[cand.ind, 1] >= all.cpts[cand.ind, 3])
      # if(li > 0) tmp <- tmp & all.cpts[cand.ind, 1] - all.cpts[li, 1] >= all.cpts[li, 3]
      # if(ri < ac + 1) tmp <- tmp & all.cpts[ri, 1] - all.cpts[cand.ind, 1] >= all.cpts[ri, 2]
      # rm.set <- c(rm.set, cand.ind[tmp])
      # 
      # if(length(current.est.cpts.ind) > 0){
      #   rm.set <- c(rm.set, cand.ind[cand.ind %in% min(current.est.cpts.ind):max(current.est.cpts.ind)])
      #   if(li_final) rm.set <- c(rm.set, cand.ind[cand.ind <= max(current.est.cpts.ind)])
      #   if(ri_final) rm.set <- c(rm.set, cand.ind[cand.ind >= min(current.est.cpts.ind)])
      # }
      # rm.set <- min(rm.set):max(rm.set)
      
      pool <- setdiff(pool, rm.set)
      cand_used[rm.set] <- TRUE
      rule.seq <- setdiff(rule.seq, rm.set)
      current <- c(pool, est.cpts.ind)
      current_cands <- is.element(cand, all.cpts[current, 1])
      ind_star <- get_comb_ind(current_cands)
      min.cost <- min.cost + adj - out$sc[nrow(out$sc), 1] + out$sc[ind_star + 1, 1]
    }
  }
  est.cpts <- sort(as.vector(est.cpts)); est.cpts.ind <- sort(as.vector(est.cpts.ind))
  
  return(list(est.cpts = est.cpts, est.cpts.ind = est.cpts.ind, min.cost = min.cost))
}
 
#' Identify the local environment for exhaustive search
#' @keywords internal
local.env <- function(j, est.cpts.ind, all.cpts, current, ac){
  li_final <- ri_final <- TRUE
  if(sum(est.cpts.ind < j)) li <- max(est.cpts.ind[est.cpts.ind < j]) else li <- 0
  if(j > 1){
    ind <- ((li + 1):(j - 1))[(li + 1):(j - 1) %in% current]
    if(length(ind) > 0 && sum(tmp <- all.cpts[j, 1] - all.cpts[ind, 1] >= pmax(all.cpts[j, 2], all.cpts[ind, 3]))) li_tmp <- max(ind[tmp]) else li_tmp <- li
    if(li_tmp > li){ li <- li_tmp; li_final <- FALSE }
  }
  
  if(sum(est.cpts.ind > j)) ri <- min(est.cpts.ind[est.cpts.ind > j]) else ri <- ac + 1
  if(j < ac){
    ind <- ((j + 1):(ri - 1))[(j + 1):(ri - 1) %in% current]
    if(length(ind) > 0 && sum(tmp <- all.cpts[ind, 1] - all.cpts[j, 1] >= pmax(all.cpts[j, 3], all.cpts[ind, 2]))) ri_tmp <- min(ind[tmp]) else ri_tmp <- ri
    if(ri_tmp < ri){ ri <- ri_tmp; ri_final <- FALSE }
  }
  
  list(li = li, li_final = li_final, ri = ri, ri_final = ri_final)
}

#' Remove duplicated from all.cpts data frame:
#' In case one change being added multiple times, choose the one
#' with smallest p-value 
#' @keywords internal
dup.merge <- function(all.cpts) {
  all.unique.cpts <- unique(all.cpts[, 1, drop=FALSE])
  out <- matrix(NA, nrow=0, ncol=ncol(all.cpts))
  for(k in all.unique.cpts){
    ind <- which(all.cpts[, 1]==k)
    ind.min <- ind[all.cpts[ind, 5]==min(all.cpts[ind, 5])]
    if(length(ind.min) > 1) ind.min <- ind.min[which.min(all.cpts[ind.min, 4])]
    out <- rbind(out, all.cpts[ind.min,])
  }
  out
}

# dup.merge <- function(all.cpts, rule='jump') {
#   all.unique.cpts <- unique(all.cpts[, 1, drop=FALSE])
#   out <- matrix(NA, nrow=0, ncol=ncol(all.cpts))
#   for(k in all.unique.cpts){
#     ind <- which(all.cpts[, 1]==k)
#     ind.min <- ind[all.cpts[ind, 4]==min(all.cpts[ind, 4])]
#     if(length(ind.min) > 1 & rule=='pval') ind.min <- ind.min[which.min(all.cpts[ind.min, 5])]
#     if(length(ind.min) > 1 & rule=='jump') ind.min <- ind.min[which.max(all.cpts[ind.min, 6])]
#     out <- rbind(out, all.cpts[ind.min,])
#   }
#   out
# }

#' Plotting the output from multiscale MOSUM procedure
#' 
#' Plotting method for S3 objects of class "multiscale.cpts".
#' @method plot multiscale.cpts
#' @param x a \code{multiscale.cpts} object
#' @param display which to be plotted against the estimated change-point locations; possible values are
#' \itemize{
#'    \item{\code{"data"}}{input time series is plotted along with the estimated piecewise constant signal}
#'    \item{\code{"significance"}}{one minus the p-values associated with the detection of change-point estimators
#'    are represented as the height of vertical lines indicating their locations}
#' }
#' @param shaded string indicating which to display as shaded areas surrounding the estimated change-point locations.
#' Poissble values are 
#' \itemize{
#'    \item{\code{"bandwidth"}}{respective detection intervals are plotted} 
#'    \item{\code{"CI"}}{bootstrap confidence intervals are plotted}
#'    \item{\code{"none"}}{none is plotted}
#' }
#' @param level,N_reps argument to be parsed to \link[mosum]{confint.multiscale.cpts}; use iff \code{shaded = "CI"}.
#' @param CI string indicating whether pointwise (\code{CI = "pw"}) or uniform (\code{CI = "unif"}) confidence intervals
#' are to be plotted; use iff \code{shaded = "CI"} 
#' @param xlab graphical parameter
#' @param ... not in use
#' @details
#' The locations of change-point estimators are plotted 
#' against the input time series and the estimated piecewise constant signal (\code{display = "data"}), or 
#' the significance of each estimator is represented by the corresponding
#' \code{1-p.value} derived from the asymptotic distribution of MOSUM test statistic (\code{display = "significance"}).
#' It also produces the rectangles representing the 
#' detection intervals (if \code{shaded = "bandwidth"}) or 
#' bootstrap confidence intervals of the corresponding change-points (if \code{shaded = "CI"})
#' around their locations.
#' @examples 
#' x <- testData(model = "blocks", seed = 1234)$x
#' mlp <- multiscale.localPrune(x)
#' par(mfrow = c(2, 1))
#' plot(mlp, display = "data", shaded = "bandwidth")
#' plot(mlp, display = "significance", shaded = "CI")
#' @importFrom grDevices rainbow
#' @importFrom graphics abline axis plot points rect segments par lines
#' @export
plot.multiscale.cpts <- function(x, display = c('data', 'significance')[1],
                                 shaded=c('CI', 'bandwidth', 'none')[1], 
                                 level=0.05, N_reps=1000, 
                                 CI = c('pw', 'unif')[1], xlab = 'Time', ...) {
  if (shaded=='bandwidth') {
    main <- 'Change-point estimates and detection intervals'
  } else if (shaded=='CI') {
    if(CI=='pw') main <- paste('Change-point estimates and pointwise ', 100*(1 - level), '% confidence intervals')
    if(CI=='unif') main <- paste('Change-point estimates and uniform ', 100*(1 - level), '% confidence intervals')
    if(length(x$cpts) > 0) if(x$do.confint) b <- x$ci else b <- confint.multiscale.cpts(x, level=level, N_reps=N_reps)
  } else if (shaded == 'none') {
    main <- 'Change-point estimates'
  } else {
    stop('shaded argument has to be either \'CI\', \'bandwidth\' or \'none\'.')
  }
  n <- length(x$x)
  if (class(x$x)=='ts') {
    x_plot <- as.numeric(time(x$x))
  } else if(class(x$x)=='timeSeries') {
    x_plot <- time(x$x)
  } else {
    x_plot <- seq_len(length(x$x))
  }
  if(length(x$cpts) > 0){
    cpts <- x$cpts.info
    cols <- rainbow(nrow(cpts),alpha=0.2)
    cols2 <- rainbow(nrow(cpts),alpha=1)
    # cols3 <- rainbow(nrow(cpts),alpha=0.2)
    newOrder <- 1:nrow(cpts) #sample(seq_len(nrow(cpts)), nrow(cpts))
    cols <- cols[newOrder]
    cols2 <- cols2[newOrder]
  
    xx <- cpts$cpts # location
    if (shaded=='bandwidth') {
      xx_l <- pmax(1, xx-cpts$G.left+1)
      xx_r <- pmin(n, xx+cpts$G.right)
    }
    if (shaded == 'CI') {
      if(CI == 'pw'){
        xx_l <- b$CI[,2]
        xx_r <- b$CI[,3]
      } else if(CI == 'unif'){
        xx_l <- b$CI[,4]
        xx_r <- b$CI[,5]
      }
    }
  }
  
  if(display == 'data'){
    brks <- c(0, x$cpts, length(x$x))
    fhat <- x$x * 0
    for(kk in 1:(length(brks) - 1)){
      int <- (brks[kk] + 1):brks[kk + 1]
      fhat[int] <- mean(x$x[int])
    }
    plot(x_plot, x$x, type='l', xlab = xlab, ylab = expression(x[t]), main = main)
    lines(x_plot, fhat, col = 'darkgray', type = 'l', lwd = 2)
    if(length(x$cpts) > 0){
      segments(x0 = x_plot[xx], y0 = par("usr")[3] - 1, y1 = par("usr")[4] + 1, col = cols2)
      if(shaded != 'none') rect(xleft = x_plot[xx_l], xright = x_plot[xx_r], ybottom = par("usr")[3] - 1, ytop = par("usr")[4] + 1, col = cols, lty = 0)
    }
  }
  if(display == 'significance'){
    y.min <- max(1 - 1.1*x$alpha, 0)  
    plot(0, type = 'n', xlim = range(x_plot), ylim = c(y.min, 1), xlab = xlab, ylab = expression(1 - p.value), main = main)
    if(length(x$cpts) > 0){
      yy <- 1 - cpts$p.value # pvalue
      points(x_plot[xx], yy)
      segments(x0 = x_plot[xx], y0 = 0, y1 = yy, col = cols2)
      if(shaded != 'none') rect(xleft = x_plot[xx_l], xright = x_plot[xx_r], ybottom = 0, ytop = yy, col = cols, lty = 0)
    }
    # if (shaded=='CI' && include_unif_CI) {
    #   # also add uniform ones
    #   xx_ll <- b$CI[,4]
    #   xx_rr <- b$CI[,5]
    #   rect(xleft=xx_ll, xright=xx_rr, ybottom=0, ytop=yy, col=cols3, lty=0)
    # }
  }
}

#' Summary of change-points estimated by multiscale MOSUM procedure
#' 
#' Summary method for objects of class \code{multiscale.cpts}
#' @method summary multiscale.cpts
#' @param object a \code{multiscale.cpts} object
#' @param ... not in use
#' @details Provide information about each estimated change-point, 
#' including the bandwidths used for its detection, associated p-value and (scaled) jump size;
#' if \code{object$do.confint=TRUE}, end points of the pointwise and uniform confidence intervals
#' are also provided.
#' @examples 
#' x <- testData(model = "mix", seed = 12345)$x
#' mlp <- multiscale.localPrune(x, do.confint = TRUE)
#' summary(mlp)
#' @export
summary.multiscale.cpts <- function(object, ...) { 
  n <- length(object$x)
  if(length(object$cpts) > 0){
    ans <- object$cpts.info
    ans$p.value <- signif(ans$p.value, 3)
    ans$jump <- round(ans$jump, 3)
  } 
  if(object$do.confint) ans <- cbind(ans, object$ci$CI[, -1, drop=FALSE])
  #cat(paste('created using mosum version ', utils::packageVersion('mosum'), sep=''))
  cat(paste('change-points estimated at alpha = ', object$alpha, ' according to ', object$criterion, '-criterion', sep=''))
  if(object$criterion=='eta') cat(paste('\n with eta = ', object$eta, sep=''))
  if(object$criterion=='epsilon') cat(paste('\n with epsilon = ', object$epsilon, ':', sep=''))
  cat('\n')
  cat('\n')
  if(length(object$cpts) > 0) print(ans, print.gap = 3) else cat('no change-point is found') 
  cat('\n')
}

#' Change-points estimated by multiscale MOSUM procedure
#' 
#' Print method for objects of class \code{multiscale.cpts}
#' @method print multiscale.cpts
#' @param x a \code{multiscale.cpts} object
#' @param ... not in use
#' @examples 
#' x <- testData(model = "mix", seed = 12345)$x
#' mlp <- multiscale.localPrune(x)
#' print(mlp)
#' @export
print.multiscale.cpts <- function(x, ...) {
  #cat(paste('created using mosum version ', utils::packageVersion('mosum'), sep=''))
  cat(paste('change-points estimated with bandwidths\n'))
  cat('  ')
  cat(x$G)
  cat(paste('\nat alpha = ', x$alpha, ' according to ', x$criterion, '-criterion', sep=''))
  if(x$criterion=='eta') cat(paste(' with eta = ', x$eta, ':', sep=''))
  if(x$criterion=='epsilon') cat(paste(' with epsilon = ', x$epsilon, ':', sep=''))
  cat('\n')
  cat('\n')
  cat('  ')
  if(length(x$cpts)==0) cat('no change-point is found') else cat(x$cpts)
  cat('\n')
}