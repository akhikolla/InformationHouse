#' Correlation skill analysis for ensemble forecasts
#' 
#' Calculate correlation between forecasts and observations for an ensemble forecast, including an adjustment for finite ensemble sizes
#'
#' @param ens a N*R matrix representing N time instances of real-valued R-member ensemble forecasts
#' @param obs a numeric vector of length N with real-valued observations
#' @param R.new positive number, can be `Inf`, ensemble size for which correlation skill should be estimated, default is NA for using the actual size R of the ensemble
#' @md
#' @return A vector with 4 entries: 
#' * `cmy`: Correlation skill of the ensemble mean forecast
#' * `cmy_adj`: Correlation skill of the ensemble mean forecast adjusted to ensemble size R.new
#' * `cxx`: Average correlation between ensemble members
#' * `cxy`: Average correlation between individual ensemble members and observation
#' @examples
#' data(eurotempforecast)
#' EnsCorr(ens, obs, R.new=Inf)
#' @seealso Corr, CorrDiff
#' @references Von Storch, Zwiers (2001): Statistical analysis in climate research. Cambridge University Press.
#'
#' Murphy (1990), Assessment of the practical utility of extended range ensemble forecasts, Q. J. R. Meteorol. Soc., 116, 89-125.
#' @export
EnsCorr = function(ens, obs, R.new=NA) {

  ## sanity checks
  stopifnot(is.matrix(ens), 
            nrow(ens) == length(obs), 
            length(R.new)==1)
  if (is.na(R.new)) {
    R.new = ncol(ens)
  }
  stopifnot(R.new > 0)


  N = length(obs)
  R = ncol(ens)
 

  ## calculate correlation of the ensemble mean
  m = rowMeans(ens, na.rm=TRUE)
  cmy = cor(obs, m)


  ## calculate average member-member correlation
  cor_mat = cor(ens, use='pairwise.complete.obs')
  cxx = mean(cor_mat[upper.tri(cor_mat)])


  ## calculate average member-obs correlation
  cxy = mean(cor(ens, obs, use='pairwise.complete.obs')) 


  ## calculate adjustement
  beta_sq = max(cxx, 0)
  adj_sq = (beta_sq * (1-1/R) + 1/R) / (beta_sq * (1-1/R.new) + 1/R.new)  
  cmy_adj = cmy * sqrt(adj_sq)
  cmy_adj = min(cmy_adj, 1)


  ## return
  ret = c(cmy, cmy_adj, cxx, cxy)
  names(ret) = c("cor_my", "cor_my_adj", "cor_xx", "cor_xy")

  return(ret)
}

