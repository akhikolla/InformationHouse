#' Calculate quantities for supplemental data
#'
#' @param CGGP CGGP object
#' @param thetaMAP theta to calculate for
#' @param ys.thisloop Supplemental ys (after transformation)
#' @param y.thisloop y (after transformation)
#' @param pw pw for y
#' @param sigma2MAP sigma2MAP if only using y, will be adjusted for ys 
#' @param only_sigma2MAP Should only sigma2MAP be returned?
#'
#' @return List with objects
## @export
#' @noRd
CGGP_internal_calc_supp_pw_sigma2_Sti <- function(CGGP, thetaMAP,
                                                  ys.thisloop,
                                                  y.thisloop,
                                                  pw, sigma2MAP,
                                                  only_sigma2MAP=TRUE) {
  if (missing(pw) || missing(sigma2MAP)) {
    if (missing(y.thisloop)) {stop("If not giving in pw and sigma2MAP, must give in y.thisloop")}
    likstuff <- CGGP_internal_calc_cholS_lS_sigma2_pw(CGGP=CGGP, y.thisloop, theta=thetaMAP)
    pw <- likstuff$pw
    sigma2MAP <- likstuff$sigma2
  }
  
  Cs = matrix(1,dim(CGGP$Xs)[1],CGGP$ss)
  for (dimlcv in 1:CGGP$d) { # Loop over dimensions
    V = CGGP$CorrMat(CGGP$Xs[,dimlcv], CGGP$xb,
                     thetaMAP[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
    Cs = Cs*V[,CGGP$designindex[,dimlcv]]
  }
  
  Sigma_t = matrix(1,dim(CGGP$Xs)[1],dim(CGGP$Xs)[1])
  for (dimlcv in 1:CGGP$d) { # Loop over dimensions
    V = CGGP$CorrMat(CGGP$Xs[,dimlcv], CGGP$Xs[,dimlcv],
                     thetaMAP[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
    Sigma_t = Sigma_t*V
  }
  
  MSE_s = list(matrix(0,dim(CGGP$Xs)[1],dim(CGGP$Xs)[1]),
               (CGGP$d+1)*(CGGP$maxlevel+1)) 
  for (dimlcv in 1:CGGP$d) {
    for (levellcv in 1:max(CGGP$uo[1:CGGP$uoCOUNT,dimlcv])) {
      MSE_s[[(dimlcv)*CGGP$maxlevel+levellcv]] =
        (-CGGP_internal_postvarmatcalc(CGGP$Xs[,dimlcv],CGGP$Xs[,dimlcv],
                                       CGGP$xb[1:CGGP$sizest[levellcv]],
                                       thetaMAP[(dimlcv-1)*CGGP$numpara +
                                                  1:CGGP$numpara],
                                       CorrMat=CGGP$CorrMat))
    }
  }
  
  for (blocklcv in 1:CGGP$uoCOUNT) {
    ME_s = matrix(1,nrow=dim(CGGP$Xs)[1],ncol=dim(CGGP$Xs)[1])
    for (dimlcv in 1:CGGP$d) {
      levelnow = CGGP$uo[blocklcv,dimlcv]
      ME_s = ME_s*MSE_s[[(dimlcv)*CGGP$maxlevel+levelnow]]
    }
    Sigma_t = Sigma_t-CGGP$w[blocklcv]*(ME_s)
  }
  yhats = Cs%*%pw
  
  
  Sti_resid = solve(Sigma_t,ys.thisloop-yhats)
  sigma2MAP = (sigma2MAP*dim(CGGP$design)[1] +
                 colSums((ys.thisloop-yhats)*Sti_resid)) / (
                   dim(CGGP$design)[1]+dim(CGGP$Xs)[1])
  
  out <- list(sigma2MAP=sigma2MAP)
  
  if (!only_sigma2MAP) {
    out$Sti = solve(Sigma_t)
    pw_adj_y = t(Cs)%*%Sti_resid
    pw_adj <- CGGP_internal_calcpw(CGGP=CGGP, y=pw_adj_y, theta=thetaMAP)
    
    out$pw_uppadj = pw-pw_adj
    out$supppw = Sti_resid
  }
  
  
  out
}

#' Calculate quantities for supplemental data only (no grid data)
#'
#' @param CGGP CGGP object
#' @param thetaMAP theta to calculate for
#' @param ys.thisloop Supplemental ys (after transformation)
#' @param only_sigma2MAP Should only sigma2MAP be returned?
#'
#' @return List with objects
## @export
#' @noRd
CGGP_internal_calc_supp_only_supppw_sigma2_Sti <- function(CGGP, thetaMAP,
                                                       ys.thisloop,
                                                       only_sigma2MAP=TRUE) {
  
  Sigma_t = matrix(1,dim(CGGP$Xs)[1],dim(CGGP$Xs)[1])
  for (dimlcv in 1:CGGP$d) { # Loop over dimensions
    V = CGGP$CorrMat(CGGP$Xs[,dimlcv], CGGP$Xs[,dimlcv], thetaMAP[(dimlcv-1)*CGGP$numpara+1:CGGP$numpara])
    Sigma_t = Sigma_t*V
  }
  
  Sti_chol <- chol(Sigma_t + diag(CGGP$nugget, nrow(Sigma_t), ncol(Sigma_t)))
  
  # Use backsolve for stability
  # supppw <- Sti %*% ys.thisloop
  # supppw is same as Sti_resid in other files
  supppw <- backsolve(Sti_chol, backsolve(Sti_chol, ys.thisloop, transpose = T))
  if (is.matrix(supppw) && ncol(supppw)==1) {supppw <- as.vector(supppw)}
  
  # sigma2MAP <- (t(ys.thisloop) %*% supppw) / nrow(CGGP$Xs)
  # if (is.matrix(sigma2MAP)) {sigma2MAP <- diag(sigma2MAP)}
  sigma2MAP = colSums(as.matrix(ys.thisloop)*as.matrix(supppw))/nrow(CGGP$Xs)
  
  out <- list()
  
  out$sigma2MAP <- sigma2MAP
  
  if (!only_sigma2MAP) {
    Sti <- chol2inv(Sti_chol)
    out$Sti <- Sti
    out$supppw <- supppw
    out$Sti_chol <- Sti_chol
  }
  
  out
}
