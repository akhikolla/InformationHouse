#' Set correlation function of CGGP object
#'
#' @param CGGP CGGP object
#' @param corr Correlation function
#'
#' @return CGGP object
#' @export
#'
#' @examples
#' obj <- CGGPcreate(3, 20, corr="matern52")
#' CGGP_internal_set_corr(obj, "gaussian")
CGGP_internal_set_corr <- function(CGGP, corr) {
  if (is.function(corr)) {
    CGGP$CorrMat <- corr
    CGGP$CorrName <- "UserDefined"
  } else if (tolower(corr) %in% c("cauchysqt")) {
    CGGP$CorrMat <- CGGP_internal_CorrMatCauchySQT
    CGGP$CorrName <- "CauchySQT"
  } else if (tolower(corr) %in% c("cauchysq")) {
    CGGP$CorrMat <- CGGP_internal_CorrMatCauchySQ
    CGGP$CorrName <- "CauchySQ"
  } else if (tolower(corr) %in% c("cauchy")) {
    CGGP$CorrMat <- CGGP_internal_CorrMatCauchy
    CGGP$CorrName <- "Cauchy"
  } else if (tolower(corr) %in% c("gaussian", "gauss", "sqexp")) {
    CGGP$CorrMat <- CGGP_internal_CorrMatGaussian
    CGGP$CorrName <- "Gaussian"
  } else if (tolower(corr) %in% c("powerexp", "pe", "powerexponential")) {
    CGGP$CorrMat <- CGGP_internal_CorrMatPowerExp
    CGGP$CorrName <- "PowerExponential"
  } else if (tolower(corr) %in% c("matern32", "m32", "m3")) {
    CGGP$CorrMat <- CGGP_internal_CorrMatMatern32
    CGGP$CorrName <- "Matern32"
  } else if (tolower(corr) %in% c("matern52", "m52", "m5")) {
    CGGP$CorrMat <- CGGP_internal_CorrMatMatern52
    CGGP$CorrName <- "Matern52"
  } else if (tolower(corr) %in% c("wendland0", "w0", "wend0", "triangle", "triangular", "tri", "triang")) {
    CGGP$CorrMat <- CGGP_internal_CorrMatWendland0
    CGGP$CorrName <- "Wendland0"
  } else if (tolower(corr) %in% c("wendland1", "w1", "wend1")) {
    CGGP$CorrMat <- CGGP_internal_CorrMatWendland1
    CGGP$CorrName <- "Wendland1"
  } else if (tolower(corr) %in% c("wendland2", "w2", "wend2")) {
    CGGP$CorrMat <- CGGP_internal_CorrMatWendland2
    CGGP$CorrName <- "Wendland2"
  } else {
    stop(paste0("corr given to CGGPcreate should be one of CauchySQT, CauchySQ,", 
                " Cauchy, Gaussian, PowerExponential, Matern32, or Matern52.\n",
                "Given value was ", corr, "."))
  }
  
  # Fix related parameters stored with CGGP
  CGGP$numpara <- CGGP$CorrMat(return_numpara=TRUE)
  CGGP$thetaMAP <- rep(0,CGGP$d*CGGP$numpara)
  CGGP$thetaPostSamples  <- matrix(2*rbeta(CGGP$d*CGGP$numpara*CGGP$numPostSamples,
                                           0.5, 0.5)-1,
                                   ncol=CGGP$numPostSamples )
  CGGP
}


#' Add more rows to CGGP object
#' 
#' After a lot of blocks have been added, the storage
#'
#' @param CGGP 
#' @param numrowstoadd 
#'
#' @return CGGP object
## @export
#' @noRd
CGGP_internal_addrows <- function(CGGP, numrowstoadd=20) {
  
  CGGP$uo <- rbind(CGGP$uo, matrix(0,numrowstoadd,ncol(CGGP$uo)))
  CGGP$ML <- nrow(CGGP$uo)
  
  # Need to get everything else upsized too
  CGGP$po = rbind(CGGP$po, matrix(0, nrow = 4 * numrowstoadd, ncol = ncol(CGGP$po))) #proposed levels tracker
  CGGP$pila = rbind(CGGP$pila, matrix(0, nrow = numrowstoadd, ncol=ncol(CGGP$pila))) #proposed immediate level ancestors
  CGGP$pala = rbind(CGGP$pala, matrix(0, nrow = numrowstoadd, ncol=ncol(CGGP$pala))) #proposedal all level ancestors
  CGGP$uala = rbind(CGGP$uala, matrix(0, nrow = numrowstoadd, ncol=ncol(CGGP$uala))) #used all level ancestors
  CGGP$pilaCOUNT = c(CGGP$pilaCOUNT, rep(0, numrowstoadd)) #count of number of pila
  CGGP$palaCOUNT = c(CGGP$palaCOUNT, rep(0, numrowstoadd)) #count of number of pala
  CGGP$ualaCOUNT = c(CGGP$ualaCOUNT, rep(0, numrowstoadd)) #count of number of uala
  CGGP$pogsize = c(CGGP$pogsize, rep(0, 4 * numrowstoadd))
  CGGP$w = c(CGGP$w, rep(0, numrowstoadd))
  
  CGGP
}

CGGP_internal_getdesignfromCGGP <- function(CGGP) {
  
  
  CGGP$gridsizes = matrix(CGGP$sizes[CGGP$uo[1:CGGP$uoCOUNT, ]], CGGP$uoCOUNT, CGGP$d)
  CGGP$gridsizest = matrix(CGGP$sizest[CGGP$uo[1:CGGP$uoCOUNT, ]], CGGP$uoCOUNT, CGGP$d)
  CGGP$gridsize = apply(CGGP$gridsizes, 1, prod)
  CGGP$gridsizet = apply(CGGP$gridsizest, 1, prod)
  CGGP$di = matrix(0, nrow = CGGP$uoCOUNT, ncol = max(CGGP$gridsize))
  CGGP$dit = matrix(0, nrow = CGGP$uoCOUNT, ncol = sum((CGGP$gridsize)))
  
  
  CGGP$design = matrix(0, nrow = sum(CGGP$gridsize), ncol = CGGP$d)
  CGGP$designindex = matrix(0, nrow = sum(CGGP$gridsize), ncol = CGGP$d)
  tv = 0
  for (blocklcv in 1:CGGP$uoCOUNT) {
    CGGP$di[blocklcv, 1:CGGP$gridsize[blocklcv]] = (tv + 1):(tv + CGGP$gridsize[blocklcv])
    for (dimlcv in 1:CGGP$d) {
      levelnow = CGGP$uo[blocklcv, dimlcv]
      if (levelnow < 1.5) {
        CGGP$design[(tv + 1):(tv + CGGP$gridsize[blocklcv]), dimlcv] = rep(CGGP$xb[1], CGGP$gridsize[blocklcv])
        CGGP$designindex[(tv + 1):(tv + CGGP$gridsize[blocklcv]), dimlcv] = rep(CGGP$xindex[1], CGGP$gridsize[blocklcv])
      } else{
        x0 = CGGP$xb[(CGGP$sizest[levelnow - 1] + 1):CGGP$sizest[levelnow]]
        xi0 = CGGP$xindex[(CGGP$sizest[levelnow - 1] + 1):CGGP$sizest[levelnow]]
        if (dimlcv < 1.5) {
          CGGP$design[(tv + 1):(tv + CGGP$gridsize[blocklcv]), dimlcv] = rep(x0, "each" = CGGP$gridsize[blocklcv] /
                                                                               CGGP$gridsizes[blocklcv, dimlcv])
          CGGP$designindex[(tv + 1):(tv + CGGP$gridsize[blocklcv]), dimlcv] = rep(xi0, "each" = CGGP$gridsize[blocklcv] /
                                                                                    CGGP$gridsizes[blocklcv, dimlcv])
          CGGP$blockassign[(tv + 1):(tv + CGGP$gridsize[blocklcv])] <- blocklcv
        }
        if (dimlcv > (CGGP$d - 0.5)) {
          CGGP$design[(tv + 1):(tv + CGGP$gridsize[blocklcv]), dimlcv] = rep(x0, CGGP$gridsize[blocklcv] /
                                                                               CGGP$gridsizes[blocklcv, dimlcv])
          CGGP$designindex[(tv + 1):(tv + CGGP$gridsize[blocklcv]), dimlcv] = rep(xi0, CGGP$gridsize[blocklcv] /
                                                                                    CGGP$gridsizes[blocklcv, dimlcv])
          CGGP$blockassign[(tv + 1):(tv + CGGP$gridsize[blocklcv])] <- blocklcv
        }
        if (dimlcv < (CGGP$d - 0.5)  && dimlcv > 1.5) {
          CGGP$design[(tv + 1):(tv + CGGP$gridsize[blocklcv]), dimlcv] = rep(rep(x0, each =
                                                                                   prod(CGGP$gridsizes[blocklcv, (dimlcv + 1):CGGP$d])), prod(CGGP$gridsizes[blocklcv, 1:(dimlcv - 1)]))
          CGGP$designindex[(tv + 1):(tv + CGGP$gridsize[blocklcv]), dimlcv] = rep(rep(xi0, each =
                                                                                        prod(CGGP$gridsizes[blocklcv, (dimlcv + 1):CGGP$d])), prod(CGGP$gridsizes[blocklcv, 1:(dimlcv - 1)]))
          CGGP$blockassign[(tv + 1):(tv + CGGP$gridsize[blocklcv])] <- blocklcv
        }
      }
    }
    
    tvv = 0
    if (blocklcv > 1.5) {
      for (ances in CGGP$uala[blocklcv, 1:CGGP$ualaCOUNT[blocklcv]]) {
        CGGP$dit[blocklcv, (tvv + 1):(tvv + CGGP$gridsize[ances])] = CGGP$di[ances, 1:CGGP$gridsize[ances]]
        tvv = tvv + CGGP$gridsize[ances]
      }
      CGGP$dit[blocklcv, (tvv + 1):(tvv + CGGP$gridsize[blocklcv])] = CGGP$di[blocklcv, 1:CGGP$gridsize[blocklcv]]
      Xset = CGGP$design[CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]], ]
      reorder = do.call(order, lapply(1:NCOL(Xset), function(kvt) Xset[, kvt]))
      CGGP$dit[blocklcv, 1:CGGP$gridsizet[blocklcv]] = CGGP$dit[blocklcv, reorder]
    } else{
      CGGP$dit[blocklcv, 1:CGGP$gridsize[blocklcv]] = CGGP$di[blocklcv, 1:CGGP$gridsize[blocklcv]]
    }
    
    tv = tv + CGGP$gridsize[blocklcv]
  }
  CGGP
}
