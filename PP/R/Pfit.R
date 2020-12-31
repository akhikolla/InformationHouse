#' Person-Fit statistics
#' 
#' Compute several person fit statistic for the 1-PL, 2-PL, 3-PL, 4-PL and PCM. 
#'
#' Please note that currently only the likelihood based LZ-Index (Drasgow, Levine, and Williams, 1985) and LZ*-Index (Snijders, 2001) are implemented. Also the INFIT-OUTIFT (Wright and Masters, 1982, 1990) statistic as well as the polytomouse version of INFIT-OUTFIT are supported. Other person fit statistics will be added soon.
#' 
#' The calculation of the person fit statistics requires the numeric response-matrix as well as an object of the fourpl-class. So first you should estimate the person parameter and afterwards calculate the person fit statistics. You could also use our PPass-function to estimate the person parameter and calculate the desired person fit simultaneously.
#' It is possible to calculate several person fit statistics at once, you only have to specify them in a vector.
#' 
#' For the Partial Credit model we currently support the infit-outfit statistic. Please submit also the numeric response-matrix as well as the estimated person parameter with an gpcm-class.
#' 
#'@param respm	      numeric response matrix
#'@param pp 		      object of the class fourpl with estimated person parameter
#'@param fitindices		character vector of desired person fit statistics c("lz","lzstar","infit","outfit")
#'@param SE		        logical: if true standard errors are computed using jackknife method
#'
#' @return list of person-fits for each person-fit statistic
#' 
#' \itemize{
#' \item the list of person-fits contains the calculated person-fit (like lz, lzstar) and also additional information like p-value or standard error if desired.
#' \item the additional information is provided after the short form of the personfit
#' \item lz (lz)
#' \item lzstar (lzs)
#' \item infit the mean-square statistic (in)
#' \item outfit the mean-square statistic (ou)
#' \item _unst: unstandardised
#' \item _se: standard error
#' \item _t: t-value
#' \item _chisq: $chi^2$-value
#' \item _df: defrees of freedom
#' \item _pv: p-value
#' }
#' 
#' @rdname pfit
#' @seealso \link{PPall}, \link{PP_4pl}, \link{PPass}
#'
#'
#'@author Jan Steinfeld
#'@references 
#'
#' Armstrong, R. D., Stoumbos, Z. G., Kung, M. T. & Shi, M. (2007). On the performance of the lz person-fit statistic.  \emph{Practical Assessment, Research & Evaluation}, \bold{12(16)}. Chicago	
#' 
#' De La Torre, J., & Deng, W. (2008). Improving Person-Fit Assessment by Correcting the Ability Estimate and Its Reference Distribution. Journal of Educational Measurement, \bold{45(2)}, 159-177.
#' 
#' Drasgow, F., Levine, M. V. & Williams, E. A. (1985) Appropriateness measurement with polychotomous item response models and standardized indices. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{38(1)}, 67--86.
#' 
#' Efron, B., & Stein, C. (1981). The jackknife estimate of variance. \emph{The Annals of Statistics}, \bold{9(3)}, 586-596.	
#' 
#' Karabatsos, G. (2003) Comparing the Aberrant Response Detection Performance of Thirty-Six Person-Fit Statistics. \emph{Applied Measurement In Education}, \bold{16(4)}, 277--298.
#' 
#' Magis, D., Raiche, G. & Beland, S. (2012) A didactic presentation of Snijders's l[sub]z[/sub] index of person fit with emphasis on response model selection and ability estimation. \emph{Journal of Educational and Behavioral Statistics}, \bold{37(1)}, 57--81.
#' 
#' Meijer, R. R. & Sijtsma, K. (2001) Methodology review: Evaluating person fit. \emph{Applied Psychological Measurement}, \bold{25(2)}, 107--135.
#' 
#' Molenaar, I. W. & Hoijtink, H. (1990) The many null distributions of person fit indices. \emph{Psychometrika}, \bold{55(1)}, 75--106. 
#' 
#' Mousavi, A. & Cui, Y. Evaluate the performance of and of person fit: A simulation study.
#' 
#' Reise, S. P. (1990). A comparison of item-and person-fit methods of assessing model-data fit in IRT.  \emph{Applied Psychological Measurement}, \bold{14(2)}, 127-137.
#' 
#' Snijders, T. B. (2001) Asymptotic null distribution of person fit statistics with estimated person parameter. \emph{Psychometrika}, \bold{66(3)}, 331--342. 
#' 
#' Wright, B. D. & Masters, G. N. (1990). Computation of OUTFIT and INFIT Statistics.  \emph{Rasch Measurement Transactions}, 3:4, 84-85.
#' 
#' Wright, B. D., & Masters, G. N. (1982). \emph{Rating Scale Analysis. Rasch Measurement.} MESA Press, 5835 S. Kimbark Avenue, Chicago, IL 60637.
#' 
#' @example ./R/.examples_pfit.R
#' @keywords Person fit, LZ-Index, Infit-Outfit
#' @export
Pfit <- function(respm,pp,fitindices,SE=FALSE) UseMethod("Pfit",object=pp)
# ------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@method Pfit fourpl
#'@export
  Pfit.fourpl <- function(respm, pp, fitindices=c("lz","lzstar","infit","outfit"),SE=FALSE){

    if(any(pp$type%in%c("eap","robust"))) stop("Only 'mle','wle' and 'map' ability estimates are supported \n")

    pfitfunctions <- list("lz" = lz,
                          "lzstar" = lzstar,
                          "infit" = Infit,
                          "outfit" = Outfit
    )
    
    fitindices <- match.arg(fitindices, several.ok = TRUE)  

    pfitfunctions_red <- pfitfunctions[names(pfitfunctions)%in%fitindices]
    
    args <- list(list("data"=respm, 
         "thetas"=pp$resPP$resPP[,"estimate"], 
         "betas"=pp$ipar$thres[2,], 
         "lowerAs"=pp$ipar$lowerA, 
         "slopes"=pp$ipar$slopes, 
         "higherAs"=pp$ipar$upperA,
         "method"=pp$type, 
         "mu"=pp$ipar$mu, 
         "sigma"=sqrt(pp$ipar$sigma2)
    ))
    
    out <- mapply(function(x,y) do.call("y",x), x=args, y=pfitfunctions_red,SIMPLIFY = FALSE)
    names(out) <- names(pfitfunctions_red)
    
    if(SE){
      if(any(fitindices=="lz")){
        lz.se <- jackknife(data = respm, pp=pp, fit="lz")
        out$lz <- cbind(out$lz,"lz_se"=lz.se)
      }
      if(any(fitindices=="lzstar")){
        lzstar.se <- jackknife(data = respm, pp=pp, fit="lzstar")
        out$lzstar <- cbind(out$lzstar,"lzs_se"=lzstar.se)
      }
    }
    
    class(out) <- append(class(out),"PPfit")
    return(out)
  }
  
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  #'@rdname pfit
  #' 
  #' @method Pfit gpcm
  #' @export
  Pfit.gpcm <- function(respm, pp, fitindices=c("infit","outfit"),SE=FALSE){
    if(SE){warning("There is no jacknife currently supported. \n")}
    if(any(pp$type%in%c("map","eap","robust"))) stop("Only 'mle' and 'wle' ability estimates are supported \n")
    
    if(!all(fitindices%in%c("infit","outfit"))){ warning("Only 'infit and outfit' are currently supported. The calculation is executed with infit outfit \n"); fitindices <- c("infit","outfit")}
    if(any(pp$ipar$slopes>1)) warning("Currently only the PCM-Modell is supported \n")
    
    pfitfunctions <- list("infit" = Infitpoly,
                          "outfit" = Outfitpoly)

     fitindices <- match.arg(fitindices, several.ok = TRUE)  

    pfitfunctions_red <- pfitfunctions[names(pfitfunctions)%in%fitindices]
    
    args <- list(list("data"=respm, 
                      "thetas"=pp$resPP$resPP[,"estimate"], 
                      "thresholds"=pp$ipar$thres, 
                      "slopes"=NULL
    ))
    
    out <- mapply(function(x,y) do.call("y",x), x=args, y=pfitfunctions_red,SIMPLIFY = FALSE)
    names(out) <- names(pfitfunctions_red)
    class(out) <- append(class(out),"PPfit")
    return(out)
  }
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
  
 # #' @rdname pfit
#  #' 
#  #' @method Pfit gpcm4pl
#  #' @export
#  Pfit.gpcm4pl <- function(respm, pp, fitindices){
#      cat("the mixed method for person fits is not yet implemented \n")
#  }
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------
