
#' Estimation of Person Parameters and calculation of Person Fit for the 1,2,3,4-PL model and the GPCM.
#'
#' PP-package has been developed to easily compute ML, WL (Warm 1989), MAP, EAP and robust estimates of person parameters for a given response matrix and given item parameters of the 1,2,3,4-PL model (Birnbaum 1968, Barton & Lord 1981) and the GPCM (Muraki 1992). It provides c++ routines which makes estimation of parameters very fast. Additional some methods to calculate the person fit are provided like the infit-outfit-, lz- and lz*-idex. Read the vignettes for getting started with this package.
#'
#'
#' @docType package
#' @name PP
#' @author Manuel Reif and Jan Steinfeld
#' @seealso \link{PPass}, \link{PP_gpcm}, \link{PP_4pl}, \link{PPall}, \link{Pfit}
#' @references 
#' 
#' Barton, M. A., & Lord, F. M. (1981). An Upper Asymptote for the Three-Parameter Logistic Item-Response Model.
#' 
#' Birnbaum, A. (1968). Some latent trait models and their use in inferring an examinee's ability. In Lord, F.M. & Novick, M.R. (Eds.), Statistical theories of mental test scores. Reading, MA: Addison-Wesley.
#' 
#' Drasgow, F., Levine, M. V. & Williams, E. A. (1985) Appropriateness measurement with polychotomous item response models and standardized indices. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{38(1)}, 67--86.
#'  
#' Muraki, Eiji (1992). A Generalized Partial Credit Model: Application of an EM Algorithm. Applied Psychological Measurement, 16, 159-176.
#'
#' Samejima, Fumiko (1993). An approximation of the bias function of the maximum likelihood estimate of a latent variable for the general case where the item responses are discrete. Psychometrika,  58, 119-138.
#' 
#' Snijders, T. B. (2001) Asymptotic null distribution of person fit statistics with estimated person parameter. \emph{Psychometrika}, \bold{66(3)}, 331--342. 
#'
#' Warm, Thomas A. (1989). Weighted Likelihood Estimation Of Ability In Item Response Theory. Psychometrika, 54, 427-450.
#' 
#' Wright, B. D. & Masters, G. N. (1990). Computation of OUTFIT and INFIT Statistics.  \emph{Rasch Measurement Transactions}, 3:4, 84-85.
#'
#' Yen, Y.-C., Ho, R.-G., Liao, W.-W., Chen, L.-J., & Kuo, C.-C. (2012). An empirical evaluation of the slip correction in the four parameter logistic models with computerized adaptive testing. Applied Psychological Measurement, 36, 75-87.
#'
#'@examples
#' set.seed(1522)
#' # intercepts
#' diffpar <- seq(-3,3,length=12)
#' # slope parameters
#' sl     <- round(runif(12,0.5,1.5),2)
#' la     <- round(runif(12,0,0.25),2)
#' ua     <- round(runif(12,0.8,1),2)
#'
#' # response matrix
#' awm <- matrix(sample(0:1,10*12,replace=TRUE),ncol=12)
#' # MLE estimation
#' res3plmle <- PP_4pl(respm = awm,thres = diffpar, slopes = sl,lowerA = la,type = "mle")
#' # WLE estimation
#' res3plwle <- PP_4pl(respm = awm,thres = diffpar, slopes = sl,lowerA = la,type = "wle")
#' # MAP estimation
#' res3plmap <- PP_4pl(respm = awm,thres = diffpar, slopes = sl,lowerA = la,type = "map")
#' 
#' # calculate person fit
#' res3plmlepfit <- Pfit(respm=awm,pp=res3plmle,fitindices=c("infit","outfit"))
#' 
#' # or estimate person parameter and calculate person fit in one step 
#' out <- PPass(respdf = data.frame(awm),thres = diffpar, items="all",
#'              mod=c("1PL"), fitindices= c("lz","lzstar","infit","outfit"))
NULL
