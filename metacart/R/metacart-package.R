#' R package for meta-CART
#'
#'In meta-analysis, heterogeneity often exists between studies.
#'To understand this heterogeneity, researchers search for study 
#'characteristics (i.e., potential moderators) that may account for 
#'the variance in study effect sizes. When multiple potential moderators 
#'are available (e.g., intervention characteristics), 
#'traditional meta-analysis methods often lack sufficient power to investigate
#'interaction effects between moderators, especially high-order interactions.
#'To solve this problem, meta-CART was proposed by integrating Classification and Regression Trees (CART)
#'into meta-analysis. The method idenfities the interaction effects between influential moderators,
#'partitions the studies into more homogeneous subgroups, and estimates summary effect size in each subgroup.
#'The fixed effect or random effects assumption can be consistently taken into account in both tree-growing process and subgroup analysis.
#'
#' @details \tabular{ll}{
#'   Package: \tab metacart\cr
#'   Type: \tab Package\cr
#'   Version: \tab 2.0-0\cr
#'   Date: \tab 2018-11-12\cr
#'   License: \tab GPL\cr
#' }
#'
#'   This method is suitable for identifying interaction effects between dichotomous,
#'   ordinal, continuous, and nominal moderators.
#'   The output of a \code{REmrt} object shows meta-CART analysis results based on the random effects model.
#'   And the output of a \code{FEmrt} object shows meta-CART analysis results based on the fixed effect model.
#'   The two objects display results for subgroup analysis including the Q-statistic and estimates for the subgroup effect sizes.
#'   Furthermore, the predict functions \code{predict.REmrt} and \code{predict.FEmrt} can be used to predict the effect size given the moderators.
#'   The plot functions \code{plot.REmrt} and \code{plot.FEmrt} show the interaction effects between identified moderators.
#'
#'   The core functions of the package are \code{\link{FEmrt}} and  \code{\link{REmrt}}.
#'
#' @author Maintainer: Xinru Li <x.li@math.leidenuniv.nl>; Contributors: Elise Dusseldorp, Kaihua Liu (supported with the plot function), Jacqueline Meulman.
#' @references Dusseldorp, E., van Genugten, L., van Buuren, S., Verheijden, M. W., & van Empelen, P. (2014). Combinations of techniques that effectively change health behavior: Evidence from meta-cart analysis.  \emph{Health Psychology, 33(12)}, 1530-1540. doi:
#'      10.1037/hea0000018.
#' @references Li, X., Dusseldorp, E., & Meulman, J. J. (2017). Meta-CART: A tool to identify interactions between moderators in meta-analysis. \emph{British Journal of Mathematical and Statistical Psychology, 70(1)}, 118-136. doi: 10.1111/bmsp.12088.
#'
#' @references Therneau, T., Atkinson, B., & Ripley, B. (2014) rpart: Recursive partitioning and regression trees. R package version, 4-1.
#' @references The articles of our own work can be found at \url{http://www.elisedusseldorp.nl/}
#' @keywords package
#' @seealso \code{\link{FEmrt}}, \code{\link{REmrt}}, \code{\link{summary.FEmrt}},\code{\link{summary.REmrt}},
#'   \code{\link{plot.FEmrt}},\code{\link{plot.REmrt}},\code{\link{predict.FEmrt}},\code{\link{predict.REmrt}}
#'
#' @docType package
#' @name metacart-package
NULL
