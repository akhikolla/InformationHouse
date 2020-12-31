#' A subset of data from Michie et al. (2009)
#'
#' The complete data consist of 101 studies reporting 122 interventions targeted at physical activity and healthy eating.
#' In this subset of the data, the interventions that include at least one
#' of the motivation-enhancing behaviour change techniques (BCTs) were selected (N = 106).
#'
#' @name  dat.BCT2009
#' @docType data
#' @details IMPORTANT: for questions about these data contact Xinru Li: x.li@math.leidenuniv.nl.
#' @references If you use these data, please refer to: Michie, S., Abraham, C., Whittington, C., McAteer, J., & Gupta, S. (2009). Effective techniques in healthy eating and physical activity interventions: a meta-regression. \emph{Health Psychology, 28(6)}, 690.
#'
#' @references An application of (a preliminary version of) meta-CART to this data set is given in: Dusseldorp, E., Van Genugten, L., van Buuren, S., Verheijden, M. W., & van Empelen, P. (2014). Combinations of techniques that effectively change health behavior: Evidence from Meta-CART analysis. \emph{Health Psychology, 33(12)}, 1530.
#'
#' @usage data(dat.BCT2009)
#'
#' @keywords data
#' @format A data frame of 106 interventions with five motivation-enhancing behavior change techniques (BCTs).
#' \itemize{
#'   \item study: The name of the intervention.
#'   \item g: The effect size of each intervention.
#'   \item vi: The sampling variance of the effect size.
#'   \item T1: Indicating whether the BCT1 "Provide information about
#'   behavior-health link" was used by the intervention. "1" for used
#'   and "0" for not used.
#'   \item T2: Indicating whether the BCT2 "Provide information on consequences"
#'    was used by the intervention. "1" for used
#'   and "0" for not used.
#'   \item T3: Indicating whether the BCT3 "Provide information about other's approval"
#'    was used by the intervention. "1" for used
#'   and "0" for not used.
#'   \item T4: Indicating whether the BCT4 "Prompt intention formation"
#'    was used by the intervention. "1" for used
#'   and "0" for not used.
#'   \item T25: Indicating whether the BCT25 " Motivational interviewing"
#'    was used by the intervention. "1" for used
#'   and "0" for not used.
#'
#'
#'   }
"dat.BCT2009"


#' A simulated meta-analytic data set with balanced pure interaction effects
#'
#' Data simulated from a true model with pure interactions between two
#' moderators: x1, x2. If x1 = 0 and x2 = 1 or x1 = 1 and x2 = 0, 
#' the true effect size is 0.50. Otherwise, the true effect size is 0.
#' @name  dat.balanced
#' @docType data
#'
#' @usage data(dat.balanced)
#'
#' @keywords simulated datasets
#' @format A data frame of 60 studies with 4 moderators
#' \itemize{
#'   \item efk: The effect size of each study expressed as Hedges' g
#'   \item vark: The sampling variance of the effect size
#'   \item x1 to x4: Four randomly generated moderators. x1, x2, and x4 are dichotomous variables,
#'   x3 is a continuous variable generated from uniform distribution.
#'   }
"dat.balanced"


#' A simulated meta-analytic data set
#'
#' Data simuated from a true model with a three-way interaction between three
#' moderators: m1, m2 and m3. If the values of the three moderators are all "B"s
#' the true effect size will be 0.80. Otherwise, the true effect size is 0.
#' @name  SimData
#' @docType data
#'
#' @usage data(SimData)
#'
#' @keywords simulated datasets
#' @format A data frame of 120 studies with 5 moderators
#' \itemize{
#'   \item efk: The effect size of each study expressed as Hedges' g
#'   \item vark: The sampling variance of the effect size
#'   \item m1 to m5: Five randomly generated moderators. m1 and m2 have two levels (A and B),
#'   whereas m3, m4 and m5 have three levels (A, B and C)
#'   }
"SimData"