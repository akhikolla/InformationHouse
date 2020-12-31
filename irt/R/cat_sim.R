###############################################################################@
########################### create_cat_design ##################################
###############################################################################@
#' Computerized Adaptive Test (CAT) Simulation Design
#' @description
#' \code{create_cat_design} is a helper function for
#' \code{\link{cat_sim}} and \code{\link{cat_sim_fast}} functions. It
#' defines the simulation design.
#'
#' Ideally, there is a design element for each item. So within this design
#' (which is a list), there are $k$ design elements for each potentially
#' administered item. Each of these sub-design elements are also a list.
#'
#' @param ip An \code{\link{Itempool-class}} object containing item parameters,
#'          content information, etc.
#'
#'          If \code{ip = NULL} this means this is an infinite item pool,
#'          where b is on demand, c = 0 and a = 1, D = 1.7.
#'
#'          If \code{true_ip} argument is \code{NULL}, this item pool will
#'          be used to generate item responses.
#' @param title A string value representing the title of this CAT design.
#' @param true_ip An \code{\link{Itempool-class}} object which holds the true
#'          values of item pool parameters that will be used to generate item
#'          responses. This is an optional argument. If it is \code{NULL}
#'          and \code{ip} is not missing, then, item responses will be
#'          generated using \code{ip}.
#'
#'          \strong{Default}: \code{NULL}
#' @param first_item_rule The method how the first item is administered.
#'          The main effect of this is to select the first item administered
#'          to an examinee. If, for example, first item is desired to be a
#'          fixed one or randomly selected from the item pool, then set that
#'          rule in \code{next_item_rule}.
#'
#'          \strong{Default}: \code{'fixed_theta'}
#'
#'          Possible values and required parameters:
#'          \describe{
#'            \item{\strong{NULL}}{If no separate first item selection rule is
#'              necessary, the first item will be selected using the
#'              \code{next_item_rule} and it's parameters \code{next_item_par}.
#'            }
#'            \item{\strong{"fixed_theta"}}{Fixed starting value.
#'
#'              Required parameters for \code{first_item_par} argument if
#'              this rule is selected:
#'              \describe{
#'                \item{theta}{The value of the initial theta estimate.}
#'              }
#'            }
#'            \item{\strong{"theta_range"}}{An initial theta estimate within
#'              \code{min_theta} and \code{max_theta} will be randomly selected.
#'
#'              Required parameters for \code{first_item_par} argument if
#'              this rule is selected:
#'              \describe{
#'                \item{min_theta}{Minimum theta value of the interval.}
#'                \item{max_theta}{Maximum theta value of the interval.}
#'              }
#'            }
#'          }
#'
#' @param first_item_par Parameters for the first item rule.
#'
#'          \strong{Default}: \code{list(theta = 0)}
#' @param next_item_rule A vector of length one or length maximum test length
#'          which is designating the next item selection rules.
#'
#'          \strong{Default}: \code{'mfi'}
#'
#'          Note that, currently, if there are testlets in an item pool and a
#'          testlet is selected for administration using one of the methods
#'          below, all items within that testlet will be administered regardless
#'          of the next item selection rule.
#'
#'          Possible values and required parameters:
#'          \describe{
#'            \item{\strong{random}}{
#'              Randomly select items from the item pool.
#'              Exposure control rules and parameters will be ignored for this
#'              selection rule.
#'
#'              Required parameters: None.
#'            }
#'            \item{\strong{mfi}}{
#'              Maximum Fisher Information.
#'
#'              Required parameters: None.
#'            }
#'            \item{\strong{mepv}}{
#'              Minimum Expected Posterior Variance.
#'
#'              Required Parameters:
#'              \describe{
#'                \item{"var_calc_method"}{
#'                  Which method to use to calculate the posterior variance.
#'                  See Equation (4) of Choi and Swartz (2009), Comparison of
#'                  CAT Criteria for Polytomous Items.
#'
#'                  Available options are:
#'
#'                  \describe{
#'                    \item{\code{"eap"}}{
#'                      Use the variance from expected a posteriori estimation.
#'                    }
#'                    \item{\code{"owen"}}{
#'                      Use the variance from Owen's Bayesian estimation.
#'                      For \code{"Rasch"}, \code{"1PL"}, \code{"2PL"},
#'                      \code{"3PL"} models this is much faster than
#'                      \code{"eap"} option above.
#'                    }
#'                  }
#'                }
#'              }
#'            }
#'            \item{\strong{b_optimal}}{
#'              Select item which has item difficulty that is close to the
#'              current ability estimate.
#'
#'              Required parameters: None.
#'            }
#'            \item{\strong{fixed}}{
#'              Administer a fixed set of items from the item pool. This is
#'              basically a linear fixed length test where the order of items
#'              are predefined. Exposure control rules and parameters will be
#'              ignored for this selection rule.
#'
#'              Required Parameters:
#'              \describe{
#'                \item{item_id}{
#'                  A vector of the item IDs that should be administered.
#'                }
#'              }
#'            }
#'          }
#' @param next_item_par A list of length one or length maximum test length
#'          that sets the parameters of next item selection rules. It can also
#'          be \code{NULL}, in which case no parameters necessary for that
#'          next item selection procedure.
#'
#'          \strong{Default}: \code{NULL}
#' @param ability_est_rule A vector of length one or length maximum test length
#'          which is designating the next item selection rules.
#'
#'          \strong{Default}: \code{"eap"}
#'
#'          Possible values and required parameters:
#'          \describe{
#'            \item{\strong{"eap"}}{
#'            Expected-a-posteriori.
#'              Required parameters:
#'              \describe{
#'                \item{prior_dist}{
#'                  Distribution of the prior distribution.
#'                  Available values:
#'
#'                  * \code{norm} for normal distribution,
#'                  * \code{unif} for uniform distribution.
#'
#'                  The default value is \code{norm}.
#'                }
#'                \item{prior_par}{
#'                  A vector of prior parameters.
#'
#'                  * For normal distribution \code{c(0, 1)}, see \code{?dnorm}
#'                  * For uniform distribution \code{c(-3, 3)}, see \code{?dunif}
#'
#'                  The default value is \code{c(0, 1)}.
#'                }
#'                \item{min_theta}{
#'                  Minimum possible value of theta. It is a lower bound.
#'
#'                  The default value is \code{-4}.
#'                }
#'                \item{max_theta}{
#'                  Maximum possible value of theta. It is an upper bound.
#'
#'                  The default value is \code{4}.
#'                }
#'                \item{no_of_quadrature}{
#'                  The number of quadrature, more specifically the number of
#'                  bins the theta range should be divided. The more bins, the
#'                  more precise (and slower) the estimates will be.
#'
#'                  The default value is \code{50}.
#'                }
#'              }
#'            }
#'            \item{\strong{"owen"}}{
#'            Owen's Bayesian Estimation
#'              Required parameters:
#'              \describe{
#'                \item{prior_mean}{Prior mean value. The default value is
#'                  \code{0}.}
#'                \item{prior_var}{Prior variance value.The default value is
#'                  \code{1}.}
#'              }
#'            }
#'            \item{\strong{"ml"}}{
#'              Maximum likelihood estimation using Newton-Raphson algorithm.
#'              If this method is used, the standard error of ability estimates
#'              are calculated using the inverse information value at this
#'              theta estimate.
#'
#'              Required parameters:
#'              \describe{
#'                \item{min_theta}{Minimum possible value of theta. It is a
#'                  lower bound. The default value is -4.
#'                }
#'                \item{max_theta}{Maximum possible value of theta. It is an
#'                  upper bound. The default value is 4.
#'                }
#'                \item{criterion}{This value determines the accuracy of
#'                  estimates. Smaller values lead more accuracy but the
#'                  speed of estimation reduces as the value of \code{criterion}
#'                  decreases. The default value is 0.001.
#'                }
#'              }
#'            }
#'            \item{\strong{"eap_ml"}}{
#'              Expected-a-posteriori until an imperfect item response string,
#'              then switch to Maximum Likelihood estimation.
#'              Required parameters:
#'              \describe{
#'                \item{prior_dist}{
#'                  Distribution of the prior distribution.
#'
#'                  Available values:
#'
#'                  \code{norm} for normal distribution,
#'
#'                  \code{unif} for uniform distribution.
#'                }
#'                \item{prior_par}{
#'                  A vector of prior parameters.
#'                  For normal distribution \code{c(0, 1)}, see \code{?dnorm}
#'                  For uniform distribution \code{c(-3, 3)}, see \code{?dunif}
#'                }
#'                \item{min_theta}{
#'                  Minimum possible value of theta. It is a lower bound.
#'                }
#'                \item{max_theta}{
#'                  Maximum possible value of theta. It is an upper bound.
#'                }
#'                \item{no_of_quadrature}{
#'                  The number of quadrature, more specifically the number of
#'                  bins the theta range should be divided. The more bins, the
#'                  more precise (and slower) the estimates will be.
#'                }
#'              }
#'            }
#'            \item{\strong{"sum_score"}}{
#'              Simple sum score.
#'              Required parameters: \code{NULL}
#'            }
#'          }
#' @param ability_est_par A list of length one or length maximum test length
#'          that sets the parameters of ability estimation rules. It can also
#'          be \code{NULL}.
#'
#'          * If \code{ability_est_rule = "eap"} then the default is
#'          \code{list(prior_dist = "norm", prior_par = list(mean = 0, sd = 2),
#'                     min_theta = -4, max_theta = 4)}
#'          * If \code{ability_est_rule = "owen"} then the default is
#'          \code{list(prior_mean = 0, prior_var = 1)}
#'
#'          If it is \code{NULL}, either no parameters necessary for that
#'          ability estimation rule or the defaults of that ability selection
#'          rule will be selected.
#'
#'          If it is a list of one, it means that the parameters will be the
#'          same throughout the test. The names of the list elements will
#'          represent the parameter types.
#'
#'          A list of lists with length of maximum test length designate
#'          different parameters for different items in the test progress.
#'
#' @param final_ability_est_rule The ability estimation method that will be
#'          used to calculate the final ability estimate. The methods and
#'          the parameters are the same as \code{ability_est_rule} and
#'          \code{ability_est_par}. Please see those for details.
#'
#'          \strong{Default}: \code{NULL}
#' @param final_ability_est_par A list of parameters that will be used
#'          for the method designated by the \code{final_ability_est_rule}.
#'
#'          \strong{Default}: \code{NULL}
#' @param termination_rule This parameter determines how CAT algorithm decides
#'          terminate the test.
#'
#'          The order of termination rules is important. The algorithm will
#'          check the rules in that order. If for example
#'          \code{termination_rule = c('min_se', 'max_item')}, first whether
#'          the SE smaller than a certain value checked and if it is smaller,
#'          then even the maximum number of items haven't been administered,
#'          test will terminate.
#'
#'          The \code{"min_item"} and \code{"max_item"} has a special property
#'          where, for \code{"min_item"}, if the number of items administered
#'          smaller than \code{min_item}, then test will not terminate
#'          regardless of whether other rules satisfied. Similarly, for
#'          \code{"max_item"}, if the number of items is larger than
#'          \code{max_item}, the test will terminate regardless of whether other
#'          conditions satisfied or not. If both \code{"min_item"} and
#'          \code{"max_item"} are in termination rules, then, test will end when
#'          both conditions satisfied, i.e. when the number of items
#'          administered is equal to or larger than \code{max_item} value in
#'          \code{termination_par}.
#'
#'          The "test length" refers to "Item" objects, i.e. individual items
#'          not testlets. For example, if an item pool has 10 testlets each
#'          having 2 items and 15 standalone items which are not within a
#'          testlet, then the test length can go up to 35 (2 x 10 + 15).
#'
#'          \strong{Default}: \code{c("min_item", "min_se", "max_item")}
#'
#'          \code{"termination_rule"} should be a vector that composed of the
#'          following termination rules:
#'
#'          \describe{
#'            \item{\code{"min_item"}}{The minimum number of items should be
#'                  satisfied. If the number of administered items are equal to
#'                  or larger than this number test ends. }
#'            \item{\code{"max_item"}}{The maximum number of items should not be
#'                  exceeded.}. If this is missing, then the item pool
#'                  size will be set as maximum length.
#'            \item{\code{"min_se"}}{If the standard error exceeds \code{min_se}
#'                  value, then the test will terminate.}
#'            \item{\code{"sprt"}}{Sequential Probability Ratio Test (SPRT).
#'            SPRT tests two hypotheses:
#'
#'            \eqn{H_0}: Examinee's ability \eqn{\hat \theta = \theta_0}
#'
#'            \eqn{H_1}: Examinee's ability \eqn{\hat \theta = \theta_1}
#'
#'            After the administration of each item, the likelihood (or
#'            log-likelihood) of the response string is calculated at
#'            \eqn{\theta_0} and \eqn{\theta_1}. The ratio of this likelihood is
#'            then compared to two decision points, \eqn{A} and \eqn{B}.
#'
#'            \deqn{LR = \frac{L(\theta = theta_1)}{\theta = theta_0}}
#'
#'            In order to calculate the lower (\eqn{A}) and upper (\eqn{B})
#'            decision points, one needs to set \eqn{\alpha} and \eqn{\beta}.
#'            \eqn{\alpha} represents the rate of false positive classification
#'            errors \eqn{(0 < \alpha < 1)}, i.e. examinees whose true
#'            classification is fail but passed at the end of test. \eqn{\beta}
#'            is the rate of false negative classification errors \eqn{(0 <
#'            \beta < 1)}, i.e. examinees whose true classification is pass but
#'            failed at the end of test. \eqn{A} and \eqn{B} can be calculated
#'            as:
#'
#'            \deqn{A = \frac{1 - \beta}{\alpha}}
#'
#'            \deqn{B = \frac{\beta}{1 - \alpha}}
#'
#'            If \eqn{LR > A}, examinee passes the test and if \eqn{LR < B}
#'            examinee fails the test. If \eqn{B < LR < A}, test continues
#'            until the maximum number of items reached (or some other test
#'            termination criteria satisfied.)
#'
#'            \code{"sprt"} termination rule needs \code{termination_par}, where
#'            the following parameters should be given in a list:
#'            \describe{
#'              \item{\code{"theta_0"}}{The highest theta value that the
#'                test developer is willing to fail an examinee. }
#'              \item{\code{"theta_1"}}{The lowest theta value that the
#'                test developer is willing to pass an examinee.}
#'              \item{\code{"alpha"}}{The rate of false positive classification
#'                errors (0 < \code{alpha} < 1), i.e. examinees whose true
#'                classification is fail but passed at the end of test.}
#'              \item{\code{"beta"}}{The rate of false negative classification
#'                errors (0 < \code{beta} < 1), i.e. examinees whose true
#'                classification is pass but failed at the end of test.}
#'              }
#'            Example: \code{termination_par = list(sprt = list(theta_0 = -.9,
#'                                                              theta_1 = -.1,
#'                                                              alpha = 0.05,
#'                                                              beta = 0.05))}
#'            }
#'          }
#' @param termination_par A list of termination rule parameters. This
#'          is a named list with length equal to the length of
#'          \code{termination_rule} argument. The names of the list elements
#'          should correspond to the elements of \code{termination_rule}
#'          argument.
#'
#'          \strong{Default}: \code{list(min_item = 10, min_se = 0.33,
#'                                       max_item = 20)}
#' @param exposure_control_rule A vector of length one or length maximum test
#'          length which is designating the next item selection rules. It can
#'          be \code{NULL} in which case there won't be any exposure control.
#'
#'          \strong{Default}: \code{NULL}, No exposure control will be imposed
#'            on item selection.
#'
#'          Possible values and required parameters:
#'          \describe{
#'            \item{\code{NULL}}{No exposure control.}
#'            \item{"randomesque"}{
#'              Select one of the most informative \code{num_items} items.
#'              \describe{
#'                \item{\code{num_items}}{The number of items to select from.}
#'              }
#'            }
#'            \item{\code{"sympson-hetter"}}{
#'              The algorithm of Sympson-Hetter exposure control is explained in
#'              Sympson and Hetter (1985).
#'
#'              This method does not require any additional
#'              "exposure_control_par" but each item/testlet should have
#'              a "misc" slot like the following
#'              \code{misc = list(sympson_hetter_k = .75)}.
#'
#'              When using 'sympson-hetter' exposure control rule, please ensure
#'              that there are sufficient number of items with
#'              'sympson_hetter_k' values 1. Otherwise, examinees might not
#'              get a complete test and an error might be raised by the
#'              simulation function.
#'            }
#'          }
#' @param exposure_control_par A list of length one or maximum test length
#'          designating the exposure control for each item. If there are no
#'          parameters it will be \code{NULL}.
#'
#'          \strong{Default}: \code{NULL}
#' @param content_bal_rule Whether a content balancing is imposed on item
#'          selection. Default value is \code{NULL}, where no content balancing
#'          will be imposed on item selection.
#'
#'          \strong{Default}: \code{NULL}
#'
#'          Possible values and required parameters:
#'          \describe{
#'            \item{\code{NULL}}{No content balancing.}
#'            \item{\strong{max_discrepancy}}{Given a target content
#'              distribution, the content with maximum discrepancy with target
#'              discrepancy will be administered.
#'
#'              Required parameters:
#'              \describe{
#'                \item{target_dist}{Target content ratios.
#'                  For example, suppose there are three content areas:
#'                  Geometry, Algebra and Arithmetic. If the plan for the test
#'                  is to include 30% Geometry items, 50% Algebra items and 20%
#'                  Arithmetic items, then, the \code{target_dist} should be:
#'                  c(Geometry = .3, Arithmetic = .2, Algebra = .5). The names
#'                  in the vector should correspond to the names of the content
#'                  areas in the item pool. \code{target_dist} should include
#'                  each content area within the item pool for it to work
#'                  properly. If the sum of the \code{target_dist} is larger
#'                  than 1, it will be converted to ratios.
#'                }
#'              }
#'            }
#'          }
#' @param content_bal_par Parameters of \code{content_bal_rule}. A list, a
#'          list of lists or \code{NULL}.
#'
#'          \strong{Default}: \code{NULL}
#'
#' @param ability_type The type of ability the test is measuring. By default
#'          it is IRT based single 'theta'.
#'          \describe{
#'            \item{\code{"theta"}}{Theta for unidimensional IRT models}
#'            \item{\code{"multi_theta"}}{Theta vector for multidimensional IRT
#'              models (Not Implemented Yet).}
#'            \item{\code{"cdm"}}{An attribute vector (Not Implemented Yet).}
#'            \item{\code{"raw_score"}}{Raw score (i.e. total score) of an
#'              examinee.}
#'            }
#'
#'          \strong{Default}: \code{"theta"}
#' @return A \code{cat_design} object that holds the test specifications of a
#'         CAT.
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @seealso \code{\link{cat_sim}}
#'
#' @references
#'   Sympson, J., & Hetter, R. D. (1985). Controlling item-exposure rates in
#'   computerized adaptive testing. 973–977.
#'
#'
#' @examples
#' ### Example Designs ###
#' # Fixed length test IRT test with ability estimation EAP-ML
#' n_items <- 30
#' ip <- itempool(data.frame(a = runif(n_items, .5, 1.5), b = rnorm(n_items)))
#' cd <- create_cat_design(ip = ip, next_item_rule = 'random',
#'                         termination_rule = 'min_item',
#'                         termination_par = list('min_item' = n_items))
#' cd
#' create_cat_design(ip = ip, next_item_rule = 'random')
#'
#'
#' n_ip <- 55
#' ip <- itempool(data.frame(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip)))
#' # Check the default:
#' create_cat_design()
#' create_cat_design(ip = ip)
#'
#' ### Termination Rule ###
#' create_cat_design(
#'   termination_rule = c('min_item', 'min_se', 'max_item'),
#'   termination_par = list(min_item = 10, min_se = .33, max_item = 20))
#'
#' cd <- create_cat_design(ip = ip, termination_rule = c('min_item', 'min_se'),
#'                         termination_par = list(min_item = 10, min_se = .33))
#'
#' ### Next Item Rule ###
#' create_cat_design(ip = ip, next_item_rule = 'random', next_item_par = NULL)
#' create_cat_design(
#'   ip = ip, termination_rule = c('min_item', 'max_item'),
#'   termination_par = list(min_item = 20, max_item = 20),
#'   next_item_rule = 'fixed',
#'   next_item_par = list(item_id = ip$id[1:20]))
#'
#' # Linear test where all of the items in the item pool administered in the
#' # same order as item pool
#' ip <- generate_ip(n = 15)
#' create_cat_design(
#'   ip = ip, termination_rule = c('max_item'),
#'   termination_par = list(max_item = 15),
#'   next_item_rule = 'fixed')
#'
#' # Generate an item pool with two testlets and three standalone items and
#' # administer first seven items as a linear test.
#' ip <- c(generate_testlet(n = 2, id = "t1"), generate_ip(n = 3),
#'         generate_testlet(n = 5, id = "t2"))
#' create_cat_design(
#'   ip = ip, termination_rule = c('max_item'),
#'   termination_par = list(max_item = 7),
#'   next_item_rule = 'fixed')
#'
#'
#' # A linear test where the item order is predefined.
#' ip1 <- itempool(data.frame(b = rnorm(5)), id = paste0("i",1:5))
#' cd <- create_cat_design(
#'   ip = ip1,
#'   next_item_rule = 'fixed',
#'   next_item_par = list(item_id = c("i3", "i2", "i4", "i5", "i1")),
#'   ability_est_rule = "eap",
#'   termination_rule = 'max_item', termination_par = list(max_item = 5))
#'
#' ### Ability Estimation Rule ###
#' create_cat_design(
#'   ability_est_rule = 'eap',
#'   ability_est_par = list(prior_dist = 'unif',
#'                          prior_par = list(min = -2, max = 2),
#'                          min_theta = -4, max_theta = 4,
#'                          no_of_quadrature = 31))
#' create_cat_design(
#'   ability_est_rule = 'ml',
#'   ability_est_par = list(min_theta = -4, max_theta = 4, criterion = 0.01))
#'
#' ### Exposure Control ###
#' create_cat_design(exposure_control_rule = 'randomesque',
#'                 exposure_control_par = list(num_items = 1))
#'
#' # 5-4-3-2-1 exposure control
#' create_cat_design(
#'   exposure_control_rule = 'randomesque',
#'   exposure_control_par = lapply(c(5:1, rep(1, 15)),
#'                                 function(x) list(num_items = x)))
#'
#' ### Content Balancing ###
#' create_cat_design(
#'   content_bal_rule = 'max_discrepancy',
#'   content_bal_par = list(target_dist = c(
#'     Geometry = .3, `Rational Numbers` = .2, Algebra = .5)))

create_cat_design <- function(
  ip = NULL,
  title = NULL,
  true_ip = NULL,
  first_item_rule = "fixed_theta",
  first_item_par = list(theta = 0),
  next_item_rule = "mfi",
  next_item_par = NULL,
  ability_est_rule = "eap",
  ability_est_par = NULL,
  final_ability_est_rule = NULL,
  final_ability_est_par = NULL,
  termination_rule = c("min_item", "min_se", "max_item"),
  termination_par = list(min_item = 10, min_se = .33, max_item = 20),
  exposure_control_rule = NULL,
  exposure_control_par = NULL,
  content_bal_rule = NULL,
  content_bal_par = NULL,
  ability_type = "theta"
  ) {

  # ip = NULL
  # title = NULL
  # true_ip = NULL
  # first_item_rule = "fixed_theta"
  # first_item_par = list(theta = 0)
  # next_item_rule = "mfi"
  # next_item_par = NULL
  # ability_est_rule = "eap"
  # ability_est_par = NULL
  # final_ability_est_rule = NULL
  # final_ability_est_par = NULL
  # termination_rule = c("min_item", "min_se", "max_item")
  # termination_par = list(min_item = 10, min_se = .33, max_item = 20)
  # exposure_control_rule = NULL
  # exposure_control_par = NULL
  # content_bal_rule = NULL
  # content_bal_par = NULL
  # ability_type = "theta"

                        #####################@###
                        ### Implemented Rules ###
                        #####################@###

  first_item_rules <- list(
    fixed_theta = list(par_names = c("theta")),
    theta_range = list(par_names = c("min_theta", "max_theta"))
    )
  # next_item_rules <- c("random", "mfi", "b_optimal", "fixed")
  next_item_rules <- list(
    random = list(par_names = NULL),
    mfi = list(par_names = NULL),
    # b_optimal = list(par_names = NULL),
    fixed = list(par_names = c("item_id")),
    mepv = list(par_names = c("var_calc_method"))
    )
  ability_est_rules <- list(
    eap = list(par_names = c("prior_dist", "prior_par", "min_theta",
                             "max_theta", "no_of_quadrature")),
    ml = list(par_names = c("min_theta", "max_theta", "criterion")),
    # "eap_ml" = list(par_names = c("prior_dist", "prior_par", "min_theta",
    #                               "max_theta", "no_of_quadrature")),
    # "sum_score" = list(par_names = NULL),
    owen = list(par_names = c("prior_mean", "prior_var")),
    sum_score = list(par_names = NULL)
    )
  exposure_control_rules <- list(
    randomesque  = list(par_names = c("num_items")),
    `sympson-hetter` = list(par_names = NULL)
    )
  # content_bal_rule <- c("max_discrepancy")
  content_bal_rules <- c()
  termination_rules <- list(
    min_item = list(par_names = "min_item"),
    max_item = list(par_names = "max_item"),
    min_se = list(par_names = "min_se"),
    sprt = list(par_names = c("theta_0", "theta_1", "alpha", "beta"))
    )
  ability_types <- c(
    # "multi_theta", "cdm", "raw_score",
    "theta"
    )


  is_single_string <- function(x) length(x) && is.character(x)

                        ####################@###
                        ### Item Pool Checks ###
                        ####################@###
  check_itempool <- function() {
    # This function check item pool for possible errors
    # ip should be an Itempool object.
    if (!is.null(ip)) {
      # Check whether 'ip' is a valid 'itempool' object
      if (!inherits(ip, "Itempool") || !is(ip, "Itempool"))
        stop("'ip' should be an 'Itempool' object. ")
      # Check validity of ip
      validObject(ip)
      # Rule: item pool should have unique id's
      if (any(duplicated(ip$id)))
        stop("Items in the item pool should have unique IDs.")
      # Rule: item pool size should not be smaller than the maximum test length.
      if (get_itempool_size(ip)["items"] < max_test_length)
        stop("Item pool size should not be smaller than the maximum test ",
             "length.")
    }

    # Check true_ip.
    # It should have the same length of ip:
    if (!is.null(true_ip)) {
      if (!inherits(true_ip, "Itempool"))
        stop("'true_ip' should be an 'Itempool' object. ")
      if (length(ip) != length(true_ip))
        stop("'true_ip' should have the same length as 'ip'.")
      if (!all(true_ip$id %in% ip$id))
        stop("Id's of 'ip' and 'true_ip' should be the same. ")
    }
    return(TRUE)
  }

                        ##########################@###
                        ### First Item Rule Checks ###
                        ##########################@###
  check_first_item_rule <- function() {
    # This function checks whether first item rule and parameters are valid.
    # First item rule can be NULL, if so, the next_item_parameter will be used.
    if (!is.null(first_item_rule)) {
      # check whether the rule is a single string
      if (!is_single_string(first_item_rule))
        stop(paste0("'first_item_rule' should be a string like:\n",
                    paste0("'", names(first_item_rules), "'", collapse = ", ")))
      if (!first_item_rule %in% names(first_item_rules))
        stop(paste0("'first_item_rule' should be a string with a value either: ",
                    paste0("'", names(first_item_rules), "'", collapse = ", "),
                    "."))
      # Check first_item_pars
      par_names <- first_item_rules[[first_item_rule]]$par_names
      if (!is.null(par_names) &&
          !(names(first_item_par) %in% par_names &&
          par_names %in% names(first_item_par))
          )
        stop(paste0("Invalid 'first_item_par'. The correct specification ",
                    "should look like this:\nfirst_item_par = list(",
                    paste0(par_names, " = ....", collapse = ", "), ")"
                    ))
    }
    return(TRUE)
  }

                        #########################@###
                        ### Next Item Rule Checks ###
                        #########################@###
  get_next_item_par_structure <- function() {
    # At the end, the next item parameter should be formed in two ways:
    # Structure "0":
    #   No next_item_par necessary, it can be NULL
    #
    # Structure "1":
    #   next_item_par = list(list(item_id = 'i3'), list(item_id = 'i2'),
    #                        list(item_id = 'i4'), list(item_id = 'i5'),
    #                        list(item_id = 'i1'))
    #
    # Structure "2":
    #   next_item_par = list(item_id = c("i3", "i2", "i4", "i5", "i1"))
    #
    # Structure "3":
    #   next_item_par = list(var_calc_method = "eap")
    #
    if (is.null(next_item_rules[[next_item_rule]]$par_names)) {
      next_item_par_structure = 0
    } else if (
      is.list(next_item_par) &&
      length(next_item_par) == max_test_length &&
      all(sapply(next_item_par, names) ==
          next_item_rules[[next_item_rule]]$par_names)
      ) {
      next_item_par_structure = 1
    } else if (
        is.list(next_item_par) &&
        length(next_item_par) == 1 &&
        names(next_item_par) == next_item_rules[[next_item_rule]]$par_names &&
        # Either the length of next_item_par is equal to the max_test_length or
        (length(next_item_par[[next_item_rules[[next_item_rule]]$par_names]]) ==
        max_test_length ||
        # or there are testlets in the item pool and the "fixed" specifies the
        # test lengths of the testlet items and standalone items.
        (
         next_item_rule == "fixed" &&
         "item_id" %in% names(next_item_par) &&
         !is.null(ip) &&
         any(sapply(ip@item_list, is, "Testlet")) &&
         ip[unique(next_item_par$item_id)]$n$items >= max_test_length &&
         (
           # Sometimes there is only one testlet so, the selection raises error
           length(unique(next_item_par$item_id)[-length(
           unique(next_item_par$item_id))]) == 0 ||
           ip[unique(next_item_par$item_id)[-length(
             unique(next_item_par$item_id))]]$n$items <= max_test_length
           )
         )
        )
      ) {
      next_item_par_structure = 2
    } else if (
      is.list(next_item_par) &&
      length(next_item_par) == 1 &&
      names(next_item_par) == next_item_rules[[next_item_rule]]$par_names &&
      length(next_item_par[[next_item_rules[[next_item_rule]]$par_names]]) == 1
      ) {
      next_item_par_structure = 3
    } else
      stop("'next_item_par' does not have an acceptable format. Please see ",
           "?create_cat_design for examples.")
    return(next_item_par_structure)
  }
  check_next_item_rule <- function() {
    # This function checks whether next item rule and parameters are valid.
    # Next item rule cannot be empty and it should be a vector of
    # valid rules.
    if (!all(next_item_rule %in% names(next_item_rules)))
      stop(paste0("next_item_rule should be a vector with elements either: ",
                  paste0("'", names(next_item_rules), "'", collapse = ", "),
                  "."))

    next_item_par_structure <- get_next_item_par_structure()

    # Make sure that the length of next_item_par and next_item_rule are the
    # same.
    if ((!is.null(next_item_par)) &&
        # check if it is a list of lists.
        all(sapply(next_item_par, is, "list")) &&
        (length(next_item_rule) != 1) &&
        (length(next_item_rule) != length(next_item_par)))
      stop("The length of next_item_rule should be equal to the length of
           next_item_par.")
    next_item_missing_par_error_text <- paste0(
      "'next_item_par' should be a list object. \n", ifelse(
        is.null(next_item_rules[[next_item_rule]]$par_names), "",
        paste0("When you specify 'next_item_rule' = '", next_item_rule,
               "', you need to add the following argument: \n ",
               "next_item_par = list(",
               paste0(next_item_rules[[next_item_rule]]$par_names,
                      collapse = " = ...,"), " = ...)\n")
          ))
    # If there should be a next item parameters, raise an error.
    if (!is.null(next_item_rules[[next_item_rule]]$par_names) &&
        is.null(next_item_par))
      stop(next_item_missing_par_error_text)
    if (!is.null(next_item_par)  && !is(next_item_par, "list"))
      stop(next_item_missing_par_error_text)

    # If next_item_rule is a vector with length larger than one, it's size should
    # be equal to the max_test_length (maximum test length)
    if ((length(next_item_rule) > 1) &&
        (length(next_item_rule) != max_test_length))
      stop("The length of next_item_rule should be equal to the maximum item
           length.")
    if (!is.null(next_item_par) &&
        all(sapply(next_item_par, is, "list")) &&
        (length(next_item_par) > 1) &&
        (length(next_item_par) != max_test_length))
      stop("The length of next_item_par should be equal to the maximum item
           length.")

    # Check the validity of Parameter values
    #
    # The item selection algorithm expects this.
    if (next_item_rule == "fixed") {
      # If the next item selection method is "fixed", there should be a
      # valid item pool (ip) argument.
      if (is.null(ip) || !is(ip, "Itempool"))
        stop("There should be a valid item pool (argument) for next item
             rule 'fixed' to work. ")
      # next_item_par should be a list of item_id's
      if (!is(next_item_par, "list"))
        stop("next_item_par should be a list object.")
      # Check which structure does the paramters obey:

      # All of the next_item_par elements should have a parameter named
      # "item_id"
      # The following check might be redundant but the error is informative so
      # keep it.
      if (any(sapply(next_item_par, "names") != "item_id") && !( # Structure 1
        # Or it can be something like: # Structure 2
        # next_item_par = list(item_id = c("i3", "i2", "i4", "i5", "i1"))
        is.list(next_item_par) &&
        length(next_item_par) == 1 &&
        names(next_item_par) == next_item_rules[[next_item_rule]]$par_names &&
        # Either the length of next_item_par is equal to the max_test_length or
        (length(next_item_par[[next_item_rules[[next_item_rule]]$par_names]]) ==
        max_test_length ||
        # or there are testlets in the item pool and the "fixed" specifies the
        # test lengths of the testlet items and standalone items.
        (any(sapply(ip@item_list, is, "Testlet")) &&
         ip[next_item_par$item_id]$n$items >= max_test_length &&
         (length(next_item_par$item_id[-length(next_item_par$item_id)])  == 0 ||
          ip[next_item_par$item_id[-length(next_item_par$item_id)]]$n$items <=
          max_test_length
         )
         )
        ))
        )
        stop(paste0(
          "If 'next_item_rule' = 'fixed', the 'next_item_par' should be ",
          "like:\n   ",
          "next_item_par = list('item_id' = c(<THE ORDERED ITEM IDs>))",
          "\nAlso, the length of the 'item_id' vector should be ", max_test_length,
          "."))
      # All elements should be unique
      if (
        (next_item_par_structure == 1 &&
         any(duplicated(sapply(next_item_par, "[[", "item_id")))) ||
        (next_item_par_structure == 2 &&
         any(duplicated(
           next_item_par[[next_item_rules[[next_item_rule]]$par_names]])))
      )
        stop("All item_id's should be unique. Please check 'next_item_par'.")
      # All of the item_id field's should be within the id's of item pool (ip)
      if ((next_item_par_structure == 1 &&
           !all(sapply(next_item_par, "[[", "item_id") %in% ip$id)) ||
          (next_item_par_structure == 2 &&
           !all(next_item_par[[
             next_item_rules[[next_item_rule]]$par_names]] %in% ip$id))
          )
        stop("All of the id's in 'item_id' field of next_item_par should be
             also in the item pool (ip) id's. Some of the 'item_id's are not
             valid.")
    }
    return(TRUE)
  }

                        #############################@###
                        ### Ability Estimation Checks ###
                        #############################@###
  check_ability_est_rule <- function(ae_rule, ae_par, final_ae = FALSE) {
    # This function checks whether ability estimation rule and parameters are
    # valid.
    # @param final_ae If TRUE, the checks will be performed for
    #   final_ability_est_rule and final_ability_est_par.
    # Ability estimation rule cannot be empty and it should be a vector of
    # valid rules
    if (!all(ae_rule %in% names(ability_est_rules)))
      stop(paste0(ifelse(final_ae, "'final_", "'"),
                  "ability_est_rule' should be a vector with elements either: ",
                  paste0("'", names(ability_est_rules), "'", collapse = ", "),
                  "."))

    # If ability_est_rule is a vector with length larger than one, it's size
    # should be equal to the max_test_length (maximum test length)
    if (!final_ae &&
        (length(ability_est_rule) > 1) &&
        (length(ability_est_rule) != max_test_length))
      stop("The length of ability_est_rule should be equal to the maximum item
           length.")

    # Make sure that the length of ability_est_rule and ability_est_par are the
    # same.
    if (!is.null(ae_par)) {
      # Make sure the ae_par is a list object
      if (!is(ae_par, "list"))
        stop(ifelse(final_ae, "'final_", "'"),
             "ability_est_par' should be a list object.")
      # check if it is a list of lists.
      par_names <- ability_est_rules[[ae_rule]]$par_names
      if (!all(sapply(ae_par, is.list))) { # if it is not a list of lists:
        if (!all(names(ae_par) %in% par_names) ||
            !all(par_names %in% names(ae_par)))
          stop(paste0("Invalid ", ifelse(final_ae, "'final_", "'"),
                      "ability_est_par'. Please provide parameter names ",
                      "like:\n'", ifelse(final_ae, "final_", ""),
                      "ability_est_par = list(",
                      paste0(par_names, collapse = " = , "), " = )"))
      } else if (!final_ae &&  # If parameters of a list of lists
                 length(ae_par) == max_test_length &&
                 # All elements of ae_par is list
                 all(sapply(ae_par, is.list))
                 ) {
        for (i in ae_par) {
          if (!all(names(i) %in% par_names) ||
              !all(par_names %in% names(i)))
            stop(paste0("Invalid ", ifelse(final_ae, "'final_", "'"),
                        "ability_est_par'. Please provide parameter names ",
                        "like:\n'", ifelse(final_ae, "final_", ""),
                        "ability_est_par = list(",
                        paste0(par_names, collapse = " = , "), " = )"))
        }
      } else {
        stop("The length of ability_est_rule should be equal to the length of
             ability_est_par.")
      }
      # if (all(sapply(ability_est_par, is, "list")) &&
      #     (length(ability_est_par) > 1) &&
      #     (length(ability_est_par) != max_test_length))
      #   stop("The length of ability_est_par should be equal to the maximum item
      #        length.")
    }
    return(TRUE)
  }

                        ###################################@###
                        ### Final Ability Estimation Checks ###
                        ###################################@###
  # check_final_ability_est_rule <- function() {
  #   # This function checks whether the final ability estimation rule and
  #   # parameters are valid.
  #   if (is.null(final_ability_est_rule)) return(TRUE)
  #   # final_ability_est_rule should be a string.
  #   if (!is_single_string(final_ability_est_rule) ||
  #       (!final_ability_est_rule %in% names(final_ability_est_rules))
  #       )
  #     stop(paste0("'final_ability_est_rule' should be a string with elements ",
  #                 "either: \n", paste0("'", names(ability_est_rules), "'",
  #                                    collapse = ", ")))
  #
  #   # if ((!is.null(ability_est_par)) &&
  #   #     # check if it is a list of lists.
  #   #     all(sapply(ability_est_par, is, "list")) &&
  #   #     (length(ability_est_rule) != 1) &&
  #   #     (length(ability_est_rule) != length(ability_est_par)))
  #   #   stop("The length of ability_est_rule should be equal to the length of
  #   #        ability_est_par.")
  #   # if (!is.null(ability_est_par)  && !is(ability_est_par, "list"))
  #   #   stop("ability_est_par should be a list object.")
  #   # if (!is.null(ability_est_par) &&
  #   #     all(sapply(ability_est_par, is, "list")) &&
  #   #     (length(ability_est_par) > 1) &&
  #   #     (length(ability_est_par) != max_test_length))
  #   #   stop("The length of ability_est_par should be equal to the maximum item
  #   #        length.")
  #   return(TRUE)
  # }
                        ###########################@###
                        ### Exposure Control Checks ###
                        ###########################@###
  check_exposure_control_rule <- function() {
    # This function checks whether exposure control rule and parameters are
    # valid.
    # Exposure control rule is either NULL or a vector of valid rules

    if (!is.null(exposure_control_rule) &&
        !all(exposure_control_rule %in% names(exposure_control_rules)))
      stop(paste0("exposure_control_rule should be a vector with elements either: ",
                  paste0("'", names(exposure_control_rules), "'", collapse = ", "),
                  "."))
    if (!is.null(exposure_control_rule))
    {
      # Make sure that the length of exposure_control_rule and
      # exposure_control_par are the same.
      if ((!is.null(exposure_control_par)) &&
          # check if it is a list of lists.
          all(sapply(exposure_control_par, is, "list")) &&
          (length(exposure_control_rule) != 1) &&
          (length(exposure_control_rule) != length(exposure_control_par)))
        stop("The length of exposure_control_rule should be equal to the ",
             "length of exposure_control_par.")
      if (!is.null(exposure_control_par) &&
          !is(exposure_control_par, "list"))
        stop("exposure_control_par should be a list object.")
      # If exposure_control_rule is a vector with length larger than one, it's
      # size should be equal to the max_test_length (maximum test length)
      if ((length(exposure_control_rule) > 1) &&
          (length(exposure_control_rule) != max_test_length))
        stop("The length of exposure_control_rule should be equal to the ",
             "maximum item length.")
      if (!is.null(exposure_control_par) &&
          all(sapply(exposure_control_par, is, "list")) &&
          (length(exposure_control_par) > 1) &&
          (length(exposure_control_par) != max_test_length))
        stop("The length of exposure_control_par should be equal to the ",
             "maximum item length.")
      ### "sympson-hetter" checks
      # If the exposure_control_rule is "sympson-hetter", then each Item/testlet
      # of the "ip" should have a parameter for "sympson-hetter" method.
      if (any(exposure_control_rule %in% "sympson-hetter")) {
        if (is.null(ip) || !is(ip, "Itempool"))
          stop("For 'sympson-hetter' exposure control rule, there should be a
               valid item pool (ip) in the arguments. ")
        for (item in ip@item_list)
          if (!"sympson_hetter_k" %in% names(item@misc) ||
              item@misc[["sympson_hetter_k"]] < 0 ||
              item@misc[["sympson_hetter_k"]] > 1)
            stop(paste0(
            "For 'sympson-hetter' exposure control rule, there should be valid
            'sympson_hetter_k' values for each item. Make sure to check for
            each item 'item$misc' and ensure that it has an element named
            'sympson_hetter_k' which is between 0 and 1. For an item, you can
            use 'add_misc(item, list(sympson_hetter_k = .75))', to add
            that value."))
        if (sum(sapply(ip, function(k) k@misc$sympson_hetter_k) == 1) < max_test_length)
          warning(paste0("When using 'sympson-hetter' exposure control rule, ",
                         "please ensure that there are at least ",
                         max_test_length, " items with 'sympson_hetter_k' ",
                         "values 1. Otherwise, examinees might not get a ",
                         "complete test and an error might be raised by ",
                         "the simulation function."), call. = FALSE)
      }
    }
    return(TRUE)
  }

                        #################################@###
                        ### Content Balancing Rule Checks ###
                        #################################@###
  check_content_balancing_rule <- function() {
    # This function checks whether content balancing rule and parameters are
    # valid.
    # Content balancing rule is either NULL or a vector of valid rules

    if (!is.null(content_bal_rule) &&
        !all(content_bal_rule %in% content_bal_rule))
      stop(paste0("content_bal_rule should be a vector with elements either: ",
                  paste0("'", content_bal_rule, "'", collapse = ", "),
                  "."))
    # Make sure that the length of content_bal_rule and exposure_control_par
    # are the same.
    if ((!is.null(content_bal_rule)) &&
        (!is.null(content_bal_par)) &&
        # check if it is a list of lists.
        all(sapply(content_bal_par, is, "list")) &&
        (length(content_bal_rule) != 1) &&
        (length(content_bal_rule) != length(content_bal_par)))
      stop("The length of content_bal_rule should be equal to the length of
           content_bal_par.")
    if (!is.null(content_bal_par)  && !is(content_bal_par, "list"))
      stop("content_bal_par should be a list object.")
    # If content_bal_rule is a vector with length larger than one, it's size
    # should be equal to the max_test_length (maximum test length)
    if ((!is.null(content_bal_rule)) && (length(content_bal_rule) > 1) &&
        (length(content_bal_rule) != max_test_length))
      stop("The length of content_bal_rule should be equal to the maximum item
           length.")
    if (!is.null(content_bal_par) &&
        all(sapply(content_bal_par, is, "list")) &&
        (length(content_bal_par) > 1) &&
        (length(content_bal_par) != max_test_length))
      stop("The length of content_bal_par should be equal to the maximum item
           length.")
    return(TRUE)
  }

                        ###########################@###
                        ### Termination Rule Checks ###
                        ###########################@###
  # Function to determine the input structure of the termination_par.
  #
  # There are two structures each element of 'termination_par' can have:
  #   Structure (1):
  #     termination_par = list(min_item = 10)
  #   Structure (2):
  #     termination_par = list(min_item = list(min_item = 10))
  #     # For "sprt" only one option is available
  #     termination_par = list(min_item = 10,
  #                            sprt = list(theta_0 = -1, theta_1 = 1,
  #                                        alpha = 0.05, beta = 0.05))
  #
  # This function determines which structure it has and returns either 1 or 2
  # or NULL if neither structures fit.
  #
  # @param tr Individual termination rule such as "min_se", "min_item"
  #   or "sprt"
  get_termination_par_structure <- function(tr) {
    tp_structure <- NULL
    # Get the parameter count
    pars <- termination_rules[[tr]]$par_names
    if (length(pars) == 1) { # either min_se, min_item or max_item
      if (is.list(termination_par[[tr]]) &&
          names(termination_par[[tr]]) == pars) {
        tp_structure = 2
      } else if (is.numeric(termination_par[[tr]])) {
        tp_structure = 1
      }
    } else if (length(pars) > 1 &&
               is.list(termination_par[[tr]]) &&
               all(pars %in% names(termination_par[[tr]]))
               ) {
      tp_structure = 2
    }
    return(tp_structure)
  }
  check_termination_rule <- function() {
    # The available termination rules:
    # Make sure that the length of termination_par and termination_rule are the
    # same.
    if (!is(termination_rule, "character")) # This includes NULL too.
      stop("termination_rule should be a vector of 'character'.")
    if (!is(termination_par, "list"))
      stop("termination_par should be a list object.")
    if (length(termination_rule) != length(termination_par))
      stop("The length of termination_rule should be equal to the length of
           termination_par.")
    # The elements of the termination_rule should be valid
    if (!all(termination_rule %in% names(termination_rules)))
      stop(paste0("termination_rule should be a vector with elements either: ",
                  paste0(names(termination_rules), collapse = ", "), "."))
    # The elements of the termination_par should be valid
    if (!all(names(termination_par) %in% names(termination_rules)))
      stop(paste0("The names of the elements of the termination_par should ",
                  "match one of the following: ",
                  paste0(termination_rules, collapse = ", "), "."))

    if (!all(names(termination_par) %in% termination_rule))
        stop("The names of the elements of termination_par should match the
              the elements of 'termination_rule'.")
    if (!"max_item" %in% termination_rule && is.null(ip))
      stop("If max_item is not found in termination_rule, then the
      max_item will be set to the size of the item pool (ip). Make sure
      'ip' is present if 'max_item' is missing. ")
    # Further check the termination_par:
    # If the length of the termination_rules$termination_rule$par_names is 1,
    # then there are two structures:
    for (tr in termination_rule) { # tr: individual termination rule
      tr_structure <- get_termination_par_structure(tr)
      if (is.null(tr_structure))
        stop(paste0("'termination_par' element ", tr, " does not have ",
                    "acceptable value. Please see ?create_cat_design."))
    }
    return(TRUE)
  }


                        ########################@###
                        ### Ability Types Checks ###
                        ########################@###
  if (is.null(ability_type) || (!ability_type %in% ability_types))
    stop(paste0("'ability_type' should be one of the following:\n",
                paste0("'", ability_types, "'", collapse = ", ")))



  ##########################################################################@###
  ######################### Start Function #################################@###
  ##########################################################################@###
  # Convert ip to Itempool object
  if (!is.null(ip) && !is(ip, "Itempool"))
    stop(paste0("\nInvalid item pool. Please provide an 'Itempool' object ",
                "for 'ip' argument. See 'itempool()' function. \n"))
    # tryCatch(
    #   ip <- itempool(ip),
    #   error = function(cond) {
    #     message(paste0("\nInvalid item pool. ip cannot be converted to ",
    #                    "an 'Itempool' object. \n"))
    #     stop(cond, call. = FALSE)
    #   })

  # Check whether termination_rule and termination_par are valid. These
  # arguments will be used for finding the max_test_length, maximum possible
  # length of the test. There are instances where this parameter cannot be set,
  # for example, when there is an infinite item pool (i.e. ip = NULL) and the
  # stopping rule is minimum standard error ("min_se") there cannot be a
  # max_test_length. For those instances, max_test_length is set to 10,000 in
  # order the test to converge.
  # max_test_length is based on individual items within testlets and standalone
  # items.
  max_test_length <- ifelse(is.null(ip), 10000, ip$n$items)
  check_termination_rule()
  # Find the maximum number of items, i.e. maximum test length:
      # This checks whether termination_par is a list of lists
  if (all(sapply(termination_par, is, "list")) &&
      # if max_item is not within the rule then set the max_item as the size of
      # the item pool
      all(sapply(termination_rule, FUN = function(x) "max_item" %in% x))
      ) {
      max_test_length <- max(sapply(termination_par,
                               FUN = function(x) x[["max_item"]]))
  } else if ("max_item" %in% names(termination_par))
    max_test_length <- termination_par[["max_item"]]

  ##########################################################################@###
  ################## Set Defaults for Parameters if Missing ################@###
  ##########################################################################@###
  if (is.null(ability_est_rule)) ability_est_rule <- "eap"

  ### ability_est_par ###
  if (is.null(ability_est_par)) {
    ability_est_par <- switch (ability_est_rule,
      "eap" = list(prior_dist = "norm", prior_par = c(0, 1),
                   min_theta = -4, max_theta = 4, no_of_quadrature = 50),
      "owen" = list(prior_mean = 0, prior_var = 1),
      "ml" = list(min_theta = -4, max_theta = 4, criterion = 0.001)
    )
  }
  ### final_ability_est_par ###
  if (!is.null(final_ability_est_rule) && is.null(final_ability_est_par)) {
    final_ability_est_par <- switch (final_ability_est_rule,
      "eap" = list(prior_dist = "norm", prior_par = c(0, 1),
                   min_theta = -4, max_theta = 4, no_of_quadrature = 50),
      "owen" = list(prior_mean = 0, prior_var = 1),
      "ml" = list(min_theta = -4, max_theta = 4, criterion = 0.001)
    )
  }

  ### next_item_par ###
  # if the next_item_rule selected is "mepv" but the next_item_par has not
  # been set, set a default value to it.
  if (length(next_item_rule) == 1 && next_item_rule == "mepv" &&
      is.null(next_item_par))
    next_item_par <- list(var_calc_method = rep("eap", max_test_length))

  # if next_item_rule is "fixed" and next_item_par = NULL and ip is not NULL
  # then use first n items in the item pool such that number of items in the
  # first n items in the item pool is just larger than the test length.
  if (next_item_rule == "fixed" && is.null(next_item_par) &&
      is(ip, "Itempool")) {
    next_item_par <- list(item_id = ip$id[1:which(sapply(1:length(ip), function(i)
      ip[1:i]$n$items) >= max_test_length)[1]])
  }

  ##########################################################################@###
  ##########################################################################@###
  cd <- list()
  # Check the title
  if (is.null(title) || (is.character(title) && length(title) == 1)) {
    cd$title <- title
  } else
    stop("Invalid 'title'. Please provide a valid string value.")

  cd$ability_type <- ability_type
  cd$max_test_length <- max_test_length

  cd$ip <- ip
  if (!is.null(true_ip)) cd$true_ip = true_ip
  # design_list$step will hold information regarding each step of the adaptive
  # test.
  cd$step <- replicate(max_test_length, list())
  names(cd$step) <- paste0(1:max_test_length)

  # Check whether next item rule makes sense. These functions needs
  # max_test_length.
  check_first_item_rule()
  check_next_item_rule()
  check_ability_est_rule(ae_rule = ability_est_rule, ae_par = ability_est_par,
                         final_ae = FALSE)
  check_ability_est_rule(ae_rule = ability_est_rule, ae_par = ability_est_par,
                         final_ae = TRUE)
  check_exposure_control_rule()
  check_content_balancing_rule()
  check_itempool()

  # -------------------------------------------------------------------------- #
  # Set the First Item criteria
  cd$first_item_rule <- first_item_rule
  cd$first_item_par <- first_item_par

  # -------------------------------------------------------------------------- #
  # Set the Next Item criteria

  # If next_item_rule is "fixed" and there are testlets in the "item_ids", then
  # repeat the testlets
  if (next_item_rule == "fixed" &&
      length(next_item_par$item_id) != max_test_length)
    next_item_par$item_id <- rep(next_item_par$item_id, times = sapply(
      ip[next_item_par$item_id]$item_list, length))

  for (i in 1:length(cd$step)) {
    ### Set next_item_rule for each step ###
    if (length(next_item_rule) > 1) {
      cd$step[[i]]$next_item_rule = next_item_rule[i]
    } else if (length(next_item_rule) == 1) {
      cd$step[[i]]$next_item_rule = next_item_rule
    } else
      stop("The length of next_item_rule should be larger than 0.")
    ### Set next_item_par for each step  ###
    # Determine the structure of the next_item_par
    next_item_par_structure <- get_next_item_par_structure()
    if (!is.null(next_item_par))
      # The following test assumes that if all elements of the next_item_par is
      # a list object, then it's length should be equal to the max_test_length.
      # In other words, there is a new parameter set for each item at each CAT
      # step. In a rare occasion, there is a possibility that this assumption
      # may not hold. All elements are a list object but the given list is
      # for all elements, and it should be replicated for max_test_length.
      if (next_item_par_structure == 1) {
        # if (length(next_item_par) != max_test_length)
        #   stop("The length of the next_item_par should equal to max
        #   test length.!")
        cd$step[[i]]$next_item_par = next_item_par[[i]]
        # If the next_item_par is something like:
        # next_item_par = list(item_id = c("i3", "i2", "i4", "i5", "i1"))
      } else if (next_item_par_structure == 2) {
        cd$step[[i]]$next_item_par[[
          next_item_rules[[next_item_rule]]$par_names]] =
          next_item_par[[next_item_rules[[next_item_rule]]$par_names]][i]
      } else if (next_item_par_structure == 3) {
        cd$step[[i]]$next_item_par[[
          next_item_rules[[next_item_rule]]$par_names]] =
          next_item_par[[next_item_rules[[next_item_rule]]$par_names]]
      } else cd$step[[i]]$next_item_par = next_item_par
  }
  # -------------------------------------------------------------------------- #
  # Set the Ability Estimation Rules
  for (i in 1:length(cd$step)) {
    if (length(ability_est_rule) > 1) {
      cd$step[[i]]$ability_est_rule = ability_est_rule[i]
    } else if (length(ability_est_rule) == 1) {
      cd$step[[i]]$ability_est_rule = ability_est_rule
    } else
      stop("The length of ability_est_rule should be larger than 0.")
    if (!is.null(ability_est_par))
      # The following test assumes that if all elements of the ability_est_par
      # is a list object, then it's length should be equal to the
      # max_test_length. In other words, there is a new parameter set for each
      # item at each CAT step.
      if (all(sapply(ability_est_par, is, "list"))) {
        cd$step[[i]]$ability_est_par = ability_est_par[[i]]
      } else cd$step[[i]]$ability_est_par = ability_est_par
  }

  # -------------------------------------------------------------------------- #
  # Set the final_ability_est_rule and final_ability_est_par
  if (is.null(final_ability_est_rule)) {
    cd$final_ability_est_rule <- NULL
    cd$final_ability_est_par <- NULL
  } else {
    cd$final_ability_est_rule <- final_ability_est_rule
    cd$final_ability_est_par <- final_ability_est_par
  }

  # -------------------------------------------------------------------------- #
  # Set the Termination Criteria
  cd$termination_rule <- termination_rule
  # All termination_par element should follow structure 2 above (see function
  # description of 'get_termination_par_structure()')
  for (tr in termination_rule) { # tr: individual termination rule
    tr_structure <- get_termination_par_structure(tr)
    if (tr_structure == 1) {
      termination_par[[tr]] <- list(termination_par[[tr]])
      names(termination_par[[tr]]) <- tr
    }
  }
  cd$termination_par <- termination_par

  # -------------------------------------------------------------------------- #
  # Set Exposure Control Rules
  for (i in 1:length(cd$step)) {
    if (!is.null(exposure_control_rule)) {
      if (length(exposure_control_rule) > 1) {
        cd$step[[i]]$exposure_control_rule = exposure_control_rule[i]
      } else if (length(exposure_control_rule) == 1) {
        cd$step[[i]]$exposure_control_rule = exposure_control_rule
      } else
        stop("The length of exposure_control_rule should be larger than 0.")
    }
    if (!is.null(exposure_control_par))
      # The following test assumes that if all elements of the
      # exposure_control_par are a list object, then it's length should be equal
      # to the max_test_length. In other words, there is a new parameter set for
      # each item at each CAT step.
      if (all(sapply(exposure_control_par, is, "list"))) {
        cd$step[[i]]$exposure_control_par = exposure_control_par[[i]]
      } else cd$step[[i]]$exposure_control_par = exposure_control_par

    # Randomesque Checks
    # If it is randomesque, there should be a parameter for it.
    if (!is.null(exposure_control_rule) &&
        cd$step[[i]]$exposure_control_rule == "randomesque" &&
        !("num_items" %in% names(cd$step[[i]]$exposure_control_par)))
      stop("For randomesque exposure control, 'num_items' parameter should be
           provided.")

  }

  # -------------------------------------------------------------------------- #
  # Set Content Balancing Rules
  for (i in 1:length(cd$step)) {
    if (!is.null(content_bal_rule)) {
      if (length(content_bal_rule) > 1) {
        cd$step[[i]]$content_bal_rule = content_bal_rule[i]
      } else if (length(content_bal_rule) == 1) {
        cd$step[[i]]$content_bal_rule = content_bal_rule
      } else
        stop("The length of content_bal_rule should be larger than 0.")
    }
    if (!is.null(content_bal_par))
      # The following test assumes that if all elements of the content_bal_par
      # is a list object, then it's length should be equal to the
      # max_test_length. In other words, there is a new parameter set for each
      # item at each CAT step.
      if (all(sapply(content_bal_par, is, "list"))) {
        cd$step[[i]]$content_bal_par = content_bal_par[[i]]
      } else cd$step[[i]]$content_bal_par = content_bal_par
  }

  # -------------------------------------------------------------------------- #
  class(cd) <- "cat_design"
  return(cd)
}


###############################################################################@
############################# cat_sim ##########################################
###############################################################################@
#' Computerized Adaptive Test (CAT) Simulation
#'
#' @description
#' \code{cat_sim} function simulates computerized adaptive test (CAT) for
#' one or more simulees. For long simulations, \code{\link{cat_sim_fast}}
#' function can be used.
#'
#' @param true_ability True ability vector to generate item responses.
#' @param cd A \code{cat_design} object that is created by function
#'          \code{create_cat_design}.
#' @param verbose This is an integer that will print the stage of the test.
#'   For example, if the value verbose = 10, a message will be printed at
#'   each tenth iteration of the cat_simulation. Default value is \code{-1},
#'   where no message will be printed. If the value is \code{0}, only the
#'   start time and end time of the simulation will be printed.
#'
#' @return If the length of \code{true_ability} vector is one a
#'   \code{"cat_output"} class output will be returned.
#'   This is a list containing following elements:
#'   \describe{
#'     \item{true_ability}{True ability (theta) value to generate item
#'       responses.}
#'     \item{est_history}{A list where each element represent a step of the
#'       CAT test. It has following elements:
#'       \describe{
#'         \item{est_before}{The estimated ability before the administration
#'           of the item. }
#'         \item{se_before}{The standard error of the estimated ability before
#'           the administration of the item. }
#'         \item{testlet}{\code{TRUE} if the item belongs to a testlet.}
#'         \item{item}{\code{\link{Item-class}} object that is administered at
#'           this step.}
#'         \item{resp}{The simulated response of the simulee for the item
#'           administered at this step using simulee's \code{true_ability}
#'           value.}
#'         \item{est_after}{The estimated ability after the administration
#'           of the item.}
#'         \item{se_after}{The standard error of the estimated ability after
#'           the administration of the item. }
#'        }
#'     }
#'   }
#'
#'   If the length of the \code{true_ability} is more than 1, a list of
#'   \code{cat_output} objects will be returned for each value of
#'   \code{true_ability}.
#'
#' @seealso \code{\link{create_cat_design}}
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @examples
#' ip <- generate_ip(n = 50)
#' # Check the default:
#' cd <- create_cat_design(ip = ip)
#' cat_sim(true_ability = rnorm(1), cd = cd)
cat_sim <- function(true_ability, cd, verbose = -1)
{
  # Make sure that cat_design is a cat_design object.
  if (!inherits(cd, "cat_design"))
    stop("'cat_design' should be a 'cat_design' object. Please run
         'create_cat_design' function.")
  # Convert true_ability to list in case it is numeric. This is to take care
  # of future expansion for possible representation of ability of an examinee
  # using multiple values (like MIRT or CDM)
  if (!is(true_ability, "list"))
    true_ability <- lapply(true_ability, function(x) x)
  started <- Sys.time()
  result <- cat_sim_cpp(true_ability, cd, verbose = as.integer(max(0, verbose)))

  # Print start and end time of the simulation
  if (verbose >= 0) {
    now <- Sys.time()
    cat(paste0("\nSimulation started at ", format(started, format = "%X"), " (",
               format(started, format = "%x"), ") and ended at ",
               ifelse(as.Date(now) == as.Date(started),
                      format(now, format = "%X"), now), ".\n"))
  }
  return(result)
}



###############################################################################@
############################# cat_sim_fast #####################################
###############################################################################@
#' Computerized Adaptive Test (CAT) Simulation (Parallel Computing)
#'
#' @description
#' \code{cat_sim_fast} function simulates computerized adaptive test (CAT) for
#' one or many simulees. This function uses parallel computing, so, for large
#' number of simulees, it might be significantly faster than
#' \code{\link{cat_sim}} function.
#'
#' @param true_ability True ability vector to generate item responses.
#' @param cd A \code{cat_design} object that is created by function
#'          \code{create_cat_design}.
#' @param verbose This is an integer that will print the stage of the test.
#'   For example, if the value verbose = 10, a message will be printed at
#'   each tenth iteration of the cat_simulation. Default value is \code{-1},
#'   where no message will be printed. If the value is \code{0}, only the
#'   start time and end time of the simulation will be printed.
#' @param n_cores an integer specifying the number of cores to be used.
#'   The value should be 1 or larger. The default is \code{NULL} where
#'   the maximum number of cores of the processor will be used.
#'
#' @return If the length of \code{true_ability} vector is one a
#'   \code{"cat_output"} class output will be returned.
#'   This is a list containing following elements:
#'   \describe{
#'     \item{true_ability}{True ability (theta) value to generate item
#'       responses.}
#'     \item{est_history}{A list where each element represent a step of the
#'       CAT test. It has following elements:
#'       \describe{
#'         \item{est_before}{The estimated ability before the administration
#'           of the item. }
#'         \item{se_before}{The standard error of the estimated ability before
#'           the administration of the item. }
#'         \item{testlet}{\code{TRUE} if the item belongs to a testlet.}
#'         \item{item}{\code{\link{Item-class}} object that is administered at
#'           this step.}
#'         \item{resp}{The simulated response of the simulee for the item
#'           administered at this step using simulee's \code{true_ability}
#'           value.}
#'         \item{est_after}{The estimated ability after the administration
#'           of the item.}
#'         \item{se_after}{The standard error of the estimated ability after
#'           the administration of the item. }
#'        }
#'     }
#'   }
#'
#'   If the length of the \code{true_ability} is more than 1, a list of
#'   \code{cat_output} objects will be returned for each value of
#'   \code{true_ability}.
#'
#' @export
#'
#' @seealso \code{\link{create_cat_design}}
#'
#' @author Emre Gonulates
#'
#' @examples
#' cd <- create_cat_design(ip = generate_ip(n = 30),
#'                         termination_rule = c('max_item'),
#'                         termination_par = list(max_item = 7))
#' cat_sim_fast(true_ability = rnorm(1), cd = cd, n_cores = 1)
#'
#' cat_sim_fast(true_ability = rnorm(2), cd = cd, n_cores = 1)
#'
cat_sim_fast <- function(true_ability, cd, verbose = -1, n_cores = NULL)
{
  # Make sure that cat_design is a cat_design object.
  if (!inherits(cd, "cat_design"))
    stop("'cat_design' should be a 'cat_design' object. Please run
         'create_cat_design' function.")
  # Convert true_ability to list in case it is numeric. This is to take care
  # of future expansion for possible representation of ability of an examinee
  # using multiple values (like MIRT or CDM)
  if (!is(true_ability, "list"))
    true_ability <- lapply(true_ability, function(x) x)
  started <- Sys.time()
  n_theta <- length(true_ability)
  max_cores <- parallel::detectCores()
  n_cores <- ifelse(is.null(n_cores), max_cores, min(max(n_cores, 1),
                                                     max_cores))
  # Use parallel processing only when there are rather large number of
  # thetas.
  if (n_theta > max_cores) {
    cl <- parallel::makeCluster(n_cores)
    result <- parallel::parLapply(cl = cl, true_ability, cat_sim_single_cpp,
                                  cd = cd)
    parallel::stopCluster(cl)
  } else
    return(cat_sim_cpp(true_ability, cd, verbose = as.integer(max(0, verbose))))

  # Print start and end time of the simulation
  if (verbose >= 0) {
    now <- Sys.time()
    cat(paste0("\nSimulation started at ", format(started, format = "%X"), " (",
               format(started, format = "%x"), ") and ended at ",
               ifelse(as.Date(now) == as.Date(started),
                      format(now, format = "%X"), now), ".\n"))
  }

  if (length(true_ability) == 1) {
    return(result[[1]]);
  } else return(result)
}



#############################################################################@
########################### calculate_exposure_cont_pars_cpp #################
#############################################################################@
# Calculate exposure control parameters
#
# This function calculates the exposure control parameters for the
# Sympson-Hetter method.
#
#
#calculate_exposure_cont_pars_cpp <- function()
