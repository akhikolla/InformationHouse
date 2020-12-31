#' Breeding Bird Survey
#'
#' This dataset is a subset of the Breeding Bird Survey.
#'
#'
#' @name birds
#' @docType data
#' @format A data frame with 13608 observations and 384 variables.
#' \describe{
#'   \item{\code{loc.id}}{Location index}
#'   \item{\code{aou530}}{Bird specie 530 presence (1) or absence (0)}
#'   \item{\code{aou540}}{Bird specie 540 presence (1) or absence (0)}
#'   \item{\code{aou590}}{Bird specie 590 presence (1) or absence (0)}
#'   \item{\code{...}}{Bird species ... presence (1) or absence (0)}
#'   \item{\code{aou4461}}{Bird specie 4461 presence (1) or absence (0)}
#'   \item{\code{aou5860}}{Bird specie 5860 presence (1) or absence (0)}
#' }
#' @source Pardieck, K.L., D.J. Ziolkowski Jr., M. Lutmerding, K. Campbell and
#' M.-A.R. Hudson. 2017. North American Breeding Bird Survey Dataset 1966 -
#' 2016, version 2016.0. U.S. Geological Survey, Patuxent Wildlife Research
#' Center. \url{https://www.pwrc.usgs.gov/bbs/RawData/}doi:10.5066/F7W0944J.
#' @keywords datasets
#' @usage data(birds)
#' @examples
#' data(birds)
#'
NULL

#' Complaints received for the **Bureau of Consumer Financial Protection** in
#'   US about financial products and services.
#'
#' Specifically in this dataset we work with only credit card complaint's for
#'   the 2015 year.
#'
#'
#' @name complaints
#' @docType data
#' @format A data frame with 17301 observations on the following 3 variables.
#'  \describe{
#'     \item{\code{Product}}{a factor with levels \code{Credit card}}
#'     \item{\code{Issue}}{a factor with levels \code{Advertising and marketing}
#'         \code{Application processing delay} \code{APR or interest rate} \code{Arbitration}
#'         \code{Balance transfer} \code{Balance transfer fee} \code{Bankruptcy}
#'         \code{Billing disputes} \code{Billing statement} \code{Cash advance}
#'         \code{Cash advance fee} \code{Closing/Cancelling account}
#'         \code{Convenience checks} \code{Credit card protection / Debt protection}
#'         \code{Credit determination} \code{Credit line increase/decrease}
#'         \code{Customer service / Customer relations} \code{Delinquent account}
#'         \code{Forbearance / Workout plans} \code{Identity theft / Fraud / Embezzlement}
#'         \code{Late fee} \code{Other} \code{Other fee} \code{Overlimit fee} \code{Payoff process}
#'         \code{Privacy} \code{Rewards} \code{Sale of account} \code{Transaction issue}
#'         \code{Unsolicited issuance of credit card}}
#'     \item{\code{Company}}{a factor variable describing the companies available}
#'   }
#' @source \url{http://catalog.data.gov/dataset/consumer-complaint-database}
#' @keywords datasets
#' @usage data(complaints)
#' @examples
#' data(complaints)
#'
NULL

#' Latitude and Longitude Fishnet dataset.
#'
#' This dataset is a subset of Fishnet.
#'
#'
#' @name fishnet
#' @docType data
#' @format A data frame with 4455 observations and 2 variables.
#' \describe{
#'   \item{\code{POINT_X}}{Longitude}
#'   \item{\code{POINT_Y}}{Latitude}
#' }
#' @keywords datasets
#' @usage data(fishnet)
#' @examples
#' data(fishnet)
#'
NULL

#' Landsat TM 5 imagery from 2010 of the Iquitos-Nauta road in the Peruvian
#' Amazon
#'
#' This data set has Binomial data from Landsat TM 5 imagery from 2010 of the
#' Iquitos-Nauta road in the Peruvian Amazon for 7 bands at 69540 locations.
#'
#'
#' @name Landsat
#' @docType data
#' @format A data frame with 69540 observations for 9 columns.
#' @source This dataset is from: Valle D, Baiser B,Woodall CW, Chazdon R
#' (2014). "Decomposing biodiversity data using the Latent Dirichlet Allocation
#' model, a probabilistic multivariate statistical method." Ecology letters,
#' 17(12), 1591-1601.
#' @keywords datasets
#' @usage data(Landsat)
#' @examples
#' data(Landsat)
#'
NULL

#' Species Presence/Absence Data
#'
#' This data set has Presence/Absence predictions for 13 species at 386
#' forested locations. It consists of species, observed presence-absence
#' values, and the probability predictions of three different models.
#'
#'
#' @name presence
#' @docType data
#' @format A data frame with 386 observations for 13 species.  Each cell
#' represents one when the specie is presented zero otherwise.
#' @source This dataset is from: Moisen, G.G., Freeman, E.A., Blackard, J.A.,
#' Frescino, T.S., Zimmerman N.E., Edwards, T.C. Predicting tree species
#' presence and basal area in Utah: A comparison of stochastic gradient
#' boosting, generalized additive models, and tree-based methods. Ecological
#' Modellng, 199 (2006) 176-187.
#' @keywords datasets
#' @usage data(presence)
#' @examples
#' data(presence)
#'
NULL

#' Daily transactions Sp500
#'
#' Daily transactions for 46 firms of the Sp500 Index in 2015.
#'
#'
#' @name sp500
#' @docType data
#' @format A data frame with 249 observations for 46 firms.  Each cell
#' represents one when some transactions ocurred and zero otherwise.
#' @keywords datasets
#' @usage data(sp500)
#' @examples
#' data(sp500)
#'
NULL

#' ID variable Latitude and Longitude Locations for Birds dataset.
#'
#' This dataset is a subset of Birds dataset with ID variable.
#'
#'
#' @name LocationsBirds
#' @docType data
#' @format A data frame with 3080 observations and 3 variables.
#' \describe{
#'   \item{\code{loc.id}}{Location ID}
#'   \item{\code{Latitude}}{Latitude}
#'   \item{\code{Longitude}}{Longitude}
#' }
#' @keywords datasets
#' @usage data(LocationsBirds)
#' @examples
#' data(LocationsBirds)
#'
NULL
