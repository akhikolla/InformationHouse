#' Instruction
#' 
#' A data frame on school instruction results.
#' 
#' 
#' @name instruction
#' @docType data
#' @format A data frame with 1190 observations on the following 13 variables.
#' \describe{
#'    \item{X}{a numeric vector}
#'    \item{girl}{a numeric vector}
#'    \item{minority}{a numeric vector}
#'    \item{mathkind}{a numeric vector}
#'    \item{mathgain}{a numeric vector}
#'    \item{ses}{a numeric vector}
#'    \item{yearstea}{a numeric vector}
#'    \item{mathknow}{a numeric vector}
#'    \item{housepov}{a numeric vector}
#'    \item{mathprep}{a numeric vector}
#'    \item{classid}{a numeric vector identifying the class within school}
#'    \item{schoolid}{a numeric vector identifying the school}
#'    \item{childid}{a numeric vector}
#' }
#' @source West, B., Welch, K. B., & Galecki, A. T. (2006). Linear mixed
#' models: a practical guide using statistical software. Chapman & Hall/CRC.
#' @keywords datasets
#' @examples
#' 
#' # The following code takes a few minutes to run.
#' # In the interest of saving CRAN's example testing time,
#' # it has been commented out. If you want to use it,
#' # just uncomment and run.
#' 
#' # data(instruction)
#' # attach(instruction)
#' 
#' # data = data.frame(
#' #   y = mathgain,
#' #   mathkind = mathkind, 
#' #   girl = girl,
#' #   minority = minority,
#' #   ses = ses, 
#' #   school = factor(schoolid), 
#' #   section = factor(classid))
#' 
#' 
#' # fit.rlme = rlme(y ~ 1 + mathkind + girl + minority + ses + (1 | school) + (1 | school:section),
#' #  data = data,
#' #  method = "gr")
#'   
#' # summary(fit.rlme)
NULL

#' rlme
#' 
#' An R package for rank-based robust estimation and prediction in random
#' effects nested models
#' 
#' \tabular{ll}{ Package: \tab rlme\cr Type: \tab Package\cr Version: \tab
#' 0.2\cr Date: \tab 2013-07-07\cr License: \tab GPL (>= 2)\cr }
#' 
#' @name rlme-package
#' @docType package
#' @author Yusuf Bilgic \email{bilgic@@geneseo.edu}, Herb Susmann
#' \email{hps1@@geneseo.edu} and Joseph McKean \email{joemckean@@yahoo.com}
#' 
#' Maintainer: Yusuf Bilgic \email{bilgic@@geneseo.edu} or
#' \email{yusuf.k.bilgic@@gmail.com}
#' @seealso \code{\link{rlme}}
#' @keywords models package
#' @import stats
#' @import graphics
#' @import Rcpp
#' @useDynLib rlme
#' @examples
#' 
#' 
#' library(rlme)
#' data(schools)
#' formula = y ~ 1 + sex + age + (1 | region) + (1 | region:school)
#' rlme.fit = rlme(formula, schools)
#' summary(rlme.fit)
#' 
NULL

#' PISA Literacy Data
#' 
#' The data in Program for International Assessment (PISA) on academic
#' proficiency in schools around the world.
#' 
#' 
#' @name schools
#' @docType data
#' @format A data frame with 334 observations on the following 6 variables.
#' \describe{
#'    \item{y}{a numeric vector indicating student literacy}
#'    \item{socio}{a numeric vector}
#'    \item{sex}{a numeric vector}
#'    \item{age}{a numeric vector}
#'    \item{region}{a numeric vector indicating four regions}
#'    \item{school}{a numeric vector indicating the schools within region}
#' }
#' @references OECD (2010). PISA 2009 Results. http://www.oecd.org/
#' @keywords datasets
#' @examples
#' 
#' #
#' # The example takes a few seconds to run, so in order to 
#' # save CRAN's testing time it has been commented out. 
#' # To run, simply uncomment and execute.
#' #
#' 
#' # data(schools)
#' # rlme.fit = rlme(y ~ 1 + sex + age + (1 | region) + (1 | region:school), 
#' #	schools, method="gr")
#' # summary(rlme.fit)
NULL