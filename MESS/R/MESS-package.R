#' Collection of miscellaneous useful and semi-useful functions
#'
#' Collection of miscellaneous useful and semi-useful functions and add-on
#' functions that enhances a number of existing packages and provides In
#' particular in relation to statistical genetics
#'
#' \tabular{ll}{ Package: \tab MESS\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2012-03-29\cr License: \tab GPL-2\cr } % ~~ An overview of
#' how to use the package, including the most important ~~ % ~~ functions ~~
#'
#' @name MESS
#' @aliases MESS-package MESS
#' @docType package
#' @useDynLib MESS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @author Claus Thorn Ekstrøm \email{claus@@rprimer.dk}\cr Maintainer: Claus Thorn Ekstrøm
#' \email{claus@@rprimer.dk}
#' @references Ekstrøm, C. (2011). The R Primer. Chapman & Hall.
#' @import utils stats graphics
#' @keywords package
NULL


#' Danish live births and deaths
#'
#' Monthly live births and deaths in Denmark from January 1901 to March 2013.
#'
#'
#' @name bdstat
#' @docType data
#' @format A data frame with 1356 observations on the following 4 variables.
#' \describe{ \item{year}{a numeric vector giving the month}
#' \item{month}{a numeric vector giving the year}
#' \item{births}{a numeric vector. The number of births for the given
#' month and year} \item{dead}{a numeric vector. The number of deaths
#' for the given month and year} }
#' @source Data were obtained from the StatBank from Danmarks Statistik, see
#' \url{http://www.statbank.dk}.
#' @keywords datasets
#' @examples
#'
#' data(bdstat)
#'
#' plot(bdstat$year + bdstat$month/13, bdstat$birth, type="l")
#'
#' # Create a table of births
#' # Remove year 2013 as it is incomplete
#' btable <- xtabs(births ~ year + month, data=bdstat, subset=(year<2013))
#'
#' # Compute yearly birth frequencies per month
#' btable.freq <- prop.table(btable, margin=1)
#'
NULL





#' Bee data. Number of different types of bees caught.
#'
#' Number of different types of bees caught in plates of different colours.
#' There are four locations and within each location there are three replicates
#' consisting of three plates of the three different colours (yellow, white and
#' blue). Data are collected at 5 different dates over the summer season. Only
#' data from one date available until data has been published.
#'
#'
#' @name bees
#' @docType data
#' @format A data frame with 72 observations on the following 7 variables.
#' \describe{ \item{Locality}{a factor with levels \code{Havreholm}
#' \code{Kragevig} \code{Saltrup} \code{Svaerdborg}. Four different localities
#' in Denmark.} \item{Replicate}{a factor with levels \code{A} \code{B}
#' \code{C}} \item{Color}{a factor with levels \code{Blue} \code{White}
#' \code{Yellow}. Colour of plates} \item{Time}{a factor with levels
#' \code{july1} \code{july14} \code{june17} \code{june3} \code{june6}. Data
#' collected at different dates in summer season. Only one day is present in
#' the current data frame until the full data has been released.}
#' \item{Type}{a factor with levels \code{Bumblebees} \code{Solitary}.
#' Type of bee.} \item{Number}{a numeric vector. The response variable
#' with number of bees catched.} \item{id}{a numeric vector. The id of
#' the clusters (each containg three plates).} }
#' @source Data were kindly provided by Casper Ingerslev Henriksen, Department
#' of Agricultural Sciences, KU-LIFE. Added by Torben Martinussen
#' <tma@@life.ku.dk>
#' @keywords datasets
#' @examples
#'
#' data(bees)
#' model <- glm(Number ~ Locality + Type*Color,
#'              family=poisson, data=bees)
#'
NULL





#' Blood clotting for 158 rats
#'
#' Blood clotting activity (PCA) is measured for 158 Norway rats from two
#' locations just before (baseline) and four days after injection of an
#' anticoagulant (bromadiolone). Normally this would cause reduced blood
#' clotting after 4 days compared to the baseline, but these rats are known to
#' possess anticoagulent resistence to varying extent. The purpose is to relate
#' anticoagulent resistence to gender and location and perhaps weight. Dose of
#' injection is, however, admistered according to weight and gender.
#'
#'
#' @name clotting
#' @docType data
#' @format A data frame with 158 observations on the following 6 variables.
#' \describe{ \item{rat}{a numeric vector} \item{locality}{a
#' factor with levels \code{Loc1} \code{Loc2}} \item{sex}{a factor with
#' levels \code{F} \code{M}} \item{weight}{a numeric vector}
#' \item{PCA0}{a numeric vector with percent blood clotting activity at
#' baseline} \item{PCA4}{a numeric vector with percent blood clotting
#' activity on day 4} }
#' @source Ann-Charlotte Heiberg, project at The Royal Veterinary and
#' Agricultural University, 1999. \cr Added by Ib M. Skovgaard <ims@@life.ku.dk>
#' @keywords datasets
#' @examples
#'
#'  data(clotting)
#'  dim(clotting)
#'  head(clotting)
#'  day0= transform(clotting, day=0, pca=PCA0)
#'  day4= transform(clotting, day=4, pca=PCA4)
#'  day.both= rbind(day0,day4)
#'  m1= lm(pca ~ rat + day*locality + day*sex, data=day.both)
#'  anova(m1)
#'  summary(m1)
#'  m2= lm(pca ~ rat + day, data=day.both)
#'  anova(m2)
#' ## Log transformation suggested.
#' ## Random effect of rat.
#' ## maybe str(clotting) ; plot(clotting) ...
#'
NULL





#' Average yearly summer air temperature for Tasiilaq, Greenland
#'
#' Average yearly summer (June, July, August) air temperature for Tasiilaq,
#' Greenland
#'
#'
#' @name greenland
#' @docType data
#' @format A data frame with 51 observations on the following 2 variables.
#' \describe{ \item{year}{year} \item{airtemp}{average air
#' temperature (degrees Celcius)} }
#' @references Aktuelt Naturvidenskab september 2010. \cr
#' http://aktuelnaturvidenskab.dk/fileadmin/an/nr-4/an4_2010gletscher.pdf
#' @source Data provided by Sebastian Mernild.\cr Originally obtained from
#' http://www.dmi.dk/dmi/index/gronland/vejrarkiv-gl.htm. \cr Added by Claus
#' Ekstrom <ekstrom@@life.ku.dk>
#' @keywords datasets
#' @examples
#'
#' data(greenland)
#' model <- lm(airtemp ~ year, data=greenland)
#' plot(greenland$year, greenland$airtemp, xlab="Year", ylab="Air temperature")
#' abline(model, col="red")
#'
NULL





#' Happiness score and tax rates for 148 countries
#'
#' Dataset on subjective happiness, tax rates, population sizes, continent, and
#' major religion for 148 countries
#'
#'
#' @name happiness
#' @docType data
#' @format A data frame with 148 observations on the following 6 variables.
#' \describe{ \item{country}{a factor with 148 levels that contain the
#' country names} \item{happy}{a numeric vector with the average
#' subject happiness score (on a scale from 0-10)} \item{tax}{a numeric
#' vector showing the tax revenue as percentage of GDP}
#' \item{religion}{a factor with levels \code{Buddhist}
#' \code{Christian} \code{Hindu} \code{Muslim} \code{None} or \code{Other}}
#' \item{continent}{a factor with levels \code{AF}, \code{AS},
#' \code{EU}, \code{NA}, \code{OC}, \code{SA}, corresponding to the continents
#' Africa, Asia, Europe, North America, Ocenaia, South American, respectively}
#' \item{population}{a numeric vector showing the population (in
#' millions)} }
#' @source Data collected by Ellen Ekstroem. \cr Population sizes are from
#' Wikipedia per August 2nd, 2012
#' \url{http://en.wikipedia.org/wiki/List_of_countries_by_population} \cr Major
#' religions are from Wikipedia per August 2nd, 2012
#' \url{http://en.wikipedia.org/wiki/Religions_by_country} \cr Tax rates are
#' from Wikipedia per August 2nd, 2012
#' \url{http://en.wikipedia.org/wiki/List_of_countries_by_tax_revenue_as_percentage_of_GDP}
#' \cr Average happiness scores are from "Veenhoven, R. Average happiness in
#' 148 nations 2000-2009, World Database of Happiness, Erasmus University
#' Rotterdam, The Netherlands". Assessed on August 2nd, 2012 at:
#' \url{http://worlddatabaseofhappiness.eur.nl/hap_nat/findingreports/RankReport_AverageHappiness.php}
#' @keywords datasets
#' @examples
#'
#' data(happiness)
#' with(happiness, symbols(tax, happy, circles=sqrt(population)/8, inches=FALSE, bg=continent))
#'
#' #
#' # Make a prettier image with transparent colors
#' #
#'
#' newcols <- rgb(t(col2rgb(palette())),
#'                alpha=100, maxColorValue=255)
#'
#' with(happiness, symbols(tax, happy, circles=sqrt(population)/8,
#'                 inches=FALSE, bg=newcols[continent],
#'                 xlab="Tax (% of GDP)", ylab="Happiness"))
#'
#' #
#' # Simple analysis
#' #
#' res <- lm(happy ~ religion + population + tax:continent, data=happiness)
#' summary(res)
#'
#'
NULL





#' Gene expression from real-time quantitative PCR
#'
#' Gene expression levels from real-time quantitative polymerase chain reaction
#' (qPCR) experiments on two different plant lines. Each line was used for 7
#' experiments each with 45 cycles.
#'
#'
#' @name qpcr
#' @docType data
#' @format A data frame with 630 observations on the following 4 variables.
#' \tabular{lll}{ \code{flour} \tab numeric \tab Fluorescence level\cr
#' \code{line} \tab factor \tab Plant lines \code{rnt} (mutant) and \code{wt}
#' (wildtype)\cr \code{cycle} \tab numeric \tab Cycle number for the
#' experiment\cr \code{transcript}\tab factor \tab Transcript used for the
#' different runs\cr }
#' @references Morant, M. et al. (2010). Metabolomic, Transcriptional, Hormonal
#' and Signaling Cross-Talk in Superroot2. \emph{Molecular Plant}. 3,
#' p.192--211.
#' @source Data provided by Kirsten Jorgensen <kij@@life.ku.dk>. \cr Added by Claus Ekstrom <ekstrom@@life.ku.dk>
#' @keywords datasets
#' @examples
#'
#' data(qpcr)
#'
#' #
#' # Analyze a single run for the wt line, transcript 1
#' #
#' run1 <- subset(qpcr, transcript==1 & line=="wt")
#'
#' model <- nls(flour ~ fmax/(1+exp(-(cycle-c)/b))+fb,
#'              start=list(c=25, b=1, fmax=100, fb=0), data=run1)
#'
#' print(model)
#'
#' plot(run1$cycle, run1$flour, xlab="Cycle", ylab="Fluorescence")
#' lines(run1$cycle, predict(model))
#'
NULL





#' Perception of points in a swarm
#'
#' Five raters were asked to guess the number of points in a swarm for 10
#' different figures (which - unknown to the raters - were each repeated three
#' times).
#'
#' The raters har approximately 10 seconds to judge each picture, and the
#' thought it was 30 different pictures. Before starting the experiment they
#' were shown 6 (unrelated) pictures and were told the number of points in each
#' of those pictures. The SAND column contains the picture id and the true
#' number of points in the swarm.
#'
#' @name rainman
#' @docType data
#' @format A data frame with 30 observations on the following 6 variables.
#' \describe{ \item{SAND}{The true number of points in the swarm. Each
#' picture is replicated thrice} \item{ME}{Ratings from judge 1}
#' \item{TM}{Ratings from judge 2} \item{AJ}{Ratings from judge
#' 3} \item{BM}{Ratings from judge 4} \item{LO}{Ratings from
#' judge 5} }
#' @source Collected by Claus Ekstrom.
#' @keywords datasets
#' @examples
#'
#' data(rainman)
#' long <- data.frame(stack(rainman[,2:6]), figure=factor(rep(rainman$SAND,5)))
#' figind <- interaction(long$figure,long$ind)
#' # Use a linear random effect model from the
#' # lme4 package if available
#' if(require(lme4)) {
#'   model <- lmer(values ~ (1|ind) + (1|figure) + (1|figind), data=long)
#' }
#'
#' #
#' # Point swarms were generated by the following program
#' #
#' \dontrun{
#' set.seed(2) # Original
#' npoints <- sample(4:30)*4
#' nplots <- 10
#' pdf(file="swarms.pdf", onefile=TRUE)
#'
#' s1 <- sample(npoints[1:nplots])
#' print(s1)
#' for (i in 1:nplots) {
#'   n <- s1[i]
#'   set.seed(n)
#'   x <- runif(n)
#'   y <- runif(n)
#'   plot(x,y, xlim=c(-.15, 1.15), ylim=c(-.15, 1.15), pch=20, axes=FALSE,
#'        xlab="", ylab="")
#' }
#' s1 <- sample(npoints[1:nplots])
#' print(s1)
#' for (i in 1:nplots) {
#'   n <- s1[i]
#'   set.seed(n)
#'   x <- runif(n)
#'   y <- runif(n)
#'   plot(y,x, xlim=c(-.15, 1.15), ylim=c(-.15, 1.15), pch=20, axes=FALSE,
#'        xlab="", ylab="")
#' }
#' s1 <- sample(npoints[1:nplots])
#' print(s1)
#' for (i in 1:nplots) {
#'   n <- s1[i]
#'   set.seed(n)
#'   x <- runif(n)
#'   y <- runif(n)
#'   plot(-x,y, xlim=c(-1.15, .15), ylim=c(-.15, 1.15), pch=20, axes=FALSE,
#'        xlab="", ylab="")
#' }
#' dev.off()
#'}
#'
NULL





#' Danish national soccer players
#'
#' Players on the Danish national soccer team. The dataset consists of all
#' players who have been picked to play on the men's senior A-team, their
#' position, date-of-birth, goals and matches.
#'
#'
#' @name soccer
#' @docType data
#' @format A data frame with 805 observations on the following 5 variables.
#' \describe{ \item{name}{a factor with names of the players}
#' \item{DoB}{a Date. The date-of-birth of the player}
#' \item{position}{a factor with levels \code{Forward} \code{Defender}
#' \code{Midfielder} \code{Goalkeeper}} \item{matches}{a numeric
#' vector. The number of A matches played by the player} \item{goals}{a
#' numeric vector. The number of goals scored by the player in A matches} }
#' @source Data collected from the player database of DBU on March 21st, 2014.
#' See \url{http://www.dbu.dk} for more information.
#' @keywords datasets
#' @examples
#'
#' data(soccer)
#'
#' birthmonth <- as.numeric(format(soccer$DoB, "%m"))
#' birthyear <- as.numeric(format(soccer$DoB, "%Y"))
#'
#'
NULL



#' Estimated life expectancy for Danish newborns
#'
#' The estimated life expectancy for newborn Danes split according to gender.
#'
#' @name lifeexpect
#' @docType data
#' @format  A data frame with 70 observations on the following 3 variables.
#'   \describe{
#'     \item{\code{year}}{a character vector} giving the calendar interval on which the estimation was based.
#'     \item{\code{male}}{a numeric vector} Life expectancy for males (in years).
#'     \item{\code{female}}{a numeric vector} Life expectancy for females (in years)
#'     \item{\code{myear}}{a numeric vector} The midpoint of the year interval
#'   }
#' @source Data collected from Danmarsk Statistik.
#' See \url{https://www.dst.dk/en} for more information.
#' @keywords datasets
#' @examples
#'
#' data(lifeexpect)
#' plot(lifeexpect$myear, lifeexpect$male)
#'
#'
NULL


#' Gene expression data from two-color dye-swap experiment
#'
#' Gene expression levels from two-color dye-swap experiment on 6 microarrays.
#' Arrays 1 and 2 represent the first biological sample (ie, the first dye
#' swap), 3 and 4 the second, and arrays 5 and 6 the third.
#'
#'
#' @name superroot2
#' @docType data
#' @format A data frame with 258000 observations on the following 5 variables.
#' \describe{ \item{color}{a factor with levels \code{green} \code{red}
#' representing the dye used for the gene expression} \item{array}{a
#' factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6}
#' corresponding to the 6 arrays} \item{gene}{a factor with 21500
#' levels representing the genes on the arrays} \item{plant}{a factor
#' with levels \code{rnt} \code{wt} for the two types of plants: runts and wild
#' type} \item{signal}{a numeric vector with the gene expression level
#' (normalized but not log transformed)} }
#' @references Morant, M. et al. (2010). Metabolomic, Transcriptional, Hormonal
#' and Signaling Cross-Talk in Superroot2. \emph{Molecular Plant}. 3,
#' p.192--211.
#' @source Data provided by Soren Bak <bak@@life.ku.dk>. \cr Added by Claus
#' Ekstrom <ekstrom@@sund.ku.dk>
#' @keywords datasets
#' @examples
#'
#' data(superroot2)
#' # Select one gene
#' g1 <- superroot2[superroot2$gene=="AT2G24000.1",]
#' model <- lm(log(signal) ~ plant + color + array, data=g1)
#' summary(model)
#'
NULL


#' Earthquakes in 2015
#'
#' Information on earthquakes worldwide in 2015 with a magnitude greater than 3 on the Richter scale. The variables are just a subset of the variables available at the source
#'
#'
#' @name earthquakes
#' @docType data
#' @format A data frame with 19777 observations on the following 22 variables.
#' \describe{ \item{time}{a factor with time of the earthquake}
#' \item{\code{latitude}}{a numeric vector giving the decimal degrees latitude. Negative values for southern latitudes}
#' \item{\code{longitude}}{a numeric vector giving the decimal degrees longitude. Negative values for western longitudes}
#' \item{\code{depth}}{Depth of the event in kilometers}
#' \item{\code{mag}}{The magnitude for the event}
#' \item{\code{place}}{a factor giving a textual description of named geographic region near to the event. }
#' \item{\code{type}}{a factor with levels \code{earthquake} \code{mining explosion} \code{rock burst}}
#'   }
#' @source \url{http://earthquake.usgs.gov/}
#' @keywords datasets
#' @examples
#'
#' data(earthquakes)
#' with(earthquakes, place[which.max(mag)])
#'
NULL



#' Non-parametric Kruskal Wallis data example
#'
#' Artificial dataset to show that the p-value obtained for the Kruskal Wallis is only valid _after_ the distributional form
#' has been checked to be the same for all groups.
#'
#' @name kwdata
#' @docType data
#' @format An artificial data frame with 18 observations in each of three groups.
#' \describe{ \item{x}{measurements for group 1}
#' \item{y}{measurements for group 2}
#' \item{z}{measurements for group 3}}
#' @source Data example found on the internet
#' @keywords datasets
#' @examples
#'
#' data(kwdata)
#' newdata <- stack(kwdata)
#' kruskal.test(values ~ ind, newdata)
#'
NULL

#' Effect of smoking on self reported health
#'
#' Effect of smoking at 45 years of age on self reported health five years later. Data are on a sample of males from the Glostrup survey.
#'
#'
#'
#' @name smokehealth
#' @docType data
#' @format A table with daily smoking categories for the rows and self reported health five years later as the columns.
#' @source Data example found on the internet but originates from Svend Kreiner
#' @keywords datasets
#' @examples
#'
#' data(smokehealth)
#' m <- smokehealth
#' m[,3] <- m[,3]+ m[,4]
#' m[4,] <- m[4,] + m[5,]
#' m <- m[1:4,1:3]
#' gkgamma(m)
#' chisq.test(m)
#'
NULL

#' Ozone concentration damage to picea spruce
#'
#' Damage scores (ordinal scale) for Picea Sitchesis shoots at two dates, at four temperatures, and 4 ozone Levels
#'
#' @name picea
#' @docType data
#' @format An artificial data frame with 18 observations in each of three groups.
#' \describe{\item{date}{a character vector giving the date}
#' \item{temp}{temperature in degrees Celcius}
#' \item{conc}{Ozone concentration at 4 different levels}
#' \item{damage}{the damage score from 0-4, higher is more damage}
#' \item{count}{The number of occurrences of this group}}
#' @source P.W. Lucas, D.A. Cottam, L.J. Sheppard, B.J. Francis (1988). "Growth
#' Responses and Delayed Winter Hardening in Sitka Spruce Following Summer
#' Exposure to Ozone," New Phytologist, Vol. 108, pp. 495-504.
#' @keywords datasets
#' @examples
#'
#' data(picea)
#'
NULL

#' Ice cream consumption and advertising
#'
#' The impact of advertizing impact, temperature, and price on ice cream consumption
#'
#' @name icecreamads
#' @docType data
#' @format A data frame with 30 observations on the following 4 variables.
#' \describe{\item{Price}{a numeric vector character vector giving the standardized price}
#' \item{Temperature}{temperature in degrees Fahrenheit}
#' \item{Consumption}{a factor with levels \code{1_low} \code{2_medium} \code{3_high}}
#' \item{Advertise}{a factor with levels \code{posters} \code{radio} \code{television}}}
#' @source Unknown origin
#' @keywords datasets
#' @examples
#'
#' data("icecreamads")
#'
NULL


#' Ammonia nitrogen found in river
#'
#' Monthly levels of ammonia nitrogen in a river over two years
#'
#' @name nh4
#' @docType data
#' @format A data frame with 120 observations on the following 3 variables.
#' \describe{\item{nh4}{The ammonia nitrogen levels (mg/l). A value of zero corresponds to a censoring, but it really is censored at <0.01}
#' \item{cens}{A logical vector indicating if the value was censored}
#' \item{year}{The year}}
#' @source Found on the internet and partly simulated
#' @keywords datasets
#' @examples
#'
#' data(nh4)
#'
NULL



#' Flu hospitalization
#'
#' Researchers in a Midwestern county tracked flu cases requiring hospitalization in those residents aged 65
#' and older during a two-month period one winter. They matched each case with 2 controls by sex and age (150 cases,
#' 300 controls). They used medical records to determine whether cases and controls had received a flu vaccine shot
#' and whether they had underlying lung disease. They wanted to know whether flu vaccination prevents hospitalization
#' for flu (severe cases of flu). Underlying lung disease is a potential confounder.
#'
#' @name matched
#' @docType data
#' @format A data frame with 450 observations on the following 4 variables.
#'   \describe{
#'     \item{\code{id}}{a numeric vector}
#'     \item{\code{iscase}}{a factor with levels \code{Control} \code{Case}}
#'     \item{\code{vaccine}}{a factor with levels \code{Not} \code{Vaccinated}}
#'     \item{\code{lung}}{a factor with levels \code{None} \code{Disease}}
#'   }
#' @source Modified from: Stokes, Davis, Koch (2000). "Categorical Data Analysis Using the SAS System," Chapter 10.
#' @keywords datasets
#' @examples
#'
#' data(matched)
#'
NULL
