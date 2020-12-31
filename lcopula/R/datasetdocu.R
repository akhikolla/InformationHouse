#' Ranks of preprocessed monthly Danube river flow measurements
#' 
#' The \code{danube} dataset contains ranks of base flow observations from the Global River Discharge project of the Oak Ridge National Laboratory Distributed Active Archive Center (ORNL DAAC), a NASA data center. The measurements are monthly average flow rate for two stations situated at Scharding (Austria) on the Inn river and at Nagymaros (Hungary) on the Danube.
#'  The data have been pre-processed to remove any time trend. Specifically, Bacigal et al. (2011) extracted the raw data, and obtain the fast Fourier transformed centered observations. The negative spectrum is retained and a linear time series model with 12 seasonal components is fitted. Residuals are then extracted and AR model fitted to the series, the selection being done based on the AIC criterion with imposed maximum order of 3 and the number of autoregressive components may differ for each series.
#'  
#' @format 
#' This data frame contains the following columns:
#' \itemize{
#'      \item{inn}{A numeric vector containing the rank of prewhitened level observations of the Inn river at Nagyramos.}
#'      \item{donau}{A numeric vector containing the rank of prewhitened level observations of the Donau river at Scharding.}
#' }
#' @name danube
#' @docType data
#' @keywords datasets
#' @source   Vorosmarty, C.J., Fekete, B.M., Tucker, B.A. (1998).   Global River Discharge, 1807--1991, V. 1.1(RivDIS). Data set.   Available from \code{http://doi.org/10.3334/ORNLDAAC/199}
#' @references
#' Bacigal, T., Jagr, V., Mesiar, R. (2011)  Non-exchangeable random variables, Archimax copulas and their fitting to real data.   \emph{Kybernetika}, \bold{47}(4), 519--531.  
#' 
#' Genest, C. and Neslehova, J. G. (2013)  Assessing and Modeling Asymmetry in Bivariate Continuous Data.  In P. Jaworski, F. Durante, and W. K. Hardle (Eds.),  Copulae in Mathematical and Quantitative Finance, Lecture Notes in Statistics, 91--114, Springer: Berlin Heidelberg. 
#'  
NULL

#' Daily air pollutant measures in Leeds
#' 
#'  The \code{airquality} data frame consists of measurements of pollutants in the city of Leeds. In the saying of Boldi and Davison (2007): 
#'  "The dataset comprise daily series of monitoring measurements of ozone levels (O3), nitrogen dioxide (NO2), nitrogen oxide (NO) and particulate matter (PM10), in the city centre of Leeds, UK, over 1994--1998. Levels of the gases are measured in parts per billion, and those of PM 10 in micrograms per cubic metre. "
#'
#' Downloaded from \url{http://www.airquality.co.uk}. (c) Crown 2015 copyright Defra via uk-air.defra.gov.uk, licenced under the Open Government Licence (OGL).
#' 
#' @references  
#' Heffernan, J., Tawn, J. (2004)  A conditional approach for multivariate extreme values (with discussion).   \emph{J. R. Stat. Soc., Ser. B Stat. Methodol.} \bold{66}(3), 497--546.
#' 
#' Sabourin, A. , Naveau, P., Fougeres, A.-L. (2013) Bayesian Model averaging for Multivariate extremes. \emph{Extremes}, \bold{16}(3), 325--350.
#' 
#'  Boldi, M.O., Davison, A.C. (2007)  A mixture model for multivariate extremes. \emph{J. R. Stat. Soc., Ser. B Stat. Methodol.} \bold{69}(2), 217--229.
#'  
#'  Cooley, D., Davis, R., Naveau, P. (2010)  The pairwise beta distribution: A flexible parametric multivariate model for extremes.   \emph{J. Multivar. Anal.} \bold{101}(9), 2103--2117.
#' @keywords datasets
#' @name airquality
#' @docType data
NULL

#' Women daily nutrient intake
#' 
#'  The \code{nutrient} data frame consists of quintuples consisting of four day measurements for intake of calcium, iron, protein, vitamin A and C from women aged 25 to 50 in the United States as part of the "Continuing Survey of Food Intakes of Individuals" program. The processed data has
#'  737 measurements from a cohort study of the United States Department of Agriculture (USDA) and it is available online from the University of Pennsylvania repository. 
#'
#' @format
#' This data frame contains the following columns:
#' \describe{
#'  \item{id}{A numeric vector containing the identification of the participant.}
#'  \item{calcium}{A numeric vector containing measured calcium.}
#'  \item{iron}{A numeric vector containing iron measurements.}
#'  \item{protein}{A numeric vector containing protein measurements.}
#'  \item{vitamin.a}{A numeric vector containing vitamin A measurements.}
#'  \item{vitamin.c}{A numeric vector containing vitamin C measurements.}
#'    }
#' @source 
#'  The survey data was processed by Dr. Andrew Wiesner for a course on multivariate statistics at The Pennsylvania State University.
#' @references
#' Genest, C., Neslehova, J., Quessy, J.-F. (2012)  Tests of symmetry for bivariate copulas.   \emph{Annals of the Institute of Statistical Mathematics}, \bold{64}(4), 811--834.
#' 
#'  Genest, C. and Neslehova, J. G. (2013)  Assessing and Modeling Asymmetry in Bivariate Continuous Data.   In P. Jaworski, F. Durante, & W. K. Hardle (Eds.),   Copulae in Mathematical and Quantitative Finance, Lecture Notes in Statistics, 91--114, Springer: Berlin Heidelberg. 
#' 
#' Li, B., Genton, M. G. (2013) Nonparametric identification of copula structures.  \emph{Journal of the American Statistical Association}, \bold{108}(502), 666--675.
#' @keywords datasets
#' @docType data
#' @name nutrient
NULL

