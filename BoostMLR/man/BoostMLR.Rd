\name{BoostMLR}
\alias{BoostMLR}
\title{
Boosting for Multivariate Longitudinal Response
}
\description{
Function jointly models the multiple longitudinal responses (referred to as multivariate longitudinal response) and multiple covariates and time from a longitudinal study using gradient boosting approach (Pande et al., 2020). 
Covariates can be time-varying or time-invariant.
Special cases include modeling of univariate longitudinal
response from a longitudinal study, and univariate or multivariate 
response from a cross-sectional study.
In all cases, responses are assumed to be continuous.
The estimated coefficient can be a function of time (referred to as time-varying coefficient in case of a longitudinal study) or a function of pre-specified covariate (in case of a longitudinal or a cross-sectional study) or fixed.
}
\usage{
  BoostMLR(x,
           tm,
           id,
           y,
           Time_Varying = TRUE,
           BS_Time = TRUE,
           nknots_t = 10,
           d_t = 3,
           All_RawX = TRUE,
           RawX_Names,
           nknots_x = 7,
           d_x = 3,
           M = 200,
           nu = 0.05,
           Mod_Grad = TRUE,
           Shrink = FALSE,
           VarFlag = TRUE,
           lower_perc = 0.25,
           upper_perc = 0.75,
           NLambda = 100,
           Verbose = TRUE,
           ...)
}
\arguments{
\item{x}{Data frame (or matrix) containing x-values (covariates).
The number of rows should match with number of rows of response \code{y}.
Covariates can be time-varying or time-invariant.
Missing values are allowed, and are ignored during estimation.
If unspecified, model will be fitted with time alone
(applicable in the situation when the interest is to
obtain an estimated mean response trajectory over time without the influence of any covariates).}
\item{tm}{Vector of time values, one entry for each row of the response \code{y}.
In case of a longitudinal study, the estimated coefficient will be a function
of \code{tm} when \code{Time_Varying} = TRUE.
If unspecified, data is assumed to be generated from a
cross-sectional study, and the relationship between
\code{y} and \code{x} can be obtained.
In case of a longitudinal or cross-sectional study, coefficient can be
a function of covariate \code{z} which is not a part of \code{x} by
using \code{z} in place of \code{tm} or it can be fixed when \code{Time_Varying} = FALSE.}
\item{id}{Vector of subject identifier with same length as the number of rows of \code{y}.
If \code{id} is unspecified along with \code{tm}, data is assumed to be generated from a cross-sectional study, and the relationship between
\code{y} and \code{x} can be obtained.}
\item{y}{Data frame (or matrix) containing the y-values (response) in case of multivariate
response or a vector of y-values in case of univariate response.}
\item{Time_Varying}{Time-varying coefficient model or a fixed coefficient model?}
\item{BS_Time}{If \code{tm} is specified, should \code{tm} is mapped using B-spline or use original scale of \code{tm}? Default is TRUE, which allows mapping of \code{tm}
using B-spline.}
\item{nknots_t}{If \command{BS_Time} = TRUE, specify number of knots for B-spline of \code{tm}.}
\item{d_t}{If \command{BS_Time} = TRUE, specify degree of polynomial for B-spline of \code{tm}.}
\item{All_RawX}{Use original scale of \code{x} or map each covariate using
B-spline? Default is TRUE, which means original scale of \code{x} is used;
if FALSE, covariates measured on continuous scale will be mapped using B-splines.}
\item{RawX_Names}{If \command{All_RawX} = FALSE, specify names of the covariates,
measured on a continuous scale, that should be used as it is without mapping
using B-spline.
Note that, even if \command{All_RawX} = FALSE, covariates not measured on a
continuous scale, such as binary, nominal, and ordinal covariates will be
used without mapping.}
\item{nknots_x}{Specify number of knots for B-spline of \code{x}. This can be a vector of length equal to the number of covariates or a scalar. 
If scalar, same value will be used for all covariates.}
\item{d_x}{Specify degree of polynomial for B-spline of \code{x}. This can be a vector of length equal to the number of covariates or a scalar. 
If scalar, same value will be used for all covariates.}
\item{M}{Number of boosting iterations.}
\item{nu}{Boosting regularization parameter.  A value from the interval (0,1].}
\item{Mod_Grad}{Use a modified gradient? Modified gradient is a special type
of gradient that is independent of the correlation coefficient.
Pande A. (2017) observed that prediction performance increases
under modified gradient.}
\item{Shrink}{Allow estimated coefficient to shrink to zero using
L1 penalization?}
\item{VarFlag}{Estimate the variance (scale parameter) and correlation parameter
for each \code{y}? Applicable for a longitudinal study. If \code{VarFlag} = FALSE,
a fixed value of scale parameter = 1 and correlation parameter = 0 is used.}
\item{lower_perc}{Lower percentile value is used to determine the lower cut-off for the distribution of parameter estimate. Applicable when \code{Shrink} = TRUE.
Refer to Pande et al. (2020) for details.}
\item{upper_perc}{Upper percentile value is used to determine the upper cut-off for the distribution of parameter estimate. Applicable when \code{Shrink} = TRUE.
Refer to Pande et al. (2020) for details.}
\item{NLambda}{Number of replications for generating distribution of parameter estimates. Applicable when \code{Shrink} = TRUE. Refer to Pande et al. (2020) for details.}
\item{Verbose}{Print the current stage of boosting iteration?}
\item{...}{Further arguments passed to or from other methods.}
%\item{dt_Add}{Use \code{dt_Add} to specify additional time-varying covariates
%for which time of measurements are different than time of measurements of the response. %\code{dt_Add} must be a list of length equal to number of time-varying
%covariates. Each element of list must have 3 columns in the order of \code{id},
%\code{tm} and \code{x}.}
%\item{Lambda_Ridge}{A constant is used for the estimation of coefficient to
%prevent the instability.}
}
\details{
This is a non-parametric approach for joint modeling of a multivariate
longitudinal response, which is based on marginal model.
Estimation is performed using gradient boosting, a generic form
of boosting (Friedman J. H., 2001). Our boosting approach is closely related to
component-wise L2 boosting with L1 penalization.
Approach can handle high dimensionalilty of covariate and response
when some of the covariates and responses are pure noise.

Approach is designed to identify covariates that affect responses
differently as different time intervals. This idea is helpful
to dissect an overall effect of covariate into different time
intervals. For example, some covariates affect response at
the beginning of the follow-up whereas others at a later stage.

Shrinking allows for early termination of boosting to prevent overfitting. Also, it provides a parsimonious model by shrinking coefficient for non-informative covariate-response pair to zero.}

\value{
\item{x}{Matrix containing x-values.}
\item{id}{Vector of subject identifier.}
\item{tm}{Vector of time values.}
\item{y}{Matrix containing y-values.}
\item{UseRaw}{Logical vector indicating indexes of covariates which are used as it is
without B-spline mapping.}
\item{x_Names}{Variable names of \code{x}.}
\item{y_Names}{Variable names of \code{y}.}
\item{M}{Number of boosting iterations. If boosting terminates before
a pre-specified \code{M}, this indicates the last boosting iteration 
before termination.}
\item{nu}{Regularization parameter.}
\item{Tm_Beta}{An estimate of the parameter beta. This consist of a list of
length equal to the number of multivariate response (denoted by L). If \code{Time_Varying} = TRUE, each element from the list represents a matrix with number of columns equal to the number of covariates and the number of rows equal to the length of \code{tm}. Each column of the matrix represents an estimate of time-varying coefficient for the given covariate. If \code{Time_Varying} = FALSE, in place of estimate of time-varying coefficient, a fixed estimate is provided similar to the estimate from a parametric model. 
The result is provided for covariates who are treated as it is (i.e., the original scale); for covariates who are mapped using B-spline, the estimates are difficult to interprete and therefore the output is NA.}
\item{mu}{Estimate of the conditional expectation of \code{y} corresponding to the M'th boosting iterations.}
\item{Error_Rate}{Training error rate for each response across 
the boosting iterations.}
\item{Variable_Select}{Indexes of important covariates that get picked-up across time and across boosting iterations. Result is shown as a matrix with M rows and H 
(number of overlapping time intervals) columns, where each element represents
index of covariate.}
\item{Response_Select}{Indexes of important responses that get picked-up across time and across boosting iterations. Result is shown as a matrix with M rows and H
columns, where each element represents index of response variable.}
\item{VarFlag}{Whether the variance (scale parameter) and correlation are estimated?}
\item{Time_Varying}{Whether estimates are time-varying or fixed?}
\item{Phi}{Matrix, having dimension M by L, representing an estimate of variance (scale parameter) for each response across the boosting iterations.}
\item{Rho}{Matrix, having dimension M by L, represent an estimate of correlation for each response across the boosting iterations.}
\item{Lambda_List}{Estimate of the lambda (the L1 penaulty parameter) for each boosting iterations. Useful for internal calculation.}
\item{Grow_Object}{Useful for internal calculation.}
%\item{Time_Add_New}{If \code{dt_Add} is specified, this corresponds to time measurements %from \code{dt_Add}, aligned with \code{tm}.}
}
\author{
Amol Pande and Hemant Ishwaran
}

\references{
  Pande A., Ishwaran H., Blackstone E.H. (2020).  Boosting for multivariate
  longitudinal response.

  Pande A., Li L., Rajeswaran J., Ehrlinger J., Kogalur U.B.,
  Blackstone E.H., Ishwaran H. (2017).  Boosted multivariate trees for
  longitudinal data, \emph{Machine Learning}, 106(2): 277--305.
  
  Pande A. (2017).  \emph{Boosting for longitudinal data}.  Ph.D. Dissertation,
  Miller School of Medicine, University of Miami.
  
  Friedman J.H. (2001). Greedy function approximation: a gradient
  boosting machine, \emph{Ann. of Statist.}, 5:1189-1232.
  
}

\seealso{
  \command{\link{updateBoostMLR}},
  \command{\link{predictBoostMLR}},
  \command{\link{simLong}}
}
\examples{

\donttest{
##-----------------------------------------------------------------
## Multivariate Longitudinal Response
##-----------------------------------------------------------------

# Simulate data involves 3 response and 4 covariates

dta <- simLong(n = 100, N = 5, rho =.80, model = 1, q_x = 0, 
                                  q_y = 0,type = "corCompSym")$dtaL

# Boosting call: Raw values of covariates, B-spline for time, 
# no shrinkage, no estimate of rho and phi

boost.grow <- BoostMLR(x = dta$features, tm = dta$time, id = dta$id, 
                            y = dta$y, M = 100, VarFlag = FALSE)

# Plot training error
plotBoostMLR(boost.grow$Error_Rate,xlab = "m",ylab = "Training Error")


##-----------------------------------------------------------------
## Laboratory data
##-----------------------------------------------------------------

data(Laboratory_Data, package = "BoostMLR")

Var_Names <- colnames(Laboratory_Data)
x_Names <- setdiff(Var_Names, c("id","time","tbili_po","creat_po"))

dta_id <- Laboratory_Data[,"id"]
dta_time <- Laboratory_Data[,"time"]
dta_x <- Laboratory_Data[,x_Names]
dta_y <- Laboratory_Data[,c("tbili_po","creat_po")]

boost.grow <- BoostMLR(x = dta_x,tm = dta_time,id = dta_id,y = dta_y,
                           Time_Varying = TRUE,BS_Time = TRUE,
                           All_RawX = TRUE,M = 10, VarFlag = TRUE)

##-----------------------------------------------------------------
## Univariate Longitudinal Response
##-----------------------------------------------------------------

# Simulate data involves 1 response and 4 covariates

dta <- simLong(n = 100, N = 5, rho =.80, model = 2, q_x = 0, 
                                  q_y = 0,type = "corCompSym")$dtaL

# Boosting call: B-spline for time and covariates, shrinkage, 
# estimate of rho and phi 

boost.grow <- BoostMLR(x = dta$features, tm = dta$time, id = dta$id, 
                          y = dta$y, M = 100, BS_Time = TRUE,
                          All_RawX = FALSE, Shrink = TRUE,VarFlag = TRUE)

# Plot training error
plotBoostMLR(boost.grow$Error_Rate,xlab = "m",ylab = "Training Error")

# Plot phi
plotBoostMLR(boost.grow$Phi,xlab = "m",ylab = "phi")

# Plot rho
plotBoostMLR(boost.grow$Rho,xlab = "m",ylab = "rho")


##-----------------------------------------------------------------
## Multivariate Longitudinal Response
##-----------------------------------------------------------------

# Simulate data involves 3 response and 4 covariates

dta <- simLong(n = 100, N = 5, rho =.80, model = 1, q_x = 0, 
                                  q_y = 0,type = "corCompSym")$dtaL

# Boosting call: Raw values of covariates, fixed parameter estimates
# instead of time varying, no shrinkage, no estimate of rho and phi

boost.grow <- BoostMLR(x = dta$features, tm = dta$time, id = dta$id, 
                       y = dta$y, M = 100,Time_Varying = FALSE,VarFlag = FALSE)

# Print parameter estimates
boost.grow$Tm_Beta

##-----------------------------------------------------------------
## Multivariate Response from Cross-sectional Data: Estimated 
## coefficient as a function of covariate
##-----------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {
data("BostonHousing")

x <- BostonHousing[,c(1:7,9:12)]
tm <- BostonHousing[,8]
id <- 1:nrow(BostonHousing)
y <- BostonHousing[,13:14]

# Boosting call: Raw values of covariates, B-spline for covariate "dis", 
# no shrinkage

boost.grow <- BoostMLR(x = x, tm = tm, id = id, y = y, M = 100,VarFlag = FALSE)

# Plot training error
plotBoostMLR(boost.grow$Error_Rate,xlab = "m",ylab = "Training Error",
                                              legend_fraction_x = 0.2)
}
##-----------------------------------------------------------------
## Univariate Response from Cross-sectional Data: Fixed estimated 
## coefficient
##-----------------------------------------------------------------

if (library("mlbench", logical.return = TRUE)) {
library(mlbench)
data("BostonHousing")

x <- BostonHousing[,1:13]
y <- BostonHousing[,14]

# Boosting call: Raw values of covariates

boost.grow <- BoostMLR(x = x, y = y, M = 100)

# Plot training error
plotBoostMLR(boost.grow$Error_Rate,xlab = "m",ylab = "Training Error",
                                               legend_fraction_x = 0.2)
}
}
}
\keyword{boosting}
