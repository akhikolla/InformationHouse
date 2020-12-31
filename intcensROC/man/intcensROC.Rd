\name{intcensROC}
\alias{intcensROC}
\title{
  Compute the ROC curves for Interval Censored Survival Data
}
\description{
 A method to compute the receiver operating characteristic (ROC) curve for the 
interval censored survival data based on a spline function based constrained 
maximum likelihood estimator. The maximization process of likelihood is 
carried out by generalized gradient projection method.
}

\usage{
  intcensROC(U, V, Marker, Delta, PredictTime, gridNumber = 500)  	
}

\arguments{
  \item{U}{An array contains left end time points of the observation time range
for the interval censored data.}
  \item{V}{An array contains right end time points of the observation time
range for the interval censored data.}
  \item{Marker}{An array contains marker levels for the samples.}
  \item{Delta}{An array of indicator for the censored type, use 1, 2, 3 for
event happened before the left bound time, within the defined time range, and
after. }
  \item{PredictTime}{A scalar indicates the predict time.}
  \item{gridNumber}{A integer for the number of gird for the ROC curve,
the default value is 500.}
}


\references{
    Wu, Yuan; Zhang, Ying. Partially monotone tensor spline estimation 
of the joint distribution function with bivariate current status data. 
Ann. Statist. 40, 2012, 1609-1636 <doi:10.1214/12-AOS1016>
}

\value{
  A \code{dataframe} with two columes
  \item{tp}{A array for true positive rate for different marker levels in the
range of 0 to 1.} 
  \item{fp}{A array for false positive rate for different marker levels in the
range of 0 to 1.}
}


\examples{
## example data
U <- runif(100, min = 0.1, max = 5)
V <- runif(100, min = 0.1, max = 5) + U
Marker <- runif(100, min = 5, max = 10)
Delta  <- sample.int(3, size = 100, replace = TRUE)
pTime <- 4

## compute the ROC curve
res <- intcensROC(U, V, Marker, Delta, pTime, gridNumber = 500)
head(res)
}

