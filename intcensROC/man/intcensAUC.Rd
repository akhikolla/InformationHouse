\name{intcensAUC}
\alias{intcensAUC}
\title{
  Compute the AUC for ROC curves for Interval Censored Survival Data
}
\description{
 A method to compute area under the curve (AUC) for the receiver operating characteristic (ROC) curve.}

\usage{
  intcensAUC(ROCdata)  	
}

\arguments{
  \item{ROCdata}{A dataframe from the function \code{intcensROC}}
}

\value{
  A scalar for AUC.
}


\examples{
## example data of ROC curve
U <- runif(100, min = 0.1, max = 5)
V <- runif(100, min = 0.1, max = 5) + U
Marker <- runif(100, min = 5, max = 10)
Delta  <- sample.int(3, size = 100, replace = TRUE)
pTime <- 4

## compute the ROC curve
res <- intcensROC(U, V, Marker, Delta, pTime, gridNumber = 500)
head(res)

##compute the AUC
auc <- intcensAUC(res)
print(auc)
}

