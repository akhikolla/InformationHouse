\name{ggas}
\alias{ggas}
\docType{data}
\title{Data from Gastrointestinal Tumor Study Group (1982)}
\description{
The Gastrointestinal Tumor Study Group (1982) compared chemotherapy with combined chemotherapy and radiation therapy, in the treatment of locally unresectable gastric cancer. Each treatment arm had 45 patients, with  two observations of the chemotherapy group and six of the combination group censored. The plot of Kaplan-Meier survival curves show that the two curves cross at around 1000 days.
}
\usage{data(ggas)}
\format{
  The data has the following information for the ninety subjects:
  \describe{
    \item{\code{time}}{the vector for the pooled lifetimes of the two groups}
    \item{\code{event}}{the numeric vector of the right-censoring indicator (event = 1, censored = 0)}
    \item{\code{group}}{the numeric vector of the group indicator (treatment = 1, control = 0)}
  }
}
\references{
Gastrointestinal Tumor Study Group (1982). A comparison of combination chemotherapy and combined modality therapy for locally advanced gastric carcinoma. Cancer.
}
\examples{
library(ClinicalTrialSummary)
data(ggas)
}
\seealso{\code{\link{ypsummary}}}
\keyword{Data}

