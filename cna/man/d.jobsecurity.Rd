\name{d.jobsecurity}

\docType{data}

\alias{d.jobsecurity}

\title{Job security regulations in western democracies}

\description{
This dataset is from Emmenegger (2011), who analyzes the determinants of 
high job security regulations in Western democracies using fsQCA. 
}

\usage{d.jobsecurity}

\format{
The data frame contains 19 rows (cases) and the following 7 columns featuring fuzzy-set factors:

\tabular{rllrll}{
[ , 1] \tab \bold{S}   \tab statism                  \tab ("1" high, "0" not high)\cr
[ , 2] \tab \bold{C}   \tab non-market coordination  \tab ("1" high, "0" not high)\cr
[ , 3] \tab \bold{L}   \tab labour movement strength \tab ("1" high, "0" not high)\cr
[ , 4] \tab \bold{R}   \tab Catholicism              \tab ("1" high, "0" not high)\cr
[ , 5] \tab \bold{P}   \tab religious party strength \tab ("1" high, "0" not high)\cr
[ , 6] \tab \bold{V}   \tab institutional veto points \tab   ("1" many, "0" not many)\cr
[ , 7] \tab \bold{JSR} \tab job security regulations \tab  ("1" high, "0" not high)
}
}

\note{
The row names are the official International Organization for Standardization 
(ISO) country code elements as specified in ISO 3166-1-alpha-2.
}

\section{Contributors}{
Thiem, Alrik: collection, documentation
}

\source{
Emmenegger, Patrick. 2011. \dQuote{Job Security Regulations in Western Democracies: 
A Fuzzy Set Analysis.} \emph{European Journal of Political Research} 50(3):336-64.
}

\keyword{datasets}
