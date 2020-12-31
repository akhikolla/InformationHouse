% $Id: projekt.Rd 80 2010-10-20 08:58:43Z Lars $
\name{projekt}
\alias{projekt}
\docType{data}
\title{
Data:  Milk producers
}

\description{
Accounting and production data for 101 milk producing farmers.
}

\usage{data(projekt)}
\format{
  A data frame with 101 observations on the following 14 variables.
  \describe{
    \item{\code{numb}}{Serial number for the milk producer}
    \item{\code{cows}}{Number of cows}
    \item{\code{vet}}{Veterinary expences}
    \item{\code{unitCost}}{Unit cost, variable cost}
    \item{\code{capCost}}{Capacity cost}
    \item{\code{fixedCost}}{Fixed cost}
    \item{\code{milkPerCow}}{Milk per cow, kg}
    \item{\code{quota}}{Milk quota}
    \item{\code{fatPct}}{Fat percent in milk}
    \item{\code{protPct}}{Protein percent in milk}
    \item{\code{cellCount}}{Cell count for milk}
    \item{\code{race}}{Race for cows, a factor with levels \code{jersey}, 
             \code{large}, and \code{mixed}}
    \item{\code{type}}{Type of production, conventional or organic, 
             a factor with levels \code{conv} \code{orga}}
    \item{\code{age}}{Age of the farmer}
  }
}

\details{
Data is a mix of accounting data and production controls.
}

\note{ Data as \code{.csv} are loaded by the command \code{data} using
\code{ read.table(..., header = TRUE, sep = ";")} such that this file
is a semicolon separated file and not a comma separated file.  }


\source{
Collected from farmers.
}


\examples{
data(projekt)
## maybe str(projekt) ; plot(projekt) ...
}
\keyword{datasets}
