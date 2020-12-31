% $Id: milkProd.Rd 114 2011-04-10 20:55:50Z Lars $
\name{milkProd}
\alias{milkProd}
\docType{data}

\title{ Data: Milk producers }

\description{
Data colected from Danish  milk producers.
}

\usage{data(milkProd)}

\format{
  A data frame with 108 observations on the following 5 variables.
  \describe{
    \item{\code{farmNo}}{farm number}
    \item{\code{milk}}{Output of milk, kg}
    \item{\code{energy}}{Energy expenses}
    \item{\code{vet}}{Veterinary expenses}
    \item{\code{cows}}{Number of cows}
  }
}

% \details{
% }

\note{ Data as \code{.csv} are loaded by the command \code{data} using
\code{ read.table(..., header = TRUE, sep = ";")} such that this file
is a semicolon separated file and not a comma separated file.  }

\source{
Accounting and buissiness check data
}

\examples{
data(milkProd)
y <- with(milkProd, cbind(milk))
x <- with(milkProd, cbind(energy, vet, cows))
}
\keyword{datasets}
