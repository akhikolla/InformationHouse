% $Id: eff.Rd 118 2011-07-22 23:28:53Z Lars $
\name{eff, efficiencies}
\alias{eff}
\alias{efficiencies}
\alias{efficiencies.default}
\alias{eff.sfa}
\alias{eff.Farrell}
\alias{efficiencies.sfa}
\alias{efficiencies.Farrell}


\title{ Calculate efficiencies for Farrell and sfa object }


\description{ Calculate efficiencies for Farrell and sfa object. For a
   sfa there are several types }


\usage{
eff( object, ... )
efficiencies( object, ... )
\method{efficiencies}{default}( object, ... )
\method{efficiencies}{Farrell}(object, type = "Farrell", ...)
\method{eff}{Farrell}(object, type = "Farrell", ...)
\method{efficiencies}{sfa}(object, type = "BC", ...)
\method{eff}{sfa}(object, type = "BC", ...)
}


\arguments{

  \item{object}{ A Farrell object returned from a DEA function like
      \link{dea}, \link{sdea}, or \link{mea} or an sfa object returned
      from the function \link{sfa}.  }

  \item{type}{ The type of efficiencies to be calculated. For a
      Farrell object the possibilities are \dQuote{Farrell} efficiency
      or \dQuote{Shephard} efficiency. For a sfa object the
      possibilities are \dQuote{BC}, \dQuote{Mode}, \dQuote{J}, or
      \dQuote{add}. }

  \item{\dots}{Further arguments ... }

}


\details{

The possible types for class \code{Farrell} (an object returned from
\code{\link{dea}} et al. are \dQuote{Farrell} and \dQuote{Shephard}.

The possible types for class \code{sfa} efficiencies are

\describe{

\item{BC}{ Efficiencies estimated by minimizing the mean square error;
     Eq. (7.21) in Bogetoft and Otto (2011, 219) and Battese and
     Coelli (1988, 392) }

\item{Mode}{ Efficiencies estimates using the conditional mode approach;
     Bogetoft and Otto (2011, 219), Jondrow et al. (1982, 235). }

\item{J}{ Efficiencies estimates using the conditional mean approach
     Jondrow et al. (1982, 235). }

\item{add}{ Efficiency in the additive model, Bogetoft and Otto (2011,
      219) }

}

}


\value{
The efficiencies are returned as an array.
}

\references{
Bogetoft and Otto; \emph{Benchmarking with DEA, SFA, and R}, Springer
2011
}

\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com} }

\note{ For the Farrell object the orientation is determined by the
 calculations that led to the object and can not be changed here.  }


\seealso{
\code{\link{dea}} and \code{\link{sfa}}.
}


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}


\keyword{ efficiency }
