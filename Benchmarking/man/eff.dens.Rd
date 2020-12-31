% $Id: eff.dens.Rd 131 2014-06-24 17:00:12Z b002961 $
\name{eff.dens}
\alias{eff.dens}
\alias{eff.dens.plot}

\title{
Estimate and plot density of efficiencies
}

\description{ A method to estimate and plot kernel estimate of
(Farrell) efficiencies taken into consideration that efficiencies are
bounded either above (input direction) or below (output direction).  }

\usage{
eff.dens(eff, bw = "nrd0")

eff.dens.plot(obj, bw = "nrd0", ..., xlim, ylim, xlab, ylab)
}


\arguments{

\item{eff}{ Either a list of (Farrell) efficiencies or a Farrell
object returned from the method \code{\link{dea}}.  }

\item{bw}{ Bandwith, look at the documentation of \code{density} for
an explanation.  }

\item{obj}{Either an array of efficiencies or a list returned from
\code{eff.dens}.}

\item{\dots}{ Further arguments to the \code{plot} method like line
type and line width.  }

\item{xlim}{ Range on the x-axis; usualy not needed, just use the
defaults. }

\item{ylim}{ Range on the x-axis; usualy not needed, just use the
defaults. }

\item{xlab}{ Label for the x-axis. }

\item{ylab}{ Label for the y-axis. }

}


\details{ The calculation is based on a reflection method (Silverman
1986, 30) using the default window kernel and defult bandwidth (window
width) in the method \code{density}.  

The method \code{eff.dens.plot} plot the density directly, and
\code{eff.dens} just estimate the numerical density, and the result
can then either be plotted by \code{plot}, corresponds to
\code{eff.dens.plot}, or by lines as an overlay on an existing plot.
}


\value{ The return from \code{eff.dens} is a list \code{list(x,y)}
with efficiencies and the corresponding density values.  }

\references{ B.W. Silverman (1986), \emph{Density Estimation for Statistics
and Data Analysis}, Chapman and Hall, London.  }


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\note{ The input efficiency is also bounded below by 0, but for normal
firms an efficiency at 0 will not happen, i.e. the boundary is not
effective, and therefore this boundary is not taken into
consideration.  }


\examples{
e <- 1 - rnorm(100)
e[e>1] <- 1
e <- e[e>0]
eff.dens.plot(e)

hist(e, breaks=15, freq=FALSE, xlab="Efficiency", main="")
den <- eff.dens(e)
lines(den,lw=2)
}

% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
