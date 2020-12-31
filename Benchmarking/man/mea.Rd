% $Id: mea.Rd 150 2015-07-06 20:48:17Z B002961 $
\name{mea}
\alias{mea}
\alias{mea.lines}


\title{ MEA multi-directional efficiency analysis }

\description{ Potential improvements PI or multi-directional
efficiency analysis.  The result is an exces value measures by the
direction.

The direction is determined by the direction corresponding to the
minimum input/maximum direction each good can be changed when they are
changed one at a time.  }


\usage{
mea(X, Y, RTS = "vrs", ORIENTATION = "in", XREF = NULL, YREF = NULL, 
    FRONT.IDX = NULL, param=NULL, TRANSPOSE = FALSE, 
    LP = FALSE, CONTROL = NULL, LPK = NULL)
mea.lines(N, X, Y, ORIENTATION="in")  
}



\arguments{

  \item{X}{ K times m matrix with K firms and m inputs as in \code{dea} }

  \item{Y}{ K times n matrix with K firms and n outputs as in \code{dea} }


  \item{RTS}{ Text string or a number defining the underlying DEA
      technology / returns to scale assumption.

\tabular{rll}{
0 \tab fdh \tab Free disposability hull, no convexity assumption \cr
1 \tab vrs \tab Variable returns to scale, convexity and free disposability \cr
2 \tab drs \tab Decreasing returns to scale, convexity, down-scaling and free disposability \cr
3 \tab crs \tab Constant returns to scale, convexity and free disposability \cr
4 \tab irs \tab Increasing returns to scale, (up-scaling, but not down-scaling), convexity and free disposability \cr
6 \tab add \tab Additivity (scaling up and down, but only with integers), and free disposability \cr
7 \tab fdh+ \tab A combination of free disposability and restricted
 or local constant return to scale
}

}

  \item{ORIENTATION}{Input efficiency "in" (1) or output efficiency
      "out" (2), and also the additional option "in-out" (0) for 
      both input and output direction. 
		 
      }

  \item{XREF}{Inputs of the firms determining the technology, defaults
      to \code{X}}

  \item{YREF}{Outputs of the firms determining the technology,
      defaults to \code{Y}}

  \item{FRONT.IDX}{Index for firms determining the technology}

  \item{param}{Possible parameters. At the moment only used for
    RTS="fdh+" to set low and high values for restrictions on lambda;
    see the section details and examples in \code{\link{dea}} for its
    use. Future versions might also use \code{param} for other
    purposes.}

  \item{TRANSPOSE}{ as in \code{dea} }

  \item{LP}{ as in \code{dea} }

  \item{CONTROL}{ as in \code{dea} }

  \item{LPK}{ as in \code{dea} }

  \item{N}{Number of firms where directional lines are to be drawn on
      an all raeady ecisting frontier plot (\link{dea.plot.frontier})}

}


\details{ Details can be found in Bogetoft and Otto (2011, 121--124).

This method is for input directional efficiency only interesting when
there are 2 or more inputs, and for output only when there are 2 or
more outputs.  }



\value{ The results are returned in a Farrell object with the
    following components.

\item{eff}{Exces value in DIRECT units of measurement, this is
    \emph{not} Farrell efficiency }

\item{lambda}{The lambdas, i.e. the weight of the peers, for each firm}

\item{objval}{The objective value as returned from the LP program;
    normally the same as eff}

\item{RTS}{The return to scale assumption as in the option \code{RTS}
     in the call}

\item{ORIENTATION}{The efficiency orientation as in the call}

\item{direct}{ A K times m|n|m+n matrix with directions for each firm:
    the number of columns depends on wether it is input, output og
    in-out orientated.}

\item{TRANSPOSE}{As in the call}

}


\references{ Peter Bogetoft and Lars Otto; \emph{Benchmarking with
   DEA, SFA, and R}; Springer 2011 }


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com} }


\note{ The calculation is done in \code{\link{dea}} after a
   calculation of the dirction that then is used in the argument
   \code{DIRECT}.  The calulation of the direction is done in a series
   LP programs, one for each good in the direction. }


\seealso{ \code{\link{dea}} and the argument \code{DIRECT}.  }


\examples{
X <- matrix(c(2, 2, 5, 10, 10, 3,    12, 8, 5, 4, 6,12), ncol=2)
Y <- matrix(rep(1,dim(X)[1]), ncol=1)

dea.plot.isoquant(X[,1], X[,2],txt=1:dim(X)[1])
mea.lines(c(5,6),X,Y)

me <- mea(X,Y)
me
peers(me)
# MEA potential saving in inputs, exces inputs
eff(me) * me$direct
me$eff *  me$direct

# Compare to traditionally Farrell efficiency
e <- dea(X,Y)
e
peers(e)
# Farrell potential saving in inputs, exces inputs
(1-eff(e)) * X
}

\keyword{ models }
\keyword{ efficiency }
