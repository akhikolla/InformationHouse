% $Id: sdea.Rd 229 2020-07-04 13:39:18Z lao $
\name{sdea}
\alias{sdea}


\title{ Super efficiency }


\description{ The method \code{sdea} calculates super-effciency and
   returns the same class of object as \code{\link{dea}}.  }



\usage{
sdea(X, Y, RTS = "vrs", ORIENTATION = "in", DIRECT = NULL, param = NULL,
     TRANSPOSE = FALSE, LP = FALSE, CONTROL = NULL)
}

\arguments{


  \item{X}{Inputs of firms to be evaluated, a K x m matrix
      of observations of K firms with m inputs (firm x input).  In case
  \code{TRANSPOSE=TRUE} the input matrix is transposed to input x
  firm. }

  \item{Y}{Outputs of firms to be evaluated, a K x n matrix
      of observations of K firms with n outputs (firm x input).  In case
  \code{TRANSPOSE=TRUE} the output matrix is transposed to output x
  firm. }

  \item{RTS}{ Text string or a number defining the underlying DEA
      technology / returns to scale assumption; the same values as for
      \code{\link{dea}}.

\tabular{rll}{

0 \tab fdh \tab Free disposability hull, no convexity assumption \cr

1 \tab vrs \tab Variable returns to scale, convexity and free
disposability \cr

2 \tab drs \tab Decreasing returns to scale, convexity, down-scaling
and free disposability \cr

3 \tab crs \tab Constant returns to scale, convexity and free
disposability \cr

4 \tab irs \tab Increasing returns to scale, (up-scaling, but not
down-scaling), convexity and free disposability \cr

5 \tab irs2 \tab Increasing returns to scale (up-scaling, but not
down-scaling), additivity, and free disposability \cr

6 \tab add \tab Additivity (scaling up and down, but only with
integers), and free disposability \cr

7 \tab fdh+ \tab A combination of free disposability and restricted
 or local constant return to scale

}

}

  \item{ORIENTATION}{Input efficiency "in" (1), output efficiency
      "out" (2), and graph efficiency "graph" (3). For use with
      \code{DIRECT}, an additional option is "in-out" (0).}

  \item{DIRECT}{Directional efficiency, \code{DIRECT} is either a
      scalar, an array, or a matrix with non-negative elements.

      If the argument is a scalar, the direction is (1,1,...,1) times
      the scalar; the value of the efficiency depends on the scalar as
      well as on the unit of measurements.

      If the argument an array, this is used for the direction for
      every firm; the length of the array must correspond to the
      number of inputs and/or outputs depending on the
      \code{ORIENTATION}.

      If the argument is a matrix then different directions are used
      for each firm. The dimensions depends on the \code{ORIENTATION},
      the number of firms must correspond to the number of firms in
      \code{X} and \code{Y}.

      \code{DIRECT} must not be used in connection with
      \code{DIRECTION="graph"}.
     }

  \item{param}{Argument is at present only used when
      \code{RTS="fdh+"}, see \code{\link{dea}} for a description.}

  \item{TRANSPOSE}{See the description in \code{\link{dea}}. }

  \item{LP}{Only for debugging, see the description in
      \code{\link{dea}}. }

  \item{CONTROL}{Possible controls to \pkg{lpSolveAPI}, see the
      documentation for that package. For examples of use see the
      function \code{\link{dea}}.}

}


\details{Super-efficiency measures are constructed by avoiding that
   the evaluated firm can help span the technology; ie. if the firm in
   qestuen is a firm on the frontier in a normal dea approach then
   this firm in super efficiency might be outside the technology set.
   }


\value{ The object returned is a Farrell object with the component
   described in \code{\link{dea}}.  The relevant components are

  \item{eff}{The efficiencies.  Note when DIRECT is used then the
      efficencies are not Farrell efficiencies but rather exces values
      in DIRECT units of measurement }

  \item{lambda}{The lambdas, i.e. the weight of the peers, for each
      firm}

  \item{objval}{The objective value as returned from the LP program;
      normally the same as eff. }

  \item{RTS}{The return to scale assumption as in the option
      \code{RTS} in the call}

  \item{ORIENTATION}{The efficiency orientation as in the call}

}



\references{ Peter Bogetoft and Lars Otto; \emph{Benchmarking with
  DEA, SFA, and R}; Springer 2011. Sect. 5.2 page 115

  P Andersen and NC Petersen; \dQuote{A procedure for ranking efficient
  units in data envelopment analysis};  \emph{Management Science} 1993
  39(10):1261--1264}


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}



\note{ Calculation of slacks for super efficiency should be done by
   using the option \code{SLACK=TRUE} in the call of the method
   \code{sdea}.  If the two phases are done in two steps as first a
   call to \code{sdea} and then a call to \code{slacks} the user must
   make sure to set the reference technology to the one corresponding
   to super-efficiency in the call to \code{slack} and this requires a
   loop with calls to \code{slack}.  }


\seealso{
\code{\link{dea}}
}


\examples{
x <- matrix(c(100,200,300,500,100,200,600),ncol=1)
y <- matrix(c(75,100,300,400,25,50,400),ncol=1)
se <- sdea(x,y)
se

# Leave out firm 3 as a determining firm of the technology set
n <- 3
dea.plot.frontier(x[-n], y[-n], txt=(1:dim(x)[1])[-n])
# Plot and label firm 3
points(x[n],y[n],cex=1.25,pch=16)
text(x[n],y[n],n,adj=c(-.75,.75))
}

\keyword{ models }
\keyword{ efficiency }
